#!/usr/bin/env python3
"""
fixfasta.py — robust auto-orientation of fungal ITS reads

Version 1.1 - June 30, 2025

By Alan Rockefeller 

Useful for building phylogenetic trees from data downloaded from Mycomap - often there are a few reversed sequences

- Accepts FASTA files on the command line (or stdin if none given)
- Repairs stray '>' symbols inside sequence lines
- Uses three conserved ITS motifs (ITS1-F, 5.8 S core, ITS4)
  with IUPAC-aware fuzzy matching (≤ 4 mismatches over 20 bp)
- Orientation logic:
    1. orientation with *more* distinct motif hits wins
    2. if tied, one with *fewer total mismatches* wins
    3. if still tied, one whose *best* hit starts earlier wins
    4. if score difference < TIE_EPS → "uncertain"
- Wraps output FASTA to 80 nt; all diagnostics go to stderr
- Reports which sequences were reversed (use -q to suppress)

Usage: 
    cat fast1.fas fast2.fas | fixfasta.py > fixed_fasta.fas
    fixfasta.py input.fasta -o output.fasta --stats
    fixfasta.py *.fasta --dry-run --verbose
    fixfasta.py input.fasta -q > output.fasta  # quiet: no warnings or reports
"""

from __future__ import annotations
import sys
import re
import fileinput
import argparse
import logging
from dataclasses import dataclass
from textwrap import wrap
from typing import List, Tuple, Optional, Iterator, Dict
from collections import defaultdict

try:
    import edlib
except ImportError:
    print("Error: edlib library is required. Install with: pip install edlib", file=sys.stderr)
    sys.exit(1)

# Constants
FWD_MOTIFS = [
    "TCCGTAGGTGAACCTGCGG",    # ITS1-F  (18 S end)
    "GCATCGATGAAGAACGCAGC",   # 5.8 S core
    "TCCTCCGCTTATTGATATGC"    # ITS4    (28 S start)
]
DEFAULT_MAX_MM = 4   # per-motif mismatch ceiling (default)
TIE_EPS = 0.3       # score delta below which we call it a tie
WRAP = 80           # output FASTA line length

# IUPAC equivalencies for edlib
IUPAC_EQUIV = [("Y", "C"), ("Y", "T"), ("R", "A"), ("R", "G"),
               ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T"),
               ("W", "A"), ("W", "T"), ("M", "A"), ("M", "C"),
               ("S", "C"), ("S", "G"), ("K", "G"), ("K", "T"),
               ("B", "C"), ("B", "G"), ("B", "T"),
               ("D", "A"), ("D", "G"), ("D", "T"),
               ("H", "A"), ("H", "C"), ("H", "T"),
               ("V", "A"), ("V", "C"), ("V", "G")]

# IUPAC complement & bitmasks
IUPAC_COMP = str.maketrans(
    "ACGTRYMKWSVHDBNacgtrymkwsvhdbn-",
    "TGCAYRKMWSBDHVNtgcayrkmwsbdhvn-"
)

# Precomputed bitmask lookup for IUPAC codes
I2M: Dict[str, int] = {
    **{b: 1 << i for i, b in enumerate("ACGT")},
    "R": 5, "Y": 10, "S": 6, "W": 9,
    "K": 12, "M": 3, "B": 14, "D": 13,
    "H": 11, "V": 7, "N": 15, "-": 0
}
I2M.update({k.lower(): v for k, v in I2M.items()})


def revcomp(seq: str) -> str:
    """Return reverse complement of a sequence."""
    return seq.translate(IUPAC_COMP)[::-1]


# Precompute reverse complements of motifs
REV_MOTIFS = [revcomp(m) for m in FWD_MOTIFS]


def compatible(a: str, b: str) -> bool:
    """Check if IUPAC symbols a and b share ≥ 1 concrete base."""
    return bool(I2M.get(a, 0) & I2M.get(b, 0))


@dataclass(frozen=True)
class Hit:
    """Represents a motif hit with mismatches and position."""
    mism: int
    pos: int


def best_hit(seq: str, motif: str, max_mm: int = DEFAULT_MAX_MM) -> Optional[Hit]:
    """
    Find best (fewest mismatches, earliest) fuzzy occurrence of motif in seq.

    Uses gapped alignment to properly handle insertions and deletions.
    """
    if len(seq) < len(motif) - max_mm:  # Sequence too short even with deletions
        return None

    seq_upper = seq.upper()
    motif_upper = motif.upper()

    # Use edlib for gapped alignment with IUPAC support
    result = edlib.align(
        motif_upper,
        seq_upper,
        mode='HW',  # Infix mode (like sliding window but with gaps)
        task='locations',
        k=max_mm,  # Maximum edit distance
        additionalEqualities=IUPAC_EQUIV
    )

    # Check if we found a match within the distance threshold
    if result['editDistance'] == -1 or result['editDistance'] > max_mm:
        return None

    # Find the best location (earliest start position for ties in edit distance)
    best_location = min(result['locations'], key=lambda loc: loc[0])

    return Hit(result['editDistance'], best_location[0])


@dataclass
class OrientationStats:
    """Statistics for orientation decision."""
    hits: List[Hit]
    total_mm: int
    best: Hit
    orientation: Optional[str] = None


def collect_stats(seq: str, motifs: List[str], max_mm: int = DEFAULT_MAX_MM, label: str = "") -> OrientationStats:
    """Collect statistics for a set of motifs."""
    hits = []
    for motif in motifs:
        hit = best_hit(seq, motif, max_mm)
        if hit:
            hits.append(hit)
    
    if not hits:
        return OrientationStats([], float('inf'), Hit(float('inf'), float('inf')))
    
    total_mm = sum(h.mism for h in hits)
    best = min(hits, key=lambda h: (h.mism, h.pos))
    
    return OrientationStats(hits, total_mm, best)


def decide_orientation(seq: str, max_mm: int = DEFAULT_MAX_MM, verbose: bool = False) -> Tuple[str, OrientationStats, OrientationStats]:
    """
    Determine sequence orientation based on motif matches.
    
    Returns: (orientation, forward_stats, reverse_stats)
    """
    fwd = collect_stats(seq, FWD_MOTIFS, max_mm, "forward")
    rev = collect_stats(seq, REV_MOTIFS, max_mm, "reverse")
    
    # No hits at all
    if not fwd.hits and not rev.hits:
        return "uncertain", fwd, rev
    
    # Different number of hits
    if len(fwd.hits) != len(rev.hits):
        ori = "forward" if len(fwd.hits) > len(rev.hits) else "reverse"
        return ori, fwd, rev
    
    # Different total mismatches
    if fwd.total_mm != rev.total_mm:
        ori = "forward" if fwd.total_mm < rev.total_mm else "reverse"
        return ori, fwd, rev
    
    # Final tie-breaker: earliest best hit (position scaled)
    score_f = fwd.best.mism * 10 + fwd.best.pos / 1000
    score_r = rev.best.mism * 10 + rev.best.pos / 1000
    
    if abs(score_f - score_r) < TIE_EPS:
        return "uncertain", fwd, rev
    
    ori = "forward" if score_f < score_r else "reverse"
    return ori, fwd, rev


def validate_sequence(seq: str, header: str) -> List[str]:
    """Validate sequence and return list of warnings."""
    warnings = []
    
    if not seq:
        warnings.append(f"Empty sequence for '{header}'")
        return warnings
    
    # Check for non-IUPAC characters
    invalid_chars = re.findall(r"[^ACGTRYMKWSVHDBNacgtrymkwsvhdbn\-]", seq)
    if invalid_chars:
        unique_chars = sorted(set(invalid_chars))
        warnings.append(f"Non-IUPAC symbols in '{header}': {', '.join(unique_chars)}")
    
    return warnings


def fasta_reader(stream, quiet: bool = False) -> Iterator[Tuple[str, str]]:
    """
    Read FASTA sequences from stream, handling stray '>' symbols.
    
    Yields: (header, sequence) tuples
    """
    header: Optional[str] = None
    sequence_parts: List[str] = []
    line_num = 0
    
    for raw_line in stream:
        line_num += 1
        line = raw_line.rstrip("\n\r")
        
        if line.startswith(">"):
            # Yield previous sequence if exists
            if header is not None:
                yield header, "".join(sequence_parts)
            
            # Start new sequence
            header = line[1:].strip() or f"<empty_header_line_{line_num}>"
            sequence_parts = []
        else:
            # Handle stray '>' symbols in sequence lines
            if ">" in line:
                parts = line.split(">")
                sequence_parts.append(parts[0].strip())
                
                # Each subsequent part becomes a new sequence
                for i, part in enumerate(parts[1:], 1):
                    logging.warning(
                        f"Line {line_num}: Stray '>' found, creating new sequence "
                        f"'{part.strip() or f'<empty_from_stray_{line_num}_{i}'}'"
                    )
                    if header is not None:
                        yield header, "".join(sequence_parts)
                    header = part.strip() or f"<empty_from_stray_{line_num}_{i}>"
                    sequence_parts = []
            else:
                # Normal sequence line
                sequence_parts.append(line.strip())
    
    # Don't forget the last sequence
    if header is not None:
        yield header, "".join(sequence_parts)


def write_fasta(header: str, seq: str, out=sys.stdout) -> None:
    """Write a FASTA sequence with proper line wrapping."""
    out.write(f">{header}\n")
    for chunk in wrap(seq, WRAP):
        out.write(chunk + "\n")


def process_file(
    input_stream,
    output_stream,
    max_mm: int = DEFAULT_MAX_MM,
    dry_run: bool = False,
    verbose: bool = False,
    stats_only: bool = False,
    report_reversed: bool = True,
    quiet: bool = False
) -> Dict[str, int]:
    """
    Process a FASTA file and return statistics.
    
    Returns: Dictionary with counts of forward, reverse, uncertain sequences
    """
    counts = defaultdict(int)
    reversed_sequences = []
    uncertain_sequences = []
    
    for header, seq in fasta_reader(input_stream, quiet=quiet):
        counts['total'] += 1
        
        # Validate sequence
        warnings = validate_sequence(seq, header)
        if not quiet:
            for warning in warnings:
                logging.warning(warning)
        
        if not seq:
            counts['empty'] += 1
            continue
        
        # Determine orientation
        orientation, fwd_stats, rev_stats = decide_orientation(seq, max_mm, verbose)
        counts[orientation] += 1
        
        # Track reversed and uncertain sequences
        if orientation == "reverse":
            reversed_sequences.append(header)
        elif orientation == "uncertain":
            uncertain_sequences.append(header)
        
        # Verbose output
        if verbose:
            logging.info(
                f"{header}: {orientation} "
                f"(fwd: {len(fwd_stats.hits)} hits, {fwd_stats.total_mm} mm; "
                f"rev: {len(rev_stats.hits)} hits, {rev_stats.total_mm} mm)"
            )
        
        # Write output (unless dry run or stats only)
        if not dry_run and not stats_only:
            if orientation == "reverse":
                write_fasta(header, revcomp(seq).upper(), output_stream)
            else:
                # Forward or uncertain - keep as is
                write_fasta(header, seq.upper(), output_stream)
        
        # Log uncertain orientations
        if orientation == "uncertain" and not quiet:
            logging.warning(f"Uncertain orientation for '{header}'")
    
    # Report reversed sequences if requested
    if report_reversed and not quiet:
        if reversed_sequences:
            print(f"\n=== Reversed sequences ({len(reversed_sequences)}) ===", file=sys.stderr)
            for header in reversed_sequences:
                print(f"  {header}", file=sys.stderr)
        
        if uncertain_sequences:
            print(f"\n=== Uncertain sequences ({len(uncertain_sequences)}) ===", file=sys.stderr)
            for header in uncertain_sequences:
                print(f"  {header}", file=sys.stderr)
    
    return dict(counts)


def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description="Auto-orient fungal ITS sequences using conserved motifs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        "files",
        nargs="*",
        help="Input FASTA files (default: stdin)"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output file (default: stdout)"
    )
    parser.add_argument(
        "-n", "--dry-run",
        action="store_true",
        help="Don't write output, just analyze"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output to stderr"
    )
    parser.add_argument(
        "-s", "--stats",
        action="store_true",
        help="Print orientation statistics"
    )
    parser.add_argument(
        "--stats-only",
        action="store_true",
        help="Only print statistics, no sequence output"
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=DEFAULT_MAX_MM,
        help=f"Maximum mismatches per motif (default: {DEFAULT_MAX_MM})"
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress all warnings and reports (quiet mode)"
    )
    
    args = parser.parse_args()
    
    # Configure logging
    if args.quiet:
        log_level = logging.ERROR  # Only show errors in quiet mode
    elif args.verbose:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING
    
    logging.basicConfig(
        level=log_level,
        format="[%(levelname)s] %(message)s",
        stream=sys.stderr
    )
    
    try:
        # Setup input/output
        input_stream = fileinput.input(args.files) if args.files else sys.stdin
        
        from contextlib import nullcontext

        output_cm = (
            open(args.output, "w")          # real file
            if args.output and not (args.dry_run or args.stats_only)
            else nullcontext(sys.stdout)    # acts as a no-op context manager
        )

        with output_cm as output_stream:
            # Process sequences
            stats = process_file(
                input_stream,
                output_stream,
                max_mm=args.max_mismatches,
                dry_run=args.dry_run,
                verbose=args.verbose,
                stats_only=args.stats_only,
                report_reversed=not args.quiet,
                quiet=args.quiet,
            )
        
        # Print statistics if requested
        if (args.stats or args.stats_only) and not args.quiet:
            print("\n=== Orientation Statistics ===", file=sys.stderr)
            print(f"Total sequences: {stats.get('total', 0)}", file=sys.stderr)
            print(f"Forward:         {stats.get('forward', 0)}", file=sys.stderr)
            print(f"Reverse:         {stats.get('reverse', 0)}", file=sys.stderr)
            print(f"Uncertain:       {stats.get('uncertain', 0)}", file=sys.stderr)
            print(f"Empty:           {stats.get('empty', 0)}", file=sys.stderr)
        
    except KeyboardInterrupt:
        logging.error("Aborted by user (KeyboardInterrupt)")
        sys.exit(130)
    except Exception as exc:
        logging.error(f"{exc.__class__.__name__}: {exc}")
        sys.exit(1)


if __name__ == "__main__":
    main()
