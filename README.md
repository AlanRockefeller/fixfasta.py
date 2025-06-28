# fixfasta.py

Version 1.0 - June 28, 2025

A tool for automatically fixing the orientation of fungal ITS sequences in FASTA files. Perfect for cleaning up sequences downloaded from databases like MycoMap or GenBank where some sequences might be reverse-complemented.

## What it does

Ever downloaded a bunch of ITS sequences only to find that some are in the wrong orientation, making your phylogenetic tree completely useless until you manually reverse-complement them? This tool automatically detects and fixes that problem by looking for conserved motifs in the ITS region. It flips sequences that are backwards, giving you a clean FASTA file ready for phylogenetic analysis.

## Features

- **Smart orientation detection** using three conserved ITS motifs (ITS1-F, 5.8S core, ITS4)
- **IUPAC-aware fuzzy matching** - handles ambiguous nucleotides and allows up to 4 mismatches
- **Automatic repair** of malformed FASTA files with stray '>' symbols
- **Detailed reporting** of which sequences were reversed (or silent operation with `-q`)
- **Fast and efficient** - processes large files quickly
- Uses orientation logic that breaks ties intelligently

## Installation

Just download the script and make it executable:

```bash
wget https://raw.githubusercontent.com/yourusername/fixfasta/main/fixfasta.py
chmod +x fixfasta.py
```

Requirements:
- Python 3.6+
- No external dependencies - uses only Python standard library

## Quick Start

Basic usage - fix orientation and see what was changed:
```bash
./fixfasta.py sequences.fasta > fixed_sequences.fasta
```

The script will report which sequences it reversed:
```
=== Reversed sequences (8) ===
  AF430271_Russula_fragilis
  KC885965_Amanita_muscaria
  MH856042_Cortinarius_sp
  ...
```

## Usage Examples

**Process multiple files:**
```bash
cat *.fasta | ./fixfasta.py > all_fixed.fasta
```

**Silent mode (no reports):**
```bash
./fixfasta.py input.fasta -q > output.fasta
```

**See detailed statistics:**
```bash
./fixfasta.py input.fasta --stats > output.fasta
```

**Verbose mode to understand decisions:**
```bash
./fixfasta.py input.fasta -v > output.fasta
```

**Dry run (analyze without modifying):**
```bash
./fixfasta.py input.fasta --dry-run
```

**Save to a specific output file:**
```bash
./fixfasta.py input.fasta -o output.fasta
```

## Command Line Options

```
-h, --help            Show help message
-o, --output          Output file (default: stdout)
-n, --dry-run        Don't write output, just analyze
-v, --verbose        Verbose output showing decision process
-s, --stats          Print orientation statistics
--stats-only         Only print statistics, no sequence output
-q, --quiet          Suppress all warnings and reports
--max-mismatches N   Maximum mismatches per motif (default: 4)
```

## How It Works

The tool uses three conserved motifs commonly found in fungal ITS sequences:

1. **ITS1-F** (TCCGTAGGTGAACCTGCGG) - found at the 18S end
2. **5.8S core** (GCATCGATGAAGAACGCAGC) - middle region
3. **ITS4** (TCCTCCGCTTATTGATATGC) - found at the 28S start

For each sequence, it:
1. Searches for these motifs in both orientations (forward and reverse-complement)
2. Counts how many motifs are found and how many mismatches each has
3. Makes a decision based on:
   - Which orientation has more motif hits
   - If tied, which has fewer total mismatches
   - If still tied, which has the earliest best hit
4. Reverse-complements sequences that are backwards

The fuzzy matching allows up to 4 mismatches per 20bp motif and understands IUPAC ambiguity codes (R, Y, S, W, K, M, etc.).

## Real-World Use Case

This tool was originally created to process ITS sequences downloaded from MycoMap for phylogenetic tree construction. It's particularly useful when combining sequences from multiple sources where orientation consistency isn't guaranteed.

Example workflow:
```bash
# Download sequences from MycoMap
# ... download process ...

# Fix orientations
./fixfasta.py mycomap.fa ncbi.fa > sequences_oriented.fas

# Now ready for MAFFT alignment, phylogeny.fr, RaxML, etc
mafft sequences_oriented.fasta > aligned.fasta
```

## Tips

- The default behavior shows you what was changed - use `-q` for silent operation in pipelines
- Use `--dry-run` first on new datasets to see what would be changed
- The tool preserves sequence names exactly 
- Handles messy FASTA files gracefully (like those with stray '>' symbols)
- All diagnostic output goes to stderr, so stdout piping remains clean

## License

This project is licensed under the MIT License:

```
MIT License

Copyright (c) 2025 Alan Rockefeller

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Contributing

Found a bug? Have a suggestion? Feel free to open an issue or submit a pull request!

https://github.com/AlanRockefeller/fixfasta.py

## Acknowledgments

Thanks to the mycological community for providing the data that made this tool necessary, and to everyone who's contributed sequences to public databases.

