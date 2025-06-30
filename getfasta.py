#!/usr/bin/env python3

"""
Mycomap FASTA downloader

By Alan Rockefeller - June 30, 2025.    Version 1.1

"""

import urllib.request
import urllib.parse
import sys
import re
import time

# Set a timeout for network requests in seconds
REQUEST_TIMEOUT = 10

if len(sys.argv) < 2:
    print("Usage: python getfasta.py <MycoBLAST-result-URL>")
    sys.exit(1)

# Helper
def count_fasta_seqs(fasta_bytes: bytes) -> int:
    """Return the number of sequences in a FASTA file (bytes)."""
    return sum(1 for line in fasta_bytes.splitlines() if line.startswith(b'>'))

# Parse URL
blast_id_match = re.search(r'r(\d+)', sys.argv[1])
if not blast_id_match:
    print("Could not find an r<digits> pattern in the URL you provided.")
    sys.exit(1)

blast_id = blast_id_match.group(1)
print(f"Downloading FASTA files for MycoBLAST ID: {blast_id}")

# HTTP opener
base_url = "https://mycomap.com/index.php"
opener = urllib.request.build_opener()
opener.addheaders = [
    ('User-Agent', 'curl/7.68.0'),
    ('Accept', '*/*')
]

total_start = time.time()

# NCBI download
params_ncbi = urllib.parse.urlencode({
    'app': 'genbank',
    'module': 'genbank',
    'controller': 'blast',
    'do': 'fasta',
    'id': blast_id
})
url_ncbi = f"{base_url}?{params_ncbi}"
t0 = time.time()

try:
    with opener.open(url_ncbi, timeout=REQUEST_TIMEOUT) as resp:
        ncbi_bytes = resp.read()
    ncbi_time = time.time() - t0

    try:
        with open(f'ncbi_{blast_id}.fasta', 'wb') as fh:
            fh.write(ncbi_bytes)

        ncbi_seqs = count_fasta_seqs(ncbi_bytes)
        print(f"NCBI downloaded in {ncbi_time:.2f}s "
              f"({len(ncbi_bytes)} bytes, {ncbi_seqs} sequences)")
    except IOError as e:
        print(f"Error writing to file ncbi_{blast_id}.fasta: {e}")

except urllib.error.URLError as e:
    print(f"Error downloading NCBI data: {e}")
except TimeoutError:
    print(f"Error: The request for NCBI data timed out after {REQUEST_TIMEOUT} seconds.")


# MycoBLAST download
params_myco = urllib.parse.urlencode({
    'app': 'genbank',
    'module': 'genbank',
    'controller': 'blast',
    'do': 'localFasta',
    'id': blast_id
})
url_myco = f"{base_url}?{params_myco}"
t0 = time.time()

try:
    with opener.open(url_myco, timeout=REQUEST_TIMEOUT) as resp:
        myco_bytes = resp.read()
    myco_time = time.time() - t0

    try:
        with open(f'myco_{blast_id}.fasta', 'wb') as fh:
            fh.write(myco_bytes)

        myco_seqs = count_fasta_seqs(myco_bytes)
        print(f"MycoBLAST downloaded in {myco_time:.2f}s "
              f"({len(myco_bytes)} bytes, {myco_seqs} sequences)")
    except IOError as e:
        print(f"Error writing to file myco_{blast_id}.fasta: {e}")

except urllib.error.URLError as e:
    print(f"Error downloading MycoBLAST data: {e}")
except TimeoutError:
    print(f"Error: The request for MycoBLAST data timed out after {REQUEST_TIMEOUT} seconds.")
