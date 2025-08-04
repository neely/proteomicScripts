"""
convert_fasta_headers.py

Converts FASTA headers into a UniProt-like format.

## INPUT FORMAT:
    >ENST00000343518.POTEH.43096 JAUMIU010000140:8603-8709
    >ENSMUST00000343518.POTEH.43096 JAUMIU010000140:8703-8809
    >TMEM129-like.CM061257:895302-901188

## OUTPUT FORMAT:
    >sp|ENST00000343518|POTEH_TABRA ENST00000343518 GN=POTEH
    >sp|ENST00000343518_2|TMEM129-like_TABRA ENST00000343518_2 GN=TMEM129-like

## RULES:
1. If header starts with "ENST" or "ENSMUST":
   - Extract **ENST/ENSMUST ID** (before first `.`).
   - Extract **GN** (between first and second `.`).
   - If duplicated, **add a counter to ENST/ENSMUST**.

2. If the header **does NOT** start with "ENST"/"ENSMUST":
   - Assign it the **last seen ENST/ENSMUST** + `-X` counter.
   - Extract **GN from the first word before `.`** (preserving `-like` if present).

## USAGE:
    python convert_fasta_headers.py input.fasta output.fasta
"""

import sys
from collections import defaultdict

def parse_header(header, enst_counters, last_enst):
    """Parse and transform a FASTA header into the desired format."""
    parts = header.lstrip('>').split()
    main_id = parts[0]  # Extract first part before spaces

    if main_id.startswith(("ENST", "ENSMUST")):
        # Extract ENST/ENSMUST ID and GN
        segments = main_id.split('.')
        if len(segments) < 3:
            raise ValueError(f"Unexpected header format: {header}")

        unique_id = segments[0]  # ENST00000343518
        gene_name = segments[1]   # POTEH

        # Track unique instances of each ENST/ENSMUST ID
        enst_counters[unique_id] += 1
        if enst_counters[unique_id] > 1:
            unique_id = f"{unique_id}_{enst_counters[unique_id]}"

        last_enst[0] = unique_id  # Store last seen ENST/ENSMUST

    else:
        # Handle headers that do NOT start with ENST/ENSMUST
        if not last_enst[0]:
            raise ValueError(f"No ENST/ENSMUST ID found before: {header}")

        # Assign last seen ENST + counter
        enst_counters[last_enst[0]] += 1
        unique_id = f"{last_enst[0]}_{enst_counters[last_enst[0]]}"

        # Extract GN from first part of ID (preserving `-like` if present)
        gene_name = main_id.split('.')[0]  # Take everything before the first "."

    # Construct new header
    new_protein_id = f"{gene_name}_TABRA"
    new_header = f"sp|{unique_id}|{new_protein_id} {unique_id} GN={gene_name}"

    return new_header

def process_fasta(input_fasta, output_fasta):
    """Reads the input FASTA, transforms headers, and writes to output."""
    enst_counters = defaultdict(int)  # Tracks ENST/ENSMUST counters
    last_enst = [None]  # Stores last seen ENST/ENSMUST ID
    sequence_lines = []

    with open(input_fasta, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                new_header = parse_header(line, enst_counters, last_enst)
                sequence_lines.append(f">{new_header}")  # ✅ Ensure ">" is present
            else:
                sequence_lines.append(line)

    # Write to output FASTA
    with open(output_fasta, 'w') as outfile:
        for line in sequence_lines:
            outfile.write(line + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_fasta_headers.py <input.fasta> <output.fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    process_fasta(input_fasta, output_fasta)
    print(f"Converted FASTA written to {output_fasta}")
