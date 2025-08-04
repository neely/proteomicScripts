"""
This script updates the headers of a RefSeq FASTA file. It reads a tab-separated value (TSV) feature table file (available on NCBI where you got the fasta), mapping protein accession numbers to gene symbols. For each sequence in the input FASTA file, the script identifies the accession number in the header and searches for a matching entry in the TSV file. If a match is found, it appends the corresponding gene symbol to the end of the FASTA header in the format "GN=<symbol>", similar to how UniProt FASTA headers are formatted. The script then writes the sequences with their modified headers to a new output FASTA file.

Usage:
To use, run the following command from your terminal in the directory containing the files:
python modify_fasta_headers.py <input.fasta> <input.tsv> <output.fasta>

- <input.fasta>: The original FASTA file.
- <input.tsv>: The feature table TSV file (will contain 'product_accession' and 'symbol' columns).
- <output.fasta>: The name for the new FASTA file with updated headers.
"""

import sys

# Read in the input files
fasta_file = sys.argv[1]
tsv_file = sys.argv[2]

# Create a dictionary to map accession numbers to symbols
accession_to_symbol = {}
with open(tsv_file) as tsv:
    # Read the header row and get the index of the product_accession and symbol columns
    header = tsv.readline().rstrip().split('\t')
    product_accession_index = header.index('product_accession')
    symbol_index = header.index('symbol')
    # Iterate through the remaining rows and populate the dictionary
    for line in tsv:
        row = line.rstrip().split('\t')
        accession = row[product_accession_index]
        symbol = row[symbol_index]
        if symbol:
            accession_to_symbol[accession] = symbol

# Process the FASTA file and write the modified headers to a new file
with open(fasta_file) as fasta, open(sys.argv[3], 'w') as output:
    accession = None
    symbol = None
    for line in fasta:
        if line.startswith('>'):  # header line
            header = line.rstrip()
            accession = header.split()[0][1:]
            if accession in accession_to_symbol:
                symbol = accession_to_symbol[accession]
                header += f' GN={symbol}'
            output.write(header + '\n')
        else:  # sequence line
            output.write(line)
