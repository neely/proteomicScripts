# -----------------------------------------------------------------------------
# INSTRUCTIONS:
# This script processes all FASTA files in the same directory where it is run.
# It is designed to convert Ensembl-style headers to a UniProt-like format.
#
# HOW TO RUN:
# 1. Place this script in a folder containing the FASTA files you want to process.
# 2. Open your terminal or command prompt.
# 3. Navigate to that folder.
# 4. Run the script with the following command:
#    python spoof_uniprot_ensembl_headers.py
#
# The script will automatically find all files with a '.fasta' extension,
# process each one, and save the output as a new file with the suffix
# '_uniprot_spoof.fasta'.
#
# -----------------------------------------------------------------------------
#
# This script takes an Ensembl-formatted FASTA file with headers like this:
# >ENSPCAT00000000047.1 ... gene_symbol:PPP1R3B description:protein phosphatase...
#
# And transforms them into a UniProt-like format:
# >sp|ENSPCAT00000000047.1|PPP1R3B_PROCA protein phosphatase 1 regulatory subunit... OS=Procavia capensis GN=PPP1R3B
#
# The script performs a two-pass process on each file:
# Pass 1: Identifies unique protein accessions and handles missing data.
#         - If a 'description' is missing, it is set to "unreported".
#         - If a 'gene_symbol' is missing, a unique identifier like 'LOC1' is generated.
#         - If a protein accession is duplicated, a '-counter' is added to make it unique.
# Pass 2: Reconstructs the new headers and writes the final FASTA file.
# The species name is hardcoded as 'Procavia capensis'.

import sys
import re
import os
import glob
from collections import Counter

# Manually defined species for the output headers
SPECIES_FULL = "Procavia capensis"
SPECIES_CODE = "PROCA"

def process_ensembl_fasta(input_fasta, output_fasta):
    """
    Reads the input Ensembl FASTA file, processes headers in a two-pass
    approach, and writes a new FASTA file with UniProt-like headers.
    """
    print(f"Processing '{input_fasta}'...")
    headers_and_sequences = []
    accession_list = []
    
    # Pass 1: Read and analyze all headers and store sequences
    try:
        with open(input_fasta, 'r') as infile:
            current_header = None
            current_sequence_lines = []
            
            for line in infile:
                if line.startswith('>'):
                    if current_header:
                        headers_and_sequences.append({
                            'header': current_header,
                            'sequence_lines': current_sequence_lines
                        })
                    
                    header_line = line.lstrip('>').strip()
                    parts = header_line.split()
                    
                    accession = parts[0]
                    accession_list.append(accession)
                    
                    current_header = {
                        'original_line': line,
                        'accession': accession,
                        'description': "unreported",
                        'gene_symbol': None
                    }

                    # Parse description and gene symbol
                    gene_symbol_match = re.search(r'gene_symbol:(\S+)', header_line)
                    if gene_symbol_match:
                        current_header['gene_symbol'] = gene_symbol_match.group(1)
                    
                    description_match = re.search(r'description:(.*?)(?: \[|$)', header_line)
                    if description_match:
                        current_header['description'] = description_match.group(1).strip()

                    current_sequence_lines = []
                else:
                    current_sequence_lines.append(line)
            
            # Append the last entry after the loop finishes
            if current_header:
                headers_and_sequences.append({
                    'header': current_header,
                    'sequence_lines': current_sequence_lines
                })

    except FileNotFoundError:
        print(f"Error: The file '{input_fasta}' was not found.", file=sys.stderr)
        return False
    
    # Handle duplicates and missing gene symbols
    accession_counts = Counter(accession_list)
    dup_counters = {}
    loc_counter = 1
    
    for entry in headers_and_sequences:
        header = entry['header']
        
        # Handle duplicate accessions
        if accession_counts[header['accession']] > 1:
            if header['accession'] not in dup_counters:
                dup_counters[header['accession']] = 0
            
            dup_counters[header['accession']] += 1
            header['accession'] = f"{header['accession']}-{dup_counters[header['accession']]}"
        
        # Handle missing gene symbol
        if not header['gene_symbol']:
            header['gene_symbol'] = f"LOC{loc_counter}"
            loc_counter += 1
    
    # Pass 2: Write the new FASTA file with formatted headers
    try:
        with open(output_fasta, 'w') as outfile:
            for entry in headers_and_sequences:
                header = entry['header']
                
                new_header = (
                    f"sp|{header['accession']}|{header['gene_symbol']}_{SPECIES_CODE} "
                    f"{header['description']} OS={SPECIES_FULL} GN={header['gene_symbol']}"
                )
                
                outfile.write(f">{new_header}\n")
                # Write each sequence line as-is to preserve original line breaks
                for seq_line in entry['sequence_lines']:
                    outfile.write(seq_line)
                
        print(f"Successfully formatted headers from '{input_fasta}' and saved to '{output_fasta}'")
        return True
    
    except Exception as e:
        print(f"Error writing to file '{output_fasta}': {e}", file=sys.stderr)
        return False

def main():
    """
    Main function to find and process all FASTA files in the current directory.
    """
    # Get all files with a .fasta extension in the current directory
    fasta_files = glob.glob('*.fasta')
    
    if not fasta_files:
        print("No .fasta files found in the current directory.")
        sys.exit(0)
    
    print(f"Found {len(fasta_files)} FASTA files to process.")
    
    for input_file in fasta_files:
        # Create a new output filename with a _uniprot_spoof.fasta suffix
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_uniprot_spoof.fasta"
        
        process_ensembl_fasta(input_file, output_file)
        
    print("\nBatch processing complete.")

if __name__ == "__main__":
    main()
