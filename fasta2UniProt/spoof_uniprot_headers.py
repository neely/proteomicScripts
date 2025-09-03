# -----------------------------------------------------------------------------
# INSTRUCTIONS:
# This script processes all FASTA files in the same directory where it is run.
#
# HOW TO RUN:
# 1. Place this script in a folder containing the FASTA files you want to process.
# 2. Open your terminal or command prompt.
# 3. Navigate to that folder.
# 4. Run the script with the following command:
#    python spoof_uniprot_headers.py
#
# The script will automatically find all files with a '.fasta' extension,
# process each one, and save the output as a new file with the suffix
# '_uniprot_spoof.fasta'.
#
# -----------------------------------------------------------------------------
#
# This script takes a FASTA file with headers in the format:
# >XP_024421550.2 lactotransferrin [Desmodus rotundus] GN=LTF
# and transforms them into a UniProt-like format:
# >rs|XP_024421550.2|LTF_DESRO lactotransferrin OS=Desmodus rotundus GN=LTF
#
# The script extracts the necessary information (accession, gene symbol,
# and species name) to construct the new formatted header. It processes
# each input FASTA file and creates a new corresponding output file.

import sys
import re
import os
import glob

def parse_and_format_header(header_line):
    """
    Parses an original RefSeq-formatted FASTA header and converts it to
    the new UniProt-like format.
    """
    # Remove the leading '>'
    header = header_line.lstrip('>')
    
    # Use regular expressions to extract key components
    # Pattern: accession, description, species name in brackets, and GN=symbol
    match = re.search(r'^(\S+)\s+(.+?)\s+\[(.+?)\]\s+GN=(.+)$', header)
    
    if not match:
        # If the expected pattern is not found, return the original header with a warning
        print(f"Warning: Could not parse header, skipping formatting for: {header_line.strip()}", file=sys.stderr)
        return header_line.strip()
        
    accession = match.group(1)
    description = match.group(2)
    species_full = match.group(3)
    gene_symbol = match.group(4)
    
    # Extract species code from full species name (e.g., "Desmodus rotundus" -> "DESRO")
    species_parts = species_full.split()
    if len(species_parts) >= 2:
        species_code = f"{species_parts[0][:3].upper()}{species_parts[1][:2].upper()}"
    else:
        # Fallback if species name is not in the expected format
        species_code = "UNKN"

    # Construct the new header in the desired format
    # Example: >rs|XP_024421550.2|LTF_DESRO lactotransferrin OS=Desmodus rotundus GN=LTF
    new_header = (
        f"rs|{accession}|{gene_symbol}_{species_code} "
        f"{description} OS={species_full} GN={gene_symbol}"
    )
    
    return f">{new_header}"

def process_fasta(input_fasta, output_fasta):
    """
    Reads the input FASTA file, transforms headers, and writes to the output file.
    """
    try:
        with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    # Process header lines
                    new_line = parse_and_format_header(line.strip())
                    outfile.write(new_line + '\n')
                else:
                    # Write sequence lines as is
                    outfile.write(line)
        print(f"Successfully formatted headers from '{input_fasta}' and saved to '{output_fasta}'")
        return True
    except FileNotFoundError:
        print(f"Error: The file '{input_fasta}' was not found.")
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
        
        process_fasta(input_file, output_file)
        
    print("\nBatch processing complete.")

if __name__ == "__main__":
    main()
