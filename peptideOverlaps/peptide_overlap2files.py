# This script is designed to compare two files containing peptide sequences and calculate 
# the number of overlapping (common) sequences between them. The script find_common_peptides.py
# will list all combinations, but this will do just 2 files. I needed it for a quick check.
# 
# Each input file is expected to contain a list of unique peptide sequences, with one sequence 
# per line. The script reads the sequences from both files, treats them as sets, and computes 
# the size of their intersection to determine the overlap. The output is in the console.
# 
# How to Run the Script from the Command Line:
# run the following from the directory containing txt files
#    python peptide_overlap.py <file1.txt> <file2.txt>
# 
#    Replace <file1.txt> and <file2.txt> with the actual names of your input files.

import sys
# python peptide_overlap.py file1.txt file2.txt


def read_peptides(file_path):
    """
    Reads peptide sequences from a file and returns a set of unique sequences.
    Each line in the file is assumed to contain a single peptide sequence.
    """
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def main():
    """
    Main function to handle command-line arguments, read peptide sequences,
    compute the overlap, and print the result.
    """
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python peptide_overlap.py <file1.txt> <file2.txt>")
        print("This script compares two files containing peptide sequences and reports the number of overlapping sequences.")
        sys.exit(1)

    # Get file names from command-line arguments
    file1, file2 = sys.argv[1], sys.argv[2]

    # Read peptide sequences from both files
    peptides1 = read_peptides(file1)
    peptides2 = read_peptides(file2)

    # Compute the intersection (overlap) of the two sets
    overlap = peptides1.intersection(peptides2)

    # Print the result
    print(f"Number of overlapping peptides: {len(overlap)}")

if __name__ == "__main__":
    main()
