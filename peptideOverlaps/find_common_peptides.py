# This script analyzes text files containg peptide sequences in a specified directory (which is
# assumed to be the directory where the script is located). It is designed to process these files,
# treating each file as a collection of peptide sequences.
# 
# The primary functions of the script are:
# 
# 1. Peptide Extraction: It reads each .txt file in the directory, skipping the first line 
# (assumed to be a header), and extracts peptide sequences from the first column of each 
# subsequent line.
# 2. Counting and Reporting: It calculates and prints the total number of unique peptides 
# found in each file.
# 3. Overlap Analysis: It determines the overlap (common peptides) between all possible 
# combinations of the input files.
# 4. Top Overlap Display: It identifies and prints the top 10 combinations of files with the 
# highest number of shared peptides.
# 5. Data Visualization: It generates a bar chart visualizing the number of peptides per file 
# and saves this plot as a PNG file named peptides_per_file.png in the script's directory.
# 
# How to Run the Script from the Command Line:
#    navigate to the directory containing txt files
#    find_common_peptides.py
# 

import os
import sys
from collections import defaultdict
from itertools import combinations

# Function to extract peptide sequences from a file
def extract_sequences(file_path):
    sequences = set()
    with open(file_path, 'r') as file:
        next(file)  # Skip header assuming first line is header
        for line in file:
            line = line.strip()
            if line:
                sequence = line.split()[0]  # Assuming the sequence is in the first column
                sequences.add(sequence)
    return sequences

if __name__ == "__main__":
    # Determine the directory containing this script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Directory containing your text files (same as script's directory)
    directory = script_dir

    # Dictionary to store sequences per file
    sequences_per_file = defaultdict(set)

    # Loop through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            file_path = os.path.join(directory, filename)
            sequences = extract_sequences(file_path)
            sequences_per_file[filename] = sequences

    # Calculate number of peptides per file
    num_peptides = {file: len(sequences) for file, sequences in sequences_per_file.items()}

    # Print results
    print("Number of peptides per file:")
    for file, num in num_peptides.items():
        print(f"{file}: {num}")
    
    # Determine overlap for all combinations of files
    overlap_combinations = defaultdict(int)
    file_list = list(sequences_per_file.keys())

    for i in range(2, len(file_list) + 1):  # Iterate over all combination sizes from 2 to len(file_list)
        for comb in combinations(file_list, i):
            overlap_size = len(set.intersection(*(sequences_per_file[file] for file in comb)))
            overlap_combinations[comb] = overlap_size

    # Sort combinations by overlap size in descending order
    sorted_combinations = sorted(overlap_combinations.items(), key=lambda x: x[1], reverse=True)

    # Print top 10 overlap combinations
    print("\nTop 10 Highest Overlap Combinations:")
    for i, (files, overlap_size) in enumerate(sorted_combinations[:10], 1):
        print(f"{i}. {files}: {overlap_size} common peptides")

    # Plotting
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.bar(range(len(num_peptides)), list(num_peptides.values()), label='Peptides per File')
    ax.set_xticks(range(len(num_peptides)))
    ax.set_xticklabels(num_peptides.keys(), rotation=45, ha='right')
    ax.set_xlabel('Files')
    ax.set_ylabel('Number of Peptides')
    ax.set_title('Peptides per File')
    ax.legend()

    # Save plot as PNG in the working directory
    plt.tight_layout()
    plt.savefig(os.path.join(script_dir, 'peptides_per_file.png'))
    plt.close()
