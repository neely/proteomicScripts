# This script uses text files in a directory that are each a list of peptides identified in a search. 
# It treats each list as a set and calculates the average sequence similarity between the peptides 
# in each pair of peptide lists by calculating sequence similarity using Levenshtein distance 
# for peptide sequences. The script was modified to use parallel processing since the number of 
# comparisons is quite high. Also, many attempts were made to include some sort of progress percent 
# (which is based on determining how many calculations will be performed by how many have been made), 
# but this isn't updated frequently as this seemed to cause too much performance hit (at scale) 
# than was benefiting. The distance matrix is also printed since other forms of clustering 
# (not simply heirarchical and ward method) may be preferred. Note this can get really big really 
# fast, and I am not sure it is useful.
#
# pip install numpy scipy matplotlib python-Levenshtein
# run the following from the directory containing txt files
# python peptide_clustering.py peptide_clustering_dendrogram.png
#
#
import os
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from Levenshtein import distance as levenshtein_distance
import argparse
from multiprocessing import Pool, cpu_count

def read_peptide_list(file_path):
    with open(file_path, 'r') as file:
        peptides = file.read().strip().split('\n')[1:]  # Ignore header
    return peptides

def compute_pairwise_distance(params):
    list1, list2 = params
    total_distance = 0
    comparisons = 0
    for pep1 in list1:
        for pep2 in list2:
            total_distance += levenshtein_distance(pep1, pep2)
            comparisons += 1
    return total_distance / comparisons if comparisons != 0 else 1

def main(output_file):
    # Step 1: Read peptide lists from files in the current directory
    peptide_lists = []
    file_names = [f for f in os.listdir('.') if f.endswith('.txt')]
    for file_name in file_names:
        peptides = read_peptide_list(file_name)
        peptide_lists.append(peptides)
    
    # Step 2: Compute the average Levenshtein distance matrix in parallel
    num_samples = len(peptide_lists)
    distance_matrix = np.zeros((num_samples, num_samples))
    total_comparisons = sum(len(list1) * len(list2) for i, list1 in enumerate(peptide_lists) for list2 in peptide_lists[i+1:])
    comparison_count = 0

    # Create a list of pairs for which we need to compute distances
    pairs = [(peptide_lists[i], peptide_lists[j]) for i in range(num_samples) for j in range(i + 1, num_samples)]
    total_pairs = len(pairs)

    # Use multiprocessing Pool to compute distances in parallel
    with Pool(cpu_count()) as pool:
        distances = []
        for i, result in enumerate(pool.imap(compute_pairwise_distance, pairs, chunksize=10)):
            distances.append(result)
            comparison_count += len(pairs[i][0]) * len(pairs[i][1])
            if i % 10 == 0 or i == total_pairs - 1:  # Reduce frequency of progress updates
                progress = (comparison_count / total_comparisons) * 100
                print(f"Progress: {progress:.2f}% ({comparison_count:,} of {total_comparisons:,})", end='\r')

    # Fill the distance matrix with the computed distances
    index = 0
    for i in range(num_samples):
        for j in range(i + 1, num_samples):
            distance_matrix[i, j] = distances[index]
            distance_matrix[j, i] = distances[index]
            index += 1

    print()  # Move to the next line after the progress bar

    # Print the distance matrix for verification
    print("Distance Matrix:")
    print(distance_matrix)

    # Print the total number of pairwise comparisons
    print("Total Number of Pairwise Comparisons:", total_comparisons)

    # Step 3: Perform hierarchical clustering
    condensed_distance_matrix = squareform(distance_matrix)
    Z = linkage(condensed_distance_matrix, method='ward')

    # Print the linkage matrix for verification
    print("Linkage Matrix:")
    print(Z)

    # Step 4: Plot the dendrogram
    plt.figure(figsize=(10, 7))
    dendrogram(Z, labels=[os.path.splitext(f)[0] for f in file_names])
    plt.title('Dendrogram of Peptide Samples with Sequence Similarity')
    plt.xlabel('Sample')
    plt.ylabel('Distance')
    plt.savefig(output_file)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster peptide lists and generate a dendrogram.")
    parser.add_argument('output_file', type=str, help="Output file for the dendrogram image.")
    args = parser.parse_args()
    
    main(args.output_file)
