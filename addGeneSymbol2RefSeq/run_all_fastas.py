import os
import glob
import gzip
import subprocess
import sys

def find_file_pairs(directory="."):
    """
    Finds pairs of gzipped protein and feature table files in a given directory.
    This will run the modify_fasta_headers.py to generate a header with GN= at the end.
    The pairs are identified by a common accession number in their filenames.
    These files are from the NCBI RefSeq FTP, since the file names contain info on refseq aasembly name.
    Assumes that when the protein.faa.gz file was downloaded, some identifier was added as a prefix.
    Ex. RefSeq_N_fur_seal_GCF*.protein.faa.gz
    This is used to name the resulting file, but the feature_table.txt.gz does not need to be renamed.
    This script is run by being placed in the directory of the pairs, along with modify_fasta_headers.py and executing:
    python run_all_fastas.py
    """
    file_pairs = {}
    
    # Get all gzipped files in the specified directory
    all_files = glob.glob(os.path.join(directory, "*.gz"))
    
    # Process each file to identify its type and accession number
    for file_path in all_files:
        filename = os.path.basename(file_path)
        
        # Check if it's a protein FASTA file
        if "_protein.faa.gz" in filename:
            accession_key = filename.split("_protein.faa.gz")[0].split("-", 2)[2]
            if accession_key not in file_pairs:
                file_pairs[accession_key] = {"protein": [], "feature": []}
            file_pairs[accession_key]["protein"].append(file_path)
            
        # Check if it's a feature table file
        elif "_feature_table.txt.gz" in filename:
            accession_key = filename.split("_feature_table.txt.gz")[0]
            if accession_key not in file_pairs:
                file_pairs[accession_key] = {"protein": [], "feature": []}
            file_pairs[accession_key]["feature"].append(file_path)
            
    return file_pairs

def process_pairs(pairs, original_script):
    """
    Iterates through the found pairs, decompresses them, runs the original
    script, and cleans up intermediate files.
    """
    for accession, files in pairs.items():
        if len(files["protein"]) == 1 and len(files["feature"]) == 1:
            protein_file_gz = files["protein"][0]
            feature_file_gz = files["feature"][0]

            print(f"Processing pair for accession: {accession}")
            
            # Decompress the files
            decompressed_protein = protein_file_gz[:-3]
            decompressed_feature = feature_file_gz[:-3]

            try:
                print(f"  Decompressing {os.path.basename(protein_file_gz)}...")
                with gzip.open(protein_file_gz, 'rb') as f_in, open(decompressed_protein, 'wb') as f_out:
                    f_out.write(f_in.read())
                
                print(f"  Decompressing {os.path.basename(feature_file_gz)}...")
                with gzip.open(feature_file_gz, 'rb') as f_in, open(decompressed_feature, 'wb') as f_out:
                    f_out.write(f_in.read())
            except Exception as e:
                print(f"Error decompressing files for {accession}: {e}")
                continue

            # Define the output filename with the new format
            try:
                # Get the filename without the path and suffix
                prefix_parts = os.path.basename(protein_file_gz).split('-', 2)[:2]
                prefix = "_".join(prefix_parts)
                output_fasta = f"{prefix}_{accession}_modified.fasta"
            except IndexError:
                print(f"  Error: Could not parse prefix from {os.path.basename(protein_file_gz)}. Using default naming.")
                output_fasta = f"{accession}_modified.fasta"

            # Run the user's original script
            try:
                print(f"  Running '{original_script}' on decompressed files...")
                subprocess.run(
                    [sys.executable, original_script, decompressed_protein, decompressed_feature, output_fasta],
                    check=True,
                    capture_output=True,
                    text=True
                )
                print(f"  Successfully created output file: {output_fasta}")
            except subprocess.CalledProcessError as e:
                print(f"  Error running script for {accession}:")
                print(f"  STDOUT: {e.stdout}")
                print(f"  STDERR: {e.stderr}")
            except FileNotFoundError:
                print(f"  Error: The script '{original_script}' was not found. Please ensure it is in the same directory.")
                
            # Clean up intermediate files
            try:
                print("  Cleaning up intermediate files...")
                os.remove(decompressed_protein)
                os.remove(decompressed_feature)
            except OSError as e:
                print(f"  Error removing intermediate files: {e}")
        else:
            print(f"Skipping accession '{accession}': Found an incomplete pair.")

def main():
    """Main function to orchestrate the workflow."""
    current_directory = "."
    print("Starting the batch processing script...")
    
    pairs = find_file_pairs(current_directory)
    
    if not pairs:
        print("No file pairs found. Please check the directory and file names.")
        return
        
    print(f"Found {len(pairs)} potential file pairs.")
    
    process_pairs(pairs, "modify_fasta_headers.py")
    
    print("\nBatch processing complete.")

if __name__ == "__main__":
    main()
