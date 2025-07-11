# Uses peptide level report from Spectronaut that includes PEP.StrippedSequence values. The script will
# generate simple txt files from these tsv's, which can be used in other peptide tools here.
# usage is 
# python process_tsv.py 
# from the directory of the tsvs

import os
import pandas as pd

def process_files(directory):
    report_lines = []

    for filename in os.listdir(directory):
        if filename.endswith(".tsv"):
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath, sep='\t')

            if 'PEP.StrippedSequence' in df.columns:
                unique_sequences = df['PEP.StrippedSequence'].drop_duplicates()
                unique_count = unique_sequences.shape[0]

                # Write unique sequences to a new txt file
                output_filepath = os.path.join(directory, filename.replace('.tsv', '.txt'))
                with open(output_filepath, 'w') as f:
                    f.write('PEP.StrippedSequence\n')
                    for sequence in unique_sequences:
                        f.write(sequence + '\n')

                # Add to report
                report_lines.append(f"{filename}: {unique_count} unique sequences")

    # Write the report
    report_filepath = os.path.join(directory, "report.txt")
    with open(report_filepath, 'w') as report_file:
        report_file.write('\n'.join(report_lines))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process TSV files to extract unique PEP.StrippedSequence values.")
    parser.add_argument("directory", nargs='?', default='.', help="Directory containing the TSV files")
    args = parser.parse_args()

    process_files(args.directory)
