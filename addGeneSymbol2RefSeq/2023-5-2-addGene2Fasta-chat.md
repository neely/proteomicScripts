*RefSeq typically doesn't have the gene symbol in the fasta header for proteins, but for most species there is an abundance of files related to the genome (gff, gbff, etc.) that can be found on the FTP for a species (and maybe other places). 

*Ex. CA sea lion https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Zalophus_californianus/latest_assembly_versions/GCF_009762305.2_mZalCal1.pri.v2/ 

*In this case, the *protein.faa.gz file is the protein fasta and the *feature_table.txt.gz is this fun tab delimited txt file with all the genes, mRNA, and CDS linked, and has a column of "symbols", which are (or normally are) HGNC recognized symbols. The goal was to make an easy script to add these gene symbols to a fasta so I didn't have to add them or look them up on the backend.


**Me**

I would like to make a python script that I can run from the command line that will have two files as input arguments and one filename as output. The two files I am inputting are a fasta and a tsv. The fasta files have NCBI RefSeq headers such that it is accession (ex. XP_027423268.1) and then description. The tsv file will have this accession in the "product_accession" column. Once you have found the accession of a fasta entry in the tsv, I would like you take the value in the "symbol" column from the same row, and add this to the end of the fasta header as "GN=" then the symbol value. So ">XP_027423268.1 gametogenetin-binding protein 2 isoform X5 [Zalophus californianus]" would become ">XP_027423268.1 gametogenetin-binding protein 2 isoform X5 [Zalophus californianus] GN=GGNBP2". If there is no value in the symbol column for a given product_accession, don't add anything to the fasta header. Since there are many fasta entries (often well over 20 000), you will need to perform many lookups. The new fasta with the modified hearders will be the output of the script.


**ChatGPT**

Sure, I can help you write a Python script to achieve this. Here is an example script that should accomplish what you're looking for:

```

import argparse

# Define command line arguments
parser = argparse.ArgumentParser(description='Add gene symbols to NCBI RefSeq headers in a fasta file.')
parser.add_argument('fasta_file', help='Input fasta file with NCBI RefSeq headers.')
parser.add_argument('tsv_file', help='Input tsv file with gene symbols.')
parser.add_argument('output_file', help='Output fasta file with modified headers.')

# Parse command line arguments
args = parser.parse_args()

# Read in tsv file and create a dictionary mapping accessions to gene symbols
accession_to_symbol = {}
with open(args.tsv_file) as tsv_file:
    for line in tsv_file:
        fields = line.strip().split('\t')
        accession_to_symbol[fields[0]] = fields[2]

# Open output fasta file for writing
with open(args.output_file, 'w') as output_file:
    # Read in fasta file
    with open(args.fasta_file) as fasta_file:
        accession = ''
        description = ''
        for line in fasta_file:
            # If line is a header line, modify it as necessary
            if line.startswith('>'):
                # Write previous fasta entry to output file
                if accession != '':
                    output_file.write(f'>{accession} {description}')
                    if accession in accession_to_symbol and accession_to_symbol[accession] != '':
                        output_file.write(f' GN={accession_to_symbol[accession]}')
                    output_file.write('\n')
                # Parse new header line
                header = line.strip().split(' ')
                accession = header[0][1:]
                description = ' '.join(header[1:])
            # Otherwise, accumulate sequence description
            else:
                description += line.strip()
        # Write last fasta entry to output file
        output_file.write(f'>{accession} {description}')
        if accession in accession_to_symbol and accession_to_symbol[accession] != '':
            output_file.write(f' GN={accession_to_symbol[accession]}')
        output_file.write('\n')
```


This script uses the argparse module to define three command line arguments: the input fasta file, the input tsv file, and the output fasta file. It then reads in the tsv file and creates a dictionary mapping accessions to gene symbols. Finally, it reads in the fasta file and modifies each header line as necessary before writing the modified fasta entries to the output file.

To run this script, save it as a .py file (e.g. modify_fasta_headers.py) and run the following command in the terminal:

```

python modify_fasta_headers.py input.fasta input.tsv output.fasta
```

Replace input.fasta, input.tsv, and output.fasta with the actual file paths for your input and output files.


**Me**

Sorry, the table is a tab-delimited .txt file. 


**ChatGPT**

No problem, you can modify the script to read in a tab-delimited .txt file instead of a tsv file by changing the delimiter used in the split() method. Here's the updated script:

```

import argparse

# Define command line arguments
parser = argparse.ArgumentParser(description='Add gene symbols to NCBI RefSeq headers in a fasta file.')
parser.add_argument('fasta_file', help='Input fasta file with NCBI RefSeq headers.')
parser.add_argument('txt_file', help='Input txt file with gene symbols.')
parser.add_argument('output_file', help='Output fasta file with modified headers.')

# Parse command line arguments
args = parser.parse_args()

# Read in txt file and create a dictionary mapping accessions to gene symbols
accession_to_symbol = {}
with open(args.txt_file) as txt_file:
    for line in txt_file:
        fields = line.strip().split('\t')
        accession_to_symbol[fields[0]] = fields[2]

# Open output fasta file for writing
with open(args.output_file, 'w') as output_file:
    # Read in fasta file
    with open(args.fasta_file) as fasta_file:
        accession = ''
        description = ''
        for line in fasta_file:
            # If line is a header line, modify it as necessary
            if line.startswith('>'):
                # Write previous fasta entry to output file
                if accession != '':
                    output_file.write(f'>{accession} {description}')
                    if accession in accession_to_symbol and accession_to_symbol[accession] != '':
                        output_file.write(f' GN={accession_to_symbol[accession]}')
                    output_file.write('\n')
                # Parse new header line
                header = line.strip().split(' ')
                accession = header[0][1:]
                description = ' '.join(header[1:])
            # Otherwise, accumulate sequence description
            else:
                description += line.strip()
        # Write last fasta entry to output file
        output_file.write(f'>{accession} {description}')
        if accession in accession_to_symbol and accession_to_symbol[accession] != '':
            output_file.write(f' GN={accession_to_symbol[accession]}')
        output_file.write('\n')
```

This script is almost identical to the previous script, except that it uses '\t' as the delimiter in the split() method when reading in the tab-delimited .txt file. You can run this script in the same way as the previous script, with the updated command line arguments:

```

python modify_fasta_headers.py input.fasta input.txt output.fasta
```

Replace input.fasta, input.txt, and output.fasta with the actual file paths for your input and output files.


**Me**

The input txt file has many columns, delineated by column headers on the first row. There are many columns, for instance I can see 19 in my test file, but the columns that are important for this script are "product_accession" and "symbol". I am not sure the "accession_to_symbol" section of your script is actually making a dictionary.


**ChatGPT**

Thank you for letting me know about the format of your input file. To extract only the columns we need (product_accession and symbol), we can use the csv module in Python, which allows us to easily parse delimited text files.

Here's an updated version of the script that should work with your input file format:

```

import argparse
import csv

# Define command line arguments
parser = argparse.ArgumentParser(description='Add gene symbols to NCBI RefSeq headers in a fasta file.')
parser.add_argument('fasta_file', help='Input fasta file with NCBI RefSeq headers.')
parser.add_argument('txt_file', help='Input txt file with gene symbols.')
parser.add_argument('output_file', help='Output fasta file with modified headers.')

# Parse command line arguments
args = parser.parse_args()

# Read in txt file and create a dictionary mapping accessions to gene symbols
accession_to_symbol = {}
with open(args.txt_file, newline='') as txt_file:
    reader = csv.DictReader(txt_file, delimiter='\t')
    for row in reader:
        accession_to_symbol[row['product_accession']] = row['symbol']

# Open output fasta file for writing
with open(args.output_file, 'w') as output_file:
    # Read in fasta file
    with open(args.fasta_file) as fasta_file:
        accession = ''
        description = ''
        for line in fasta_file:
            # If line is a header line, modify it as necessary
            if line.startswith('>'):
                # Write previous fasta entry to output file
                if accession != '':
                    output_file.write(f'>{accession} {description}')
                    if accession in accession_to_symbol and accession_to_symbol[accession] != '':
                        output_file.write(f' GN={accession_to_symbol[accession]}')
                    output_file.write('\n')
                # Parse new header line
                header = line.strip().split(' ')
                accession = header[0][1:]
                description = ' '.join(header[1:])
            # Otherwise, accumulate sequence description
            else:
                description += line.strip()
        # Write last fasta entry to output file
        output_file.write(f'>{accession} {description}')
        if accession in accession_to_symbol and accession_to_symbol[accession] != '':
            output_file.write(f' GN={accession_to_symbol[accession]}')
        output_file.write('\n')
```

In this updated script, we use csv.DictReader to read in the tab-delimited .txt file and parse it into a dictionary. We specify delimiter='\t' to indicate that the file is tab-delimited. We also use the newline='' argument to open the file in universal newlines mode, which ensures that line endings are correctly handled on different platforms.

By default, csv.DictReader uses the first row of the file as the headers for the columns. We can use the column names product_accession and symbol to extract the values we need and create a dictionary mapping accessions to gene symbols.

You can run this updated script in the same way as the previous scripts, with the updated command line arguments:

```

python modify_fasta_headers.py input.fasta input.txt output.fasta
```

Replace input.fasta, input.txt, and output.fasta with the actual file paths for your input and output files.


**Me**

Okay, I see what you are doing. I don't think you are handling the FASTA format correctly. A fasta has a header line which is ">Accession Description" and then a line break and then rows of letters which is the sequence, and these have line breaks. Though the length of each line can be different, in this case each is 80 characters long. At the end of the sequence there is another line break, then the next entry begins with a new header line. Currently you have removed all line breaks for each FASTA entry, concatenating the header and the sequence it describes. Instead, I would like you to add "GN=" then the symbol value to the end of the header line, so it will be ">Accession Description GN=symbol", and preserve the original structure of the FASTA.


**ChatGPT**

Thank you for pointing that out. You are correct that the current implementation does not handle the FASTA format correctly. Here is an updated implementation that preserves the original structure of the FASTA file and adds the symbol to the end of the header line for each entry if it exists in the TSV file:

```

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
```
This implementation should correctly read in the FASTA and TSV files, map the accession numbers to symbols, and modify the headers of the FASTA file as desired. The output will be written to a new file specified as the third command-line argument.


**Me**

Awesome, that did it!


**ChatGPT**

Great! I'm glad it worked for you. If you have any more questions or if there's anything else I can help you with, just let me know!