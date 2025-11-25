# fasta2UniProt Scripts

This folder contains various helper scripts with the general goal of making non-UniProt FASTA into quasi-UniProt FASTA. This is important for search tools such as FragPipe and DIA-NN.


## convert_Bat1K_fasta_headers.py

The first script, named **convert_Bat1K_fasta_headers.py**, was created to make an in-house one-off protein FASTA look like a UniProt FASTA. This situation is described in the file itself.


## spoof_uniprot_headers.py

The second script, **spoof_uniprot_headers.py**, is built off the [addGeneSymbol2RefSeq](https://github.com/neely/proteomicScripts/tree/main/addGeneSymbol2RefSeq) output. This script takes a FASTA file with headers in the format:
```
>XP_024421550.2 lactotransferrin [Desmodus rotundus] GN=LTF
```

and transforms them into a UniProt-like format:

```
>rs|XP_024421550.2|LTF_DESRO lactotransferrin OS=Desmodus rotundus GN=LTF
```

NOTE: I may have it placing 'rs' to indicate RefSeq and not confuse sequences with sp (swiss prot) or tr (trembl), but it turns out some parsers are actually looking for either sp or tr (like FragPipe) therefore change this >rs| to >sp| or just check in the output.


Overall, the **spoof_uniprot_headers.py** script extracts the necessary information (accession, gene symbol, and species name) to construct the new formatted header. It processes each input FASTA file and creates a new corresponding output file.

This script processes all FASTA files in the same directory where it is run.

#### HOW TO RUN spoof_uniprot_headers.py:
1. Place this script in a folder containing the FASTA file(s) you want to process.
2. Open your terminal or command prompt.
3. Navigate to that folder.
4. Run the script with the following command:
```python spoof_uniprot_headers.py```

<br>
The script will automatically find all files with a '.fasta' extension, process each one, and save the output as a new file with the suffix '_uniprot_spoof.fasta'.

<br>
## ensembl_spoof_uniprot_headers.py

The third script, **ensembl_spoof_uniprot_headers.py** is the same concept but takes in Ensembl formatted files. This is rare you would need to use an Ensembl file instead of a RefSeq (or even UniProt for that matter), so maybe don't worry about this unless you have hyrax data. The usage is the same as **spoof_uniprot_headers.py**.

This **ensembl_spoof_uniprot_headers.py** script takes an Ensembl-formatted FASTA file with headers like this:
```
>ENSPCAT00000000047.1 ... gene_symbol:PPP1R3B description:protein phosphatase...
```

And transforms them into a UniProt-like format:
```
>sp|ENSPCAT00000000047.1|PPP1R3B_PROCA protein phosphatase 1 regulatory subunit... OS=Procavia capensis GN=PPP1R3B
```

The script performs a two-pass process on each file:

Pass 1: Identifies unique protein accessions and handles missing data.
- If a 'description' is missing, it is set to "unreported".
- If a 'gene_symbol' is missing, a unique identifier like 'LOC1' is generated.
- If a protein accession is duplicated, a '-counter' is added to make it unique.

Pass 2: Reconstructs the new headers and writes the final FASTA file.

In this case, the species name is hardcoded as 'Procavia capensis'.