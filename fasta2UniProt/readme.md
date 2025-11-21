# fasta2UniProt Scripts

This folder contains various helper scripts with the general goal of making non-UniProt FASTA into quasi-UniProt FASTA. This is important for search tools such as FragPipe and DIA-NN.

The first script, named generically convert_fasta_headers.py, was created to make an in-house one-off protein FASTA look like a UniProt FASTA. This situation is described in the file itself.

The second script, spoof_uniprot_headers.py, is built off the [addGeneSymbol2RefSeq](https://github.com/neely/proteomicScripts/tree/main/addGeneSymbol2RefSeq) output. This script takes a FASTA file with headers in the format:
```
>XP_024421550.2 lactotransferrin [Desmodus rotundus] GN=LTF
```

and transforms them into a UniProt-like format:

```
>rs|XP_024421550.2|LTF_DESRO lactotransferrin OS=Desmodus rotundus GN=LTF
```

NOTE: I have it placing 'rs' to indicate RefSeq and not confuse sequences with sp (swiss prot) or tr (trembl), but it turns out some parsers are actually looking for either sp or tr (like DIA-NN) therefore change this >rs| to >sp| .


Overall, the script extracts the necessary information (accession, gene symbol, and species name) to construct the new formatted header. It processes each input FASTA file and creates a new corresponding output file.

This script processes all FASTA files in the same directory where it is run.

#### HOW TO RUN:
1. Place this script in a folder containing the FASTA file(s) you want to process.
2. Open your terminal or command prompt.
3. Navigate to that folder.
4. Run the script with the following command:
```python spoof_uniprot_headers.py```

<br>
The script will automatically find all files with a '.fasta' extension, process each one, and save the output as a new file with the suffix '_uniprot_spoof.fasta'.