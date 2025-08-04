# addGeneSymbol2RefSeq

This script updates the headers of a RefSeq FASTA file. It reads a tab-separated value (TSV) feature table file (available on NCBI where you got the fasta), mapping protein accession numbers to gene symbols. For each sequence in the input FASTA file, the script identifies the accession number in the header and searches for a matching entry in the TSV file. If a match is found, it appends the corresponding gene symbol to the end of the FASTA header in the format "GN=<symbol>", similar to how UniProt FASTA headers are formatted. The script then writes the sequences with their modified headers to a new output FASTA file.

Usage:
To use, run the following command from your terminal in the directory containing the files:
python modify_fasta_headers.py <input.fasta> <input.tsv> <output.fasta>

- <input.fasta>: The original FASTA file.
- <input.tsv>: The feature table TSV file (will contain 'product_accession' and 'symbol' columns).
- <output.fasta>: The name for the new FASTA file with updated headers.


Example files:
There are two ways to retrieve the needed files. First you can use [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets) to search for the desired species, in this example [Zalophus californianus](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9704/). When you click the download button you will have various file types to download, including the protein (FASTA). Now you would think this is where the feature table is, and one day maybe it will, but it is actually over in a different location. Now that you know the species name and that it has a protein FASTA, navigate to the [Eukaryotice genomes annotated page](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/all/), and find the species. There will be four links on the right column, FTP, B, AR, and GDV. Though the AR (annotation report) is great, select [FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9704/101). This will bring you to a directory structure that is for the current release (you will recognize this GCF "name" from the NCBI Datasets page. If you need to get a different release you can navigate up a directory, but in this case we will click on the directory with the name of the current [annotation](https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9704/101/GCF_009762305.2_mZalCal1.pri.v2/), and will see a feature_table.txt.gz file to download. This is the file we need. You could also download the protein.fasta here which is nice since it keeps names in the files, helping traceability. 

You should have these two files now, which are also included here:
GCF_009762305.2_mZalCal1.pri.v2_feature_table.txt
GCF_009762305.2_mZalCal1.pri.v2_protein.fasta