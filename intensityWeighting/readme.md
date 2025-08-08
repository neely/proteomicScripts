# Intensity-Weighting Identifications


Typically, in bottom-up data-dependent mass spec proteomics data an identification rate is determined using the number of peptide spectral matches (PSMs), which can be divided by the total number of acquired spectra to give the identification rate. How many spectra are identified is a relationship between the sample prep (digestion and cleanup), the separation, the instruments mass accuracy, ion transfer efficiency, and fragmentation, and the settings used when searching the data. On this last point, this can include fasta choices, which even in humans can lead to 10 to 20 % change in identified PSMs.


The following scripts were developed to investigate if weighting the identified spectra using the precursor intensity or sum of MS2 ions. For instance, searching data generated from one crow species (Mariana Crow) using a fasta from another crow species (American Crow), will have lower identification rates than using a species-specific fasta. But peptide homology or diversity is not uniform across the proteome between species, and so we wondered how highly homologous and highly abundant proteins might show more nuanced effects than simply PSM identification rate. This is why these scripts were built.


These scripts use mgf files from msconvert, but any converter should work as long as their is ion intensity data in the document. For this example, msconvert was used. The mgf is used to get a scan number, precursor mass, precursor abundance (which may not be the actual precursor intensity), and fragment ion abundance. A more effective method would be to integrate the XIC for each precursor, but this was more computationally heavy and difficult than we wanted. The mzID used is from Mascot, but mzID from any search algorithm should work. This allowed assigning an "identification" to scan number.


Finally, there is an added "experimental" function to extract intensity of given fragment ion masses (targets.csv). This could be used for glycans, contaminants, and other diagnostic ions. This would allow quantification of how much of a given fragment ion is present, maybe indicating how many peptides are glycosylated (as an example). The attempt was to try and mimic the functionality of [mzsniffer](https://github.com/wfondrie/mzsniffer) that evaluates a contaminant list on the MS1-level, but do so on the MS2.

The usage of these files is:
process_mgf_mzid.py + mgf + optional targets csv + mzid = "csv of all MSMS scans with identified column and targets columns"
process_mgf.py + mgf + optional targets csv = "csv of all MSMS scans and targets columns"
output-mgfmzid-stats.py + "csv of all MSMS scans with identified column" = text summary of weighted and count identification rates

Ex. commands with included files
```
python process_mgf_mzid.py 2020-1-16_Crow-Brain-truncated.mgf targets.csv MarianaCrow-F003830.mzid MarianaCrow.csv
python output-mgfmzid-stats.py MarianaCrow.csv MarianaCrow.txt
python process_mgf.py 2020-1-16_Crow-Brain-truncated.mgf MarianaCrow-mgfonly.csv targets.csv
```
