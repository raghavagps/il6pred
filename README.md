# il6pred: A method for predicting and desiging IL-6 inducing peptides.
Several methods have been developed for the prediction of the antigenic regions for subunit vaccines designing. Interleukin-6 (IL-6) is a rapidly produced proinflammatory cytokine generated as an immune response in various infections and tissue injuries. The elevated level of IL-6 causes cytokine release syndrome(CRS) in severe COVID-19 patients. Based on existing knowledge, we develop an in silico tool that allows the user to predict, scan, and map the IL-6 inducing/non-inducing peptides.
## Reference
Dhall A, Patiyal S, Sharma N, Usmani SS, Raghava GPS (2020) Computer-aided prediction and design of IL-6 inducing peptides: IL-6 plays a crucial role in COVID-19. <a href="https://pubmed.ncbi.nlm.nih.gov/33034338/"> Brief Bioinform. 22(2):936-945 </a>.

## Web Server
https://webs.iiitd.edu.in/raghava/il6pred/
# Installation
```
git clone https://github.com/raghavagps/il6pred
change dir il6pred
unzip RF_model.zip
python il6pred.py -i peptide.fa
```
# Introduction
IL6Pred is developed for predicting, desiging and scanning interleukin-6 inducing peptides. More information on IL6Pred is abvailble from its web server http://webs.iiitd.edu.in/raghava/il6pred . This page provide information about stnadalone version of IL6Pred. Please read/cite following paper for complete information including algorithm behind IL6Pred.

**Models:** In this program, one model has been incorporated for predicting interleukin-6 inducing peptides. The model is trained on IL-6 inducing and non-inducing peptides.

**Modules/Jobs:** This program implement three modules (job types); i) Predict: for predictin interleukin-6 inducing peptides, ii) Design: for generating all possible peptides and computing interleukin-6 inducing potential (score) of peptides, iii) Scan: for creating all possible overlapping peptides of given length (window) and computing interleukin-6 inducing potential (score) of these overlapping peptides.

**Minimum USAGE:** Minimum usage is "python il6pred.py -i peptide.fa" where peptide.fa is a input fasta file. This will predict the interleukin-6 inducing potential of sequence  in fasta format. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma seperated variables).

**Full Usage:** Following is complete list of all options, you may get these options by "python il6pred.py -h" 
```
il6pred.py [-h] -i INPUT [-o OUTPUT] [-j {1,2,3}] [-t THRESHOLD]
                  [-w {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}]
                  [-d {1,2}]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -j {1,2,3}, --job {1,2,3}
                        Job Type: 1:predict, 2:design and 3:scan, by default 1
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.11
  -w {5,6,7,..,25}, --winleng
                        Window Length: 5 to 25 (scan mode only), by default 10
  -d {1,2}, --display {1,2}
                        Display: 1:Interleukin-6 inducing peptide, 2: All peptides, by defaIL6Pred
```

**Input File:** It allow users to provide input in two format; i) FASTA format (standard) and ii) Simple Format. In case of simple format, file should have one one peptide sequence in a single line in single letter code (eg. peptide.seq). Please note in case of predict and design module (job) length of peptide should be upto 25 amino acids, if more than 25, program will take first 25 residues. In case of scan module, minimum length of protein/peptide sequence should be more than equal to window length (pattern), see peptide.fa . Please note program will ignore peptides having length less than 8 residues (e.g., protein.fa).

**Output File:** Program will save result in CSV format, in case user do not provide output file name, it will be stored in outfile.csv.

**Threshold:** User should provide threshold between 0 and 1, please note score is propotional to interleukin-6 inducing potential of peptide.


IL6Pred Package Files
=======================
It contantain following files, brief descript of these files given below
```
INSTALLATION  	: Installations instructions

LICENSE       	: License information

README.md     	: This file provide information about this package

RF_model : Model file required for running Model 2

il6pred.py 	: Main python program 

outfile.csv	: Example output file in csv format

peptide.fa	: Example file contain peptide sequenaces in FASTA format

peptide.seq	: Example file contain peptide sequenaces in simple format

protein.fa	: Example file contain protein sequenaces in FASTA format 

Data            : This folder contains the files required to run the in-built python scripts.
```
# Address for contact
In case of any query please contact
```
Prof. G. P. S. Raghava, Head Department of Computational Biology,            
Indraprastha Institute of Information Technology (IIIT), 
Okhla Phase III, New Delhi 110020 ; Phone:+91-11-26907444; 
Email: raghava@iiitd.ac.in  Web: http://webs.iiitd.edu.in/raghava/
```
