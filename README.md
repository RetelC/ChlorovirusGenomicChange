##############################
## Scripts belonging to manuscript: 
## Strong selection and high mutation supply characterise 
## experimental Chlorovirus evolution
##############################
## Name: Cas Retel
## Last update: 2021.08.25
## Contact: casretel@gmail.com  |  philine.feulner@eawag.ch
##############################

The python script below takes a .fasta input file containing a reference
genome sequence, and generates a basic .vcf file containing all 
possible single nucleotide variations possible: 
- createVCFWithAllPossSNPs.py

The files below contain R functions, should be available and 
invoked using R::source() in order to run the R analysis scripts: 
- processFbsVCF.R
- plotAfs.R
- gg_multiplot.R
- calcGenDiv.R
- calcSEffective.R

The R analysis script below reads in a .vcf file as returned by freebayes, 
and applies variant filtering steps as described in the manuscript: 
- virus2020-1_filter_vcfs.R

The R analysis script below reads in the filtered .vcf file as returned by 
virus2020-1_filter_vcfs.R, and executes all of the visualisation and 
analysis steps discussed in the manuscript (visualise genomic variation
along the genome and the distribution of predicted phenotypic effects; 
analyse the correlations between predicted phenotypic effect, 
repeatability and frequency change; calculate population genetic diversity; 
and analyse the variation in the three 'variation hotspots'): 
- virus2020-1_analyse_genvar.R 


