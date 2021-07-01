#!/bin/python
###################
## This python script takes a .fasta input file containing a reference
## genome sequence, and generates a basic .vcf file containing all 
## possible single nucleotide variations possible. 
## This .vcf file can be used as input file for snpeff, to assess the 
## relative frequencies of theoretically possible mutations. 
## required input arguments: 
##  - file containing minimum header 
##    (lines ##fileformat, ##fileDate, ##reference, list of contigs, ##INFO)
##  - .fasta file with reference sequence (scaffold names start with '>'
##  - output filename
###################
## Cas Retel, cas.retel@eawag.ch, 2020.01.30
## ! Does not do any checks on validity of the reference .fasta file ! 
## ! Has not been tested for reference genomes consisting of more
## than one scaffold !
###################
from sys import argv
import re 

## process input arguments
script, file_header, file_ref, file_out = argv

print "\n###### CREATING AN ARTIFICAL .VCF FILE ######"
print "--- Inputs: %s and %s ---" %(file_header, file_ref)
print "--- Writing to: %s ---\n" %(file_out)

## write lines in file_header to file_out
obj_header = open(file_header, 'r')
obj_out = open(file_out, 'w')
for line in obj_header:
	obj_out.write(line)

obj_header.close()

## create a new line containing all necessary information, append to file_out
scaffold = ""
bases = ['A', 'C', 'G', 'T']
ref = open(file_ref, 'r')
for line in ref: 
	## if it starts with '>', update scaffold name
	if re.search('^>', line): 
		scaffold = line.lstrip('>').rstrip('\n')
		print "--- Starting with %s ---" %(scaffold)
	## otherwise, loop over every character (base) in the line (sequence), 
	## select the three bases it is not equal to, write those to a file
	else:
		pos = 1
		for refbase in line.rstrip('\n'):
			# print "--- next base is: %s ---" %(refbase)
			for otherbase in bases: 
				if otherbase != refbase: 
					appendstring = "%s\t%s\t.\t%s\t%s\t10000\t.\tTYPE=snp;technology.illumina=1\tGT:DP:RO:AO\t0/1:200:100:100\n" %(scaffold, pos, refbase, otherbase)
					obj_out.write(appendstring)
			pos += 1
		
print "###### DONE ######"
ref.close()
obj_out.close()



# python /Users/reteladmin/Documents/HVInt/Project3_chlorovirusmuts/scripts/createVCFWithAllPossSNPs.py minimum_header.vcf ~/Documents/HVInt/Chlorella/ref/PBCV1.fasta pbcv1_allpossiblemuts.vcf