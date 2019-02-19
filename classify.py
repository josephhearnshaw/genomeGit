#!/usr/bin/env python

# PYTHON SCRIPT TO DETERMINE FILE FORMAT (FASTA, VCF, GFF)

#Make the imports
import sys

#Intiate variables:
comments=False
header=False

#Open the input file
with open(str(sys.argv[1]),"r") as input_file:
	#Read the first line of the file
	line=input_file.readline()
	#If the first line starts with >, it is a genome dataset
	if(line[0]==">"):
		print("Genome")
		sys.exit()
	#Otherwise it is necessary to loop through the lines of the file to determine if it is vcf or gff
	for line in input_file:
		#Only consider the line if it is not a comment
		if (line[0]!="#" and line[0]!="@"):
			#Split the line into its fields
			line=line.split()
			#If the number of fields is different than 0, determine which fields are numbers (field 2 in vcf; fields 4 and 5 in gff)
			if(len(line)!=0):
				#If it is alignment file
				if(len(line)>=11 and line[1].isdigit() and line[3].isdigit() and line[4].isdigit() and line[3].isdigit() and line[7].isdigit()):
					print("Alignment")
					sys.exit()
				#If it is annotation file
				elif(len(line)>=9 and line[3].isdigit() and line[4].isdigit()):
					print("Annotation")
					sys.exit()
				#If it is variants file
				elif(len(line)>=8 and line[1].isdigit()):
					print("Variants")
					sys.exit()
				#Otherwise it is an unrecognised format
				else:
					print("*** ERROR. UNABLE TO RECOGNISE THE FORMAT OF THE PROVIDED FILE. ONLY FASTA, GFF, VCF OR SAM FILES WILL BE ACCEPTED. ***")
					sys.exit()
			#Otherwise it is an unrecognised dataset
			else:
				print("*** ERROR. UNABLE TO RECOGNISE THE FORMAT OF THE PROVIDED FILE. ONLY FASTA, GFF, VCF OR SAM FILES WILL BE ACCEPTED. ***")
				sys.exit()
	#If the loop ends, it means that the file was not classified properly (is not fasta, gff or vcf)
	print("*** ERROR. UNABLE TO RECOGNISE THE FORMAT OF THE PROVIDED FILE. ONLY FASTA, GFF OR VCF FILES WILL BE ACCEPTED. ***")
