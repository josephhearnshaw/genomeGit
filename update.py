#!/usr/bin/env python

#	INSTRUCTIONS: The script will update the coordinates of the variants and annotation files of a given genome assembly to match those of a new one. Arguments:
#	update_dependent_datasets.py -new_file- -threads- -file_size- -template_length-

#Make the imports
import os
import sys
import datetime
from subprocess import Popen
from update_functions import detect_updates
from reconstruct_functions import reconstruct_dataset
from ObtainAlignment_functions import obtain_alignment
from StoreAlignment import store_variables,load_variables,obtain_alignment_pickle
from InterpretAlignment_functions import interpret_alignment

#Load arguments: new_assembly threads size tlength
new_assembly=str(sys.argv[1])
number_threads=int(sys.argv[2])
file_size=str(sys.argv[3])
template_length=int(sys.argv[4])

##########
#	PART 0.DETERMINE FILES THAT NEED TO BE UPDATED
##########

#Determine the files to be updated. ToUpdate={dataset:[[filename.extension,directory,size],[],...]}
ToUpdate=detect_updates("./RepoMap.txt")
#If the there are no files to update, inform the user
if (ToUpdate=="empty"):
	print("No datasets to be updated were detected.\n")
	sys.exit()
#Inform the user of the files to be updated
else:
	print("Number of threads selected: "+str(number_threads)+"\nMaximum insert size: "+str(template_length)+"\nThe following files that need to be updated have been detected:\n")
	for dataset in ToUpdate.keys():
		if(dataset!="Genome"):
			print("\t-Files detected in the "+dataset+" dataset\n")
			#If there are files, print them
			if(len(ToUpdate[dataset])!=0):
				for subfile in ToUpdate[dataset]:
					#If the size of the file is lower than 1MB
					if(subfile[2]=="1"):
						print("\t\t--"+subfile[0]+" (<1 MB)\n")
					#Otherwise is bigger than 1MB
					else:
						print("\t\t--"+subfile[0]+" ("+subfile[2]+" MB)\n")
			#If there are no files in this dataset, inform the user
			else:
				print("\t\t--No files detected in the "+dataset+" dataset.\n")

#Inform the user that the update will start now
if(ToUpdate["Genome"][0][2]=="1"):
	if(file_size=="1"):
		print("\n*** NOW STARTING UPDATE OF THE REPOSITORY: "+ToUpdate["Genome"][0][0]+" (<1 MB) ---> "+os.path.basename(new_assembly)+" (<1 MB)***\n")
	else:
		print("\n*** NOW STARTING UPDATE OF THE REPOSITORY: "+ToUpdate["Genome"][0][0]+" (<1 MB) ---> "+os.path.basename(new_assembly)+" ("+file_size+" MB)***\n")
else:
	if(file_size=="1"):
		print("\n*** NOW STARTING UPDATE OF THE REPOSITORY: "+ToUpdate["Genome"][0][0]+" ("+ToUpdate["Genome"][0][2]+" MB) ---> "+os.path.basename(new_assembly)+" (<1 MB)***\n")
	else:
		print("\n*** NOW STARTING UPDATE OF THE REPOSITORY: "+ToUpdate["Genome"][0][0]+" ("+ToUpdate["Genome"][0][2]+" MB) ---> "+os.path.basename(new_assembly)+" ("+file_size+" MB)***\n")

##########
#	PART 1. RECONSTRUCTION OF REPOSITORY DATA
##########

#Inform the user
print("\n\t*PART I. RECONSTRUCTION OF REPOSITORY DATA.*")
print("\t"+str(datetime.datetime.now()))
#Create temporal directory
os.mkdir("./temporary_directory")
#Reconstruct the necessary datasets. First the old genome.
reconstruct_dataset(size=60,directory="./Genome",output_file="./temporary_directory/genome_old.fa",mode="Genome",seqID="0",region="0")
#Reconstruct all the files related to the variants and annotation datasets
for dataset in ToUpdate.keys():
	if(dataset!="Genome" and len(ToUpdate[dataset])!=0):
		#Loop through the files of the dataset and reconstruct the file in the temporary directory with the same original name
		for subfile in ToUpdate[dataset]:
			reconstruct_dataset(size=1,directory=subfile[1],output_file="./temporary_directory/"+subfile[0],mode=dataset,update=True,seqID="0",region="0")
			#Create an appropiate tabix library for the file. If it is annotation or alignment there are two pseudo files
			if(dataset=="Annotation" or dataset=="Alignment"):
				ShellCommand=Popen('(grep "^#" ./temporary_directory/'+subfile[0]+'_A; grep -v "^#" ./temporary_directory/'+subfile[0]+'_A | sort -V -k1,1 -k2,2n) | bgzip > ./temporary_directory/'+subfile[0]+'_A.gz; tabix -p vcf ./temporary_directory/'+subfile[0]+'_A.gz;',shell=True).wait()
				ShellCommand=Popen('(grep "^#" ./temporary_directory/'+subfile[0]+'_B; grep -v "^#" ./temporary_directory/'+subfile[0]+'_B | sort -V -k1,1 -k2,2n) | bgzip > ./temporary_directory/'+subfile[0]+'_B.gz; tabix -p vcf ./temporary_directory/'+subfile[0]+'_B.gz;',shell=True).wait()
			elif(dataset=="Variants"):
				ShellCommand=Popen('(grep "^#" ./temporary_directory/'+subfile[0]+'_A; grep -v "^#" ./temporary_directory/'+subfile[0]+'_A | sort -V -k1,1 -k2,2n) | bgzip > ./temporary_directory/'+subfile[0]+'_A.gz; tabix -p vcf ./temporary_directory/'+subfile[0]+'_A.gz;',shell=True).wait()
			#Otherwise there was an error
			else:
				print("***INTERNAL ERROR*** DATASET NOT RECOGNIZED: "+dataset)

##########
#	PART 2. OBTAIN THE ALIGNMENT INFORMATION
##########

# Obtain the information of the alignment between both assemblies.
# tabix queries:	tabix ./temporary_directory/filename.ext_AB.gz chr:x-y > ./temporary_directory/filename.ext_chr:x_y_AB
#
#		{tabix_query:["compare","_A/_B//","Annotation/Variants/Alignment",filename.ext,oldID,query_outfile,sub_updated,sub_discarded,[old_block_start,old_block_stop,new_block_start,new_block_stop,block_modifications],[snps...]]}
#		{tabix_query:["reversed","_A/_B//","Annotation/Variants/Alignment",filename.ext,oldID,query_outfile,sub_updated,sub_discarded,length]}
#		{tabix_query:"omitted"}
#		{tabix_query:"identical"}
#
# OldNewID_Dict	{oldID:newID,...}
# file_crack {finalfilename:[sub_updated/discarded...]}
#
#		finalfilename:	updated_filename.ext	/	discarded_filename.ext
#		sub_updated/discarded:	updated_filename.ext_chr:x_y_AB	/	discarded_filename.ext_chr:x_y_AB
#
# summary_Dict {oldID:[status,length]}

#Determine the alignment pickle
alignment_pickle=obtain_alignment_pickle("./temporary_directory/genome_old.fa",new_assembly)
#If there is an alignment already stored in the repository
if(os.path.isdir(alignment_pickle)):
	#Inform the user
	print("\n\t*PART II. OBTAINING GENOME ALIGNMENT: STORED ALIGNMENT DETECTED, NOW LOADING SAVED DATA.*")
	print("\t"+str(datetime.datetime.now()))
	#Load the variables stored in the pickle [tabix_queries,OldNewID_Dict,alignment_pickle,summary_Dict,file_crack]
	variables=load_variables(alignment_pickle+"/pickle")
	tabix_queries=variables[0]
	OldNewID_Dict=variables[1]
	alignment_pickle=variables[2]
	summary_Dict=variables[3]
	file_crack=variables[4]
#Otherwise it is necessary to obtain the alignment
else:
	#Inform the user
	print("\n\t*PART II. OBTAINING GENOME ALIGNMENT: NO STORED ALIGNMENT DETECTED, NOW CREATING NEW ALIGNMENT*")
	print("\t"+str(datetime.datetime.now()))
	#Obtain the variables
	variables=obtain_alignment(old_assembly="./temporary_directory/genome_old.fa",new_assembly=new_assembly,directory="./temporary_directory",threads=number_threads,ToUpdate=ToUpdate,alignment_pickle=alignment_pickle)
	#Store the pickle [tabix_queries,OldNewID_Dict,alignment_pickle,summary_Dict,file_crack]
	tabix_queries=variables[0]
	OldNewID_Dict=variables[1]
	alignment_pickle=variables[2]
	summary_Dict=variables[3]
	file_crack=variables[4]
	store_variables(variables=variables,alignment_pickle=alignment_pickle)

##########
#	PART 3. START THE INTERPRETATION OF THE ALIGNMENT INFORMATION
##########

#Inform the user
print("\n\t*PART III. INTERPRETATION OF THE ALIGNMENT INFORMATION AND CREATION OF UPDATED FILES")
print("\t"+str(datetime.datetime.now()))
#Interpret the information contained in the delta_dict and obtain the updated files.
interpret_alignment(queries=tabix_queries,oldnew=OldNewID_Dict,threads=number_threads,ToUpdate=ToUpdate,tlength=template_length,filecrack=file_crack)

#Inform the user the update is completed
print("\n\t***UPDATE COMPLETED: NOW PARSING THE GENOME DATASET***")
