#!/usr/bin/env python

#	Arguments: commitA commitB threads merge

#Make the imports
import sys
import os
import pickle
import shutil
import subprocess
from subprocess import Popen
from StoreAlignment import load_variables,obtain_alignment_pickle,store_variables
from reconstruct_functions import reconstruct_dataset
from parse_functions import parse_delta_file,parse_snp_file,parse_coords_file
from ObtainAlignment_functions import obtain_alignment
from update_functions import detect_updates

#Load arguments: commit_A commit_B threads report
commit_A=str(sys.argv[1])
commit_B=str(sys.argv[2])
number_threads=str(sys.argv[3])
merge=int(sys.argv[4])

#Create a temporary directory
if(os.path.isdir("./.git/info/temporary_directory")):
	shutil.rmtree("./.git/info/temporary_directory")
os.mkdir("./.git/info/temporary_directory")
#Obtain the commit messages
print("\nNow producing a summary review between commits "+commit_A+" and "+commit_B+". Warning: If you have added any data into the repository but did not commit, these changes will be lost.\n")
message_A=subprocess.check_output("git log --format=%B -n 1 "+commit_A, shell=True)
message_B=subprocess.check_output("git log --format=%B -n 1 "+commit_B, shell=True)
#Reconstruct the first version of the assembly
ShellCommand=Popen("git checkout "+str(sys.argv[1])+" 2> /dev/null",shell=True).wait()
reconstruct_dataset(size=60,directory="./Genome",output_file="./.git/info/temporary_directory/assembly_1.fa",mode="Genome")
#Reconstruct the second version of the assembly
ShellCommand=Popen("git checkout "+str(sys.argv[2])+" 2> /dev/null",shell=True).wait()
reconstruct_dataset(size=60,directory="./Genome",output_file="./.git/info/temporary_directory/assembly_2.fa",mode="Genome")
#Checkout back to the master
ShellCommand=Popen("git checkout master 2> /dev/null",shell=True).wait()

#Obtain alignment pickle of both assemblies
alignment_pickle=obtain_alignment_pickle("./.git/info/temporary_directory/assembly_1.fa","./.git/info/temporary_directory/assembly_2.fa")
#If the alignemnt containing the information between these two assemblies already exists, parse the information in the alignment files
if(os.path.isdir(alignment_pickle)):
	#Inform the user
	print("A stored alignment has been found. Now loading data.\n")
	#Load the OldNewID dictionary
	variables=load_variables(alignment_pickle+"/pickle")
	tabix_queries=variables[0]
	OldNewID_Dict=variables[1]
	alignment_pickle=variables[2]
	summary_Dict=variables[3]
#Otherwise it is required to perform the alignment between both assemblies
else:
	#Inform the user
	print("No stored alignment was found. Now creating genome alignment. This will take some time.\n")
	#Detect the updates
	ToUpdate=detect_updates("./RepoMap.txt")
	#Obtain the variables
	variables=obtain_alignment(old_assembly="./.git/info/temporary_directory/assembly_1.fa",new_assembly="./.git/info/temporary_directory/assembly_2.fa",directory="./.git/info/temporary_directory",threads=number_threads,ToUpdate=ToUpdate,alignmentpickle=alignment_pickle)
	#Store the pickle [tabix_queries,OldNewID_Dict,alignment_pickle,summary_Dict,file_crack]
	store_variables(variables=variables,alignment_pickle=alignment_pickle)
	#Load the variables
	tabix_queries=variables[0]
	OldNewID_Dict=variables[1]
	alignment_pickle=variables[2]
	summary_Dict=variables[3]
	file_crack=variables[4]

#If the user wants information about the merged sequences
if(merge):
	#Create an empty merge dictionary {newID:[oldID,oldID]...}
	merge_dict={}
	#Initiate counts
	mergeA_count=0
	mergeB_count=0
	#Parse the delta file {oldID:[length_A,length_B,score]...}
	delta_dict=parse_delta_file(alignment_pickle+"/Filtered.delta")
	#Loop through the keys of the oldnew dictionary
	for oldID in OldNewID_Dict.keys():
		#If the newID is already a key in the merge dictionary, add it to the list
		if(OldNewID_Dict[oldID] in merge_dict.keys()):
			merge_dict[OldNewID_Dict[oldID]].append(oldID)
		#Otherwise create a new list
		else:
			merge_dict[OldNewID_Dict[oldID]]=[oldID]
	#Inform the user
	print("The --merge option has been detected: a merged sequences report will be produced.\n")
	#Open a report file
	report_file=open("../GenomeGit_Merge.txt","w")
	#Write info
	report_file.write("\nNow producing a merge sequences report between commits "+commit_A+" and "+commit_B+". Warning: If you have added any data into the repository but did not commit, these changes will be lost.\n")
	report_file.write("\n\t###\n\t#\tVERSION A: "+message_A.rstrip()+" ("+commit_A+")\n\t###\n")
	report_file.write("\n\t###\n\t#\tVERSION B: "+message_B.rstrip()+" ("+commit_B+")\n\t###\n\n")
	#Loop through the newIDs
	for newID in merge_dict.keys():
		#Only if the newID is a product of a merge sequence
		if(len(merge_dict[newID])>1):
			#Update the counts
			mergeA_count+=len(merge_dict[newID])
			mergeB_count+=1
			#Write the information in the file
			report_file.write("\t------\n\t|\n\t| Sequence "+newID+" ("+delta_dict[merge_dict[newID][0]][1]+" nt) is a result of "+str(len(merge_dict[newID]))+" merged sequences:")
			#Loop through the oldIDs
			for oldID in merge_dict[newID]:
				#Print the info
				report_file.write("\n\t|\n\t|\t- "+oldID+ "("+delta_dict[oldID][0]+" nt)")
			report_file.write("\n\t|\n\t------\n")
	#Write the end of the report
	report_file.write("\n***End of report.")
	#Close the file
	report_file.close()
	#Inform the user
	print("Review completed, merged sequences report stored in GenomeGit_Merge.txt file.\n"+str(mergeA_count)+" sequences present in "+message_A.rstrip()+" experienced merging.\n"+str(mergeB_count)+" sequences present in "+message_B.rstrip()+" are a product of two or more sequences merging.")
			
#Otherwise produce a normal report
else:
	#Parse the coords file {oldID:[contigs]...}
	coords_dict=parse_coords_file(alignment_pickle+"/summary.coords")
	#Parse the snp file	{oldID:[insertions,deletetions,substitutions]}
	snp_dict=parse_snp_file(alignment_pickle+"/Filtered.snp")
	#Parse the delta file {oldID:[length_A,length_B,score]...}
	delta_dict=parse_delta_file(alignment_pickle+"/Filtered.delta")
	#Create a table with the information about the sequences {oldID:[table lines]}
	table={}
	#Initiate empty variables
	identical_count=0
	deleted_count=0
	reversed_count=0
	snp_count=0
	#Add information regarding the deleted, reversed and identical sequences
	for oldID in summary_Dict.keys():
		#If it is identical sequence
		if(summary_Dict[oldID][0]=="identical"):
			table[oldID]=["\t------\n\t|\n\t| Sequence: "+oldID+" ("+summary_Dict[oldID][1]+" nt) --> "+OldNewID_Dict[oldID]+" ("+summary_Dict[oldID][1]+" nt) (identical sequences).\n\t|\n\t------"]
			identical_count+=1
		#If it is a deleted sequence
		elif(summary_Dict[oldID][0]=="deleted"):
			table[oldID]=["\t------\n\t|\n\t| Sequence: "+oldID+" ("+summary_Dict[oldID][1]+" nt) was deleted from the assemlby.\n\t|\n\t------"]
			deleted_count+=1
		#Otherwise it was a reversed sequence
		else:
			table[oldID]=["\t------\n\t|\n\t| Sequence: "+oldID+" ("+summary_Dict[oldID][1]+" nt) --> "+OldNewID_Dict[oldID]+" ("+summary_Dict[oldID][1]+" nt) (reversed sequences).\n\t|\n\t------"]
			reversed_count+=1
	#Loop through the keys in the delta_dict and add the corresponding rows to the table
	for oldID in delta_dict.keys():
		#Add the snps if there are any (if the oldID is a key in snp_dict)
		if oldID in snp_dict.keys():
			table[oldID]=["\t------\n\t|\n\t| Sequence: "+oldID+" ("+delta_dict[oldID][0]+" nt) --> "+OldNewID_Dict[oldID]+" ("+delta_dict[oldID][1]+" nt)"]
			table[oldID].append("\t| Score: "+delta_dict[oldID][2])
			table[oldID].append("\t| SNPs:\n\t|\n\t|\t-Insertions: "+str(snp_dict[oldID][0])+"\n\t|\t-Deletions: "+str(snp_dict[oldID][1])+"\n\t|\t-Substitutions: "+str(snp_dict[oldID][2])+"\n\t|")
			table[oldID].append("\t| Contig distribution:\n\t|")
			#Calculate SNPs
			snp_count+=snp_dict[oldID][0]+snp_dict[oldID][1]+snp_dict[oldID][2]
		#Otherwise complete the table with 0s
		else:
			table[oldID]=["\t------\n\t|\n\t| Sequence: "+oldID+" ("+delta_dict[oldID][0]+" nt) --> "+OldNewID_Dict[oldID]+" ("+delta_dict[oldID][1]+" nt)"]
			table[oldID].append("\t| Score: "+delta_dict[oldID][2])
			table[oldID].append("\t| SNPs:\n\t|\n\t|\t-Insertions: 0\n\t|\t-Deletions: 0\n\t|\t-Substitutions: 0\n\t|")
			table[oldID].append("\t| Contig distribution:\n\t|")
		#Add information regarding contigs. Initiate variables
		i=0
		coverage_quer=0
		coverage_ref=0
		#Loop through the contig dictionary
		for contig in coords_dict[oldID]:
			#Add the contig info to the table
			i+=1
			contig=contig.split("\t")
			table[oldID].append("\t|\t-Contig "+str(i)+": "+contig[0]+"-"+contig[1]+" --> "+contig[2]+"-"+contig[3]+"\t ("+contig[6]+"% sim.)")
			#Sum coverage
			coverage_ref=coverage_ref+float(contig[7])
			coverage_quer=coverage_quer+float(contig[8])
		#Append the closing of the table
		table[oldID].append("\t|\n\t| Total coverage:\n\t|\n\t|\t-Version A: "+str(coverage_ref)+"%\n\t|\t-Version B: "+str(coverage_quer)+"%")
		table[oldID].append("\t|\n\t------")

	#Print the table into the output file. 
	with open("../GenomeGit_Diff.txt","w") as output_file:
		#Print initial info
		output_file.write("\nNow producing a summary review between commits "+commit_A+" and "+commit_B+". Warning: If you have added any data into the repository but did not commit, these changes will be lost.\n")
		output_file.write("\n\t###\n\t#\tVERSION A: "+message_A.rstrip()+" ("+commit_A+")\n\t###\n")
		output_file.write("\n\t###\n\t#\tVERSION B: "+message_B.rstrip()+" ("+commit_B+")\n\t###\n")
		#Loop through the oldIDs
		for oldID in table.keys():
			#Print all the lines related to this sequence
			for line in table[oldID]:
				output_file.write("\n"+line)
			output_file.write("\n")
		output_file.write("\n\n***End of report.")
	#Close the file
	output_file.close()
	#Print a quick summary in the prompt
	print("Review completed, differences stored in GenomeGit_Diff.txt file.\nNumber of reversed sequences: "+str(reversed_count))
	print("Number of identical sequences: "+str(identical_count))
	print("Number of deleted sequences: "+str(deleted_count))
	print("Total number of SNPs detected: "+str(snp_count))

#Remove the temporary directory
shutil.rmtree("./.git/info/temporary_directory")
