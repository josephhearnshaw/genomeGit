#!/usr/bin/env python

# PYTHON FUNCTIONS REQUIRED IN THE REPOSITORY AUTO-UPDATES PART I: ALIGNMENT OBTENTION.

#make the imports
import re
import hashlib
import os
import datetime
import sys
from subprocess import Popen
from threading import Thread
from multiprocessing import Process
from StoreAlignment import obtain_alignment_pickle

#Create global variables. They will be used for the threads
processed=0
processing=0
tabix_queries={}
alignment_pickle=""

#Create a obtain_SHA1 function that that returns a SHA-1 hash corresponding with the input string.
def obtain_SHA1(input_string):
	#Calculate the SHA-1 hash of both input sequences
	return str(hashlib.sha1(input_string).hexdigest())

#Create a function to start a given function and its arguments in a subprocess in the shell. It will need to update the number of processes being processed and the number already processed
def call_batch(function,argument_list):
	global processing
	global processed
	p = Process(target=function, args=(argument_list,))
	p.start()
	p.join()
	processing-=1
	processed+=1

#Create a function parse_assembly to obtain a dictionary with the following characteristics: {seqKey:[oldID,[line,line...],hash,hash-1]}
def parse_assembly(assembly_file):
	#Create output dictionary
	assembly_dictionary={}
	#Open the asembly file
	with open(assembly_file, "r") as input_file:
		#Loop through the lines
		sequence_count=0
		for line in input_file:
			line=line.rstrip()
			#If the line starts with >, it is a new sequence
			if(line[0]==">"):
				#The key is the seqID
				current_key=line[1:]
				#Add the first two values (the seqID and an empty string for the sequence)
				assembly_dictionary[current_key]=[line[1:],[]]
			#Otherwise the line corresponds with the nucleotides of the sequence
			else:
				assembly_dictionary[current_key][1].append(line)
		#After looping through all the lines of the file, determine the SHA-1 of the sequences in the dictionary
		for key in assembly_dictionary:
			#Forward SHA1
			assembly_dictionary[key].append(obtain_SHA1("".join(assembly_dictionary[key][1])))
			#Reverse SHA1
			assembly_dictionary[key].append(obtain_SHA1("".join(assembly_dictionary[key][1][::-1])))
	#Close the file
	input_file.close()
	#Return the assembly dictionary
	return assembly_dictionary

#Create a parse snp file to store the information contained in a mummer snp file in form of a dictionary {SeqID:"oldIndex oldChar newChar newIndex","...",}
def parse_snp_file(file_path):
	#Initate variables
	snp_dict={}
	#Open the file and loop through the lines
	with open(file_path,"r") as snp_file:
		for line in snp_file:
			#Split the line fields
			line=line.split("\t")
			#If the seqID is already a key in the dictionary, append the new line to the list
			if(line[10] in snp_dict.keys()):
				snp_dict[line[10]].append(" ".join([line[0],line[1],line[2],line[3],line[10]]))
			#Otherwise create a new list
			else:
				snp_dict[line[10]]=[" ".join(line[0:4])]
	#The snp file can be now closed and deleted
	snp_file.close()
	return snp_dict

#Create a calculate score function to return the score of a given alignment
def calculate_score(list_BlocksModifications,lengthA,lengthB):
	#Initiate points and fragment_number variables
	coherence=0
	block_number=1
	#Calculate length factor
	length_factor=2-(float(lengthA)/float(lengthB)+float(lengthB)/float(lengthA))
	#Save the first block as the prevblock
	prevblock=list_BlocksModifications[0].split()
	#Loop through the blocks and modifications
	for block in list_BlocksModifications[1:]:
		#Make sure it is a block and not a modification index
		if(" " in block):
			block=block.split()
			#Add one to the count of fragments
			block_number+=1
			if(int(block[0])>int(prevblock[0]) and int(block[2])>int(prevblock[2])):
				coherence+=1
			elif(int(block[0])>int(prevblock[0]) and int(block[2])<int(prevblock[2])):
				coherence-=1
	#Determine the score of this alignment
	score=float(coherence)/float(block_number)+length_factor
	#Return the obtained score
	return score

#Create a function filter_delta to filter the results in the multi-delta so that it only contains the best alignment of each query sequence
def filter_delta(delta_file,out_file):
	#Initiate variables
	Unfiltered_DeltaDict={}	#{>alignment:[blocks and modifications]}
	Filtered_DeltaDict={}	#{old_seqID:[>alignment,[blocks and modifications],score]}
	#Create a delta dictionary out of the contents of the unfiltered delta file 
	with open(delta_file,"r") as input_file:
		#The first line corresponds with the files compared. It must be preserved
		first_line=input_file.readline().rstrip()
		second_line=input_file.readline().rstrip()
		#Loop through the lines
		for line in input_file:
			line=line.rstrip()	
			#if the line starts with >, it is a new alignment
			if (line[0]==">"):
				current_key=line
				Unfiltered_DeltaDict[current_key]=[]
			#Otherwise, include it corresponds with a block of the current alignment or a modification index (to be included as well).
			else:
				Unfiltered_DeltaDict[current_key].append(line)
	#Close the delta file
	input_file.close()
	#Loop through the alignments in the dictionary and delete those that do not pass a score treshold. Preserve only the best of each query.
	#score_treshold=-0.5	TRESHOLD NO LONGER ON USE
	for alignment in Unfiltered_DeltaDict.keys():
		#Split the alignment
		alignment=alignment.split()
		#calculate the score of the current alginment
		score=calculate_score(list_BlocksModifications=Unfiltered_DeltaDict[" ".join(alignment)],lengthB=alignment[2],lengthA=alignment[3])
		#print("OldID: "+alignment[0]+" NewID: "+alignment[1])	#THIS IS FOR TESTING
		#print(score)	#THIS IS FOR TESTING
		#Only consider the alignment if the score is higher than the threshold. TRESHOLD NO LONGER ON USE.
		#if (float(score)>=float(score_treshold)):
		#If the current alignment is better than the one already stored, substitute it
		try:
			if(score>Filtered_DeltaDict[alignment[0]][2]):
				Filtered_DeltaDict[alignment[0]]=[" ".join(alignment),Unfiltered_DeltaDict[" ".join(alignment)],score]
		#If it is the first alignment of this old seqID, add it to the filtered dictionary
		except KeyError:
			Filtered_DeltaDict[alignment[0]]=[" ".join(alignment),Unfiltered_DeltaDict[" ".join(alignment)],score]
	#Write the contents of Filtered_delta into a new Filtered delta file
	with open(out_file,"w") as output_file:
		#The first line corresponds with the first line of the unfiltered delta
		output_file.write(first_line+"\n")
		output_file.write(second_line+"\n")
		#Loop through the keys in the Filtered_DeltaDict
		for old_seqID in Filtered_DeltaDict.keys():
			#Write the alignment header. Add an extra field at the end indicating the score.
			output_file.write(Filtered_DeltaDict[old_seqID][0]+" "+str(Filtered_DeltaDict[old_seqID][2])+"\n")
			#Write the blocks and modifications
			for BlockModifications in Filtered_DeltaDict[old_seqID][1]:
				output_file.write(BlockModifications+"\n")
	#Close the delta file
	output_file.close()

#Create a obtain_alignemnt function to perform the alignemtns for the given assemblies.
def obtain_alignment(old_assembly,new_assembly,directory,threads,ToUpdate,alignment_pickle):
	#Create the pickle directory
	os.makedirs(alignment_pickle)
	#Initiate variables:
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
	tabix_queries={}
	OldNewID_Dict={}
	summary_Dict={}
	file_crack={}
	#Create an entry in the file_crack for each subfile in ToUpdate
	for dataset in ToUpdate.keys():
		#If the dataset is variants
		if(dataset=="Variants"):
			#Loop through the subfiles of the dataset and create a new key for each subfile with an empty list as a value. 
			for subfile in ToUpdate[dataset]:
				file_crack["./temporary_directory/updated_"+subfile[0]+"_A"]=[]
		#Only if the dataset is annotation or variants, do it for B as well
		elif(dataset!="Genome"):
			for subfile in ToUpdate[dataset]:
				file_crack["./temporary_directory/updated_"+subfile[0]+"_A"]=[]
				file_crack["./temporary_directory/updated_"+subfile[0]+"_B"]=[]
	#Parse the genome assemblies into dictionaries with the following characteristics: {seqKey:[oldID,[line,line...],hash,hash-1]}
	OldAssembly_Dict=parse_assembly(old_assembly)
	NewAssembly_Dict=parse_assembly(new_assembly)
	#All vs all alignment must be carried out, and then filtering of the resulting delta file
	#Loop through the keys in the OldAssembly_Dict
	for oldID in OldAssembly_Dict.keys():
		#Loop through the keys in NewAssembly_Dict
		for newID in NewAssembly_Dict.keys():
			#Check if the SHA1 is identical in the NewAssembly_Dict
			if(OldAssembly_Dict[oldID][2]==NewAssembly_Dict[newID][2]):
				#Add the info to the summary_Dict
				summary_Dict[oldID]=["identical",str(len("".join(OldAssembly_Dict[oldID][1])))]
				#Add it to the tabix queries dict as an identical sequence. Loop through the files of each dataset.
				for dataset in ToUpdate.keys():
					#If the dataset is not the genome dataset, create the tabix query
					if(dataset!="Genome"):
						#Loop through the files of the dataset  {dataset:[[filename.extension,directory,size],[...]]}
						for subfile in ToUpdate[dataset]:
							#Create the query for subfile A
							query="tabix ./temporary_directory/"+subfile[0]+"_A.gz "+oldID+":0 > ./temporary_directory/"+subfile[0]+"_"+oldID+":0_A"
							tabix_queries[query]="identical"
							#Add the information into the file crack. No need to add a discarded file since the sequence is identical
							file_crack["./temporary_directory/updated_"+subfile[0]+"_A"].append("./temporary_directory/"+subfile[0]+"_"+oldID+":0_A")
							#Only Alignment and annotation files need subfile B as well
							if(dataset=="Alignment" or dataset=="Annotation"):
								#Repeat with B
								query="tabix ./temporary_directory/"+subfile[0]+"_B.gz "+oldID+":0 > ./temporary_directory/"+subfile[0]+"_"+oldID+":0_B"
								tabix_queries[query]="identical"
								file_crack["./temporary_directory/updated_"+subfile[0]+"_B"].append("./temporary_directory/"+subfile[0]+"_"+oldID+":0_B")
				#It is identical but it their IDs might be different, so add it to the translate dict
				OldNewID_Dict[oldID]=newID
				#No need to keep looking for the equivalent sequence, break
				break
			#Check if the old reverse SHA1 is identical to the new forward SHA1
			elif(NewAssembly_Dict[newID][2]==OldAssembly_Dict[oldID][3]):
				#Add the info to the summary_Dict
				summary_Dict[oldID]=["reversed",str(len("".join(OldAssembly_Dict[oldID][1])))]
				#Add it to the tabix queries dict as an identical sequence. Loop through the files of each dataset.
				for dataset in ToUpdate:
					#If the dataset is not the genome dataset, create the tabix query
					if(dataset!="Genome"):
						#Loop through the files of the dataset
						for subfile in ToUpdate[dataset]:
							#Create the query for subfile A
							query="tabix ./temporary_directory/"+subfile[0]+"_A.gz "+oldID+":0 > ./temporary_directory/"+subfile[0]+"_"+oldID+":0_A"
							tabix_queries[query]=["reversed","_A",dataset,subfile[0],oldID,"./temporary_directory/"+subfile[0]+"_"+oldID+":0_A","./temporary_directory/updated_"+subfile[0]+"_"+oldID+":0_A","./temporary_directory/discarded_"+subfile[0]+":0_A",str(len("".join(OldAssembly_Dict[oldID][1])))]
							#Add the information into the file crack. No need to add a discarded file since the sequence is identical reversed.
							file_crack["./temporary_directory/updated_"+subfile[0]+"_A"].append("./temporary_directory/updated_"+subfile[0]+"_"+oldID+":0_A")
							#Only Alignment and annotation files need subfile B as well
							if(dataset=="Alignment" or dataset=="Annotation"):
								#Do the same with B
								query="tabix ./temporary_directory/"+subfile[0]+"_B.gz "+oldID+":0 > ./temporary_directory/"+subfile[0]+"_"+oldID+":0_B"
								tabix_queries[query]=["reversed","_B",dataset,subfile[0],oldID,"./temporary_directory/"+subfile[0]+"_"+oldID+":0_B","./temporary_directory/updated_"+subfile[0]+"_"+oldID+":0_B","./temporary_directory/discarded_"+subfile[0]+":0_B",str(len("".join(OldAssembly_Dict[oldID][1])))]
								#Add the information into the file crack. No need to add a discarded file since the sequence is identical reversed.
								file_crack["./temporary_directory/updated_"+subfile[0]+"_B"].append("./temporary_directory/updated_"+subfile[0]+"_"+oldID+":0_B")
				#Add Ids to the transalte dict
				OldNewID_Dict[oldID]=newID
				#No need to keep looking for the equivalent sequence, break
				break
	#Create a multifasta file for each assembly with those seqIDs that were not reversed or identical
	#First create for the old assembly
	with open(directory+"/Compare_OldAssembly.fa","w") as output_file:
		for oldID in OldAssembly_Dict.keys(): 
			if not oldID in OldNewID_Dict.keys():
				output_file.write(">"+oldID+"\n"+"".join(OldAssembly_Dict[oldID][1])+"\n")
	output_file.close()
	#Now with the new assemlbly
	with open(directory+"/Compare_NewAssembly.fa","w") as output_file:
		for newID in NewAssembly_Dict.keys(): 
			if not newID in OldNewID_Dict.values():
				output_file.write(">"+newID+"\n"+"".join(NewAssembly_Dict[newID][1])+"\n")
	output_file.close()
	#Run nucmer with the resulting multifastas. Use threads.
	print("\n\t\t - Running nucmer "+str(datetime.datetime.now()))
	ShellCommand=Popen("nucmer --forward --mum --threads="+str(threads)+" -p "+directory+"/ToFilter "+directory+"/Compare_OldAssembly.fa "+directory+"/Compare_NewAssembly.fa",shell=True).wait()
	#Filter the delta file
	print("\n\t\t - Filtering resulting alignment to obtain equivalent sequences acros versions "+str(datetime.datetime.now()))
	filter_delta(delta_file=directory+"/ToFilter.delta",out_file=alignment_pickle+"/Filtered.delta")
	#Run show-snps with the resulting filtered delta file
	print("\n\t\t - Identification of SNPs in sequences "+str(datetime.datetime.now()))
	ShellCommand=Popen("show-snps -H -T "+alignment_pickle+"/Filtered.delta > "+alignment_pickle+"/Filtered.snp",shell=True).wait()
	#Create a dictionary of snps in the filtered snp file {SeqID:"oldIndex oldChar newChar newIndex","...",}
	snp_dict=parse_snp_file(alignment_pickle+"/Filtered.snp")
	#Loop through the lines of the filtered delta file
	print("\n\t\t - Alignment finished. Now creating tabix queries "+str(datetime.datetime.now()))
	with open(alignment_pickle+"/Filtered.delta","r") as delta_file:
		#Ommit the first line
		delta_file.readline()
		#Initate variables
		prev_snp=0
		prev_end="1"
		#Loop through the lines
		for line in delta_file:
			#If the line starts with >, it is a new alignment
			if(line[0]==">"):
				line=line.split()
				current_oldID=line[0][1:]
				#Add the new and old Id to the OldNewID_Dict
				OldNewID_Dict[current_oldID]=line[1]
				#Re-start the previous end index and prev_snp
				prev_end="1"
				prev_snp=0
				#If the seqID has snps, load the new snp list
				if(current_oldID in snp_dict.keys()):
					snp_list=snp_dict[current_oldID]
				#Otherwise create an empty list
				else:
					snp_list=[]
			#Otherwise add the block (if it has whitespaces it is a block) and its snps into the delta_dict
			elif(" " in line):
				#Split the line
				line=line.rstrip().split()
				#Determine the current snps
				number_snp=int(line[4])
				current_snp=snp_list[int(prev_snp):int(number_snp+prev_snp)]
				prev_snp=number_snp+prev_snp
				#Add it to the tabix queries dict as a to be compared sequence. Loop through the files of each dataset.
				for dataset in ToUpdate:
					#If the dataset is not the genome dataset, create the tabix query
					if(dataset!="Genome"):
						#Loop through the files of the dataset
						for subfile in ToUpdate[dataset]:
							#Add the query A
							query="tabix ./temporary_directory/"+subfile[0]+"_A.gz "+current_oldID+":"+line[0]+"-"+line[1]+" > ./temporary_directory/"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_A"
							tabix_queries[query]=["compare","_A",dataset,subfile[0],current_oldID,"./temporary_directory/"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_A","./temporary_directory/updated_"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_A","./temporary_directory/discarded_"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_A",line[0:4],current_snp]
							#Add the to the file crack the updated and discarded files
							file_crack["./temporary_directory/updated_"+subfile[0]+"_A"].append("./temporary_directory/updated_"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_A")
							#Create a query B for alignment and annotation.
							if(dataset!="Variants"):
								#Do the same with B
								query="tabix ./temporary_directory/"+subfile[0]+"_B.gz "+current_oldID+":"+line[0]+"-"+line[1]+" > ./temporary_directory/"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_B"
								tabix_queries[query]=["compare","_B",dataset,subfile[0],current_oldID,"./temporary_directory/"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_B","./temporary_directory/updated_"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_B","./temporary_directory/discarded_"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_B",line[0:4],current_snp]
								#Add the to the file crack the updated and discarded files
								file_crack["./temporary_directory/updated_"+subfile[0]+"_B"].append("./temporary_directory/updated_"+subfile[0]+"_"+current_oldID+":"+line[0]+"_"+line[1]+"_B")
				#The previous end is now equal to the current region end
				prev_end=str(line[1])
	#Close the delta file
	delta_file.close()
	#If the oldID is not in the OldNew_Dict, it must be absent from the new assembly
	for old_key in OldAssembly_Dict.keys():
		if not(old_key in OldNewID_Dict.keys()):
			#Add the info to the summary_Dict
			summary_Dict[OldAssembly_Dict[old_key][0]]=["deleted",str(len("".join(OldAssembly_Dict[old_key][1])))]
	#Run show coords
	ShellCommand=Popen("show-coords -c -T -H "+alignment_pickle+"/Filtered.delta > "+alignment_pickle+"/summary.coords",shell=True).wait()
	#Return all the information in form of a list  [tabix_queries,OldNewID_Dict,alignment_pickle,summary_Dict,file_crack]
	return [tabix_queries,OldNewID_Dict,alignment_pickle,summary_Dict,file_crack]
