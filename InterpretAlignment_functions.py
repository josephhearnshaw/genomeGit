#!/usr/bin/env python

# PYTHON FUNCTIONS REQUIRED IN THE REPOSITORY AUTO-UPDATES PART II: INTERPRETATION OF ALIGNEMNT

#Make the imports
import os
import sys
import time
import datetime
import shutil
from subprocess import Popen
from threading import Thread
from multiprocessing import Process
from parse_functions import parse_dataset

#Create global variables
tabix_queries={}
OldNewID_Dict={}
number_threads=1
processed=0
processing=0
updated={}
discarded={}

#Create a weld_files function to join back together files contained in filecrack dictionary {finalfilename:[subfile,]}
def weld_cracks(file_crack):
	#Loop through the keys of the dictionary
	for key in file_crack.keys():
		#Create a final file
		with open(key,"w") as final_file:
			#Loop through the values of the key
			for value in file_crack[key]:
				#Open the subfile in read mode and append all its lines to the final file
				with open(value,"r") as crack:
					for line in crack:
						final_file.write(line)
				#Close the file
				crack.close()
		#Close the file
		final_file.close()

#Create a function to eliminate repeated barcodes in sub SAM files (this can mess arround with the merge process)
def prepare_barcode(infile,outfile,dataset):
	#Initiate the index and open the input file
	old_index="0"
	with open(infile,"r") as input_file:
		#Open the output file and loop through the lines of the input file
		output_file=open(outfile,"w")
		#If annotation/alignment
		if(dataset!="Variants"):
			for line in input_file:
				line=line.split("\t")
				new_index=line[2]
				#If the new index is equal to the old one, skip the line
				if(new_index==old_index):
					continue
				#Otherwise include the line in the output file
				else:
					output_file.write("\t".join(line))
				old_index=new_index
		#Otherwise it is variants
		else:
			for line in input_file:
				line=line.split("\t")
				new_index=line[3]
				#If the new index is equal to the old one, skip the line
				if(new_index==old_index):
					continue
				#Otherwise include the line in the output file
				else:
					output_file.write("\t".join(line))
				old_index=new_index
		#Close all files
		output_file.close()
	input_file.close()

#Create a funtion to merge the two SAM/GFF/VCF subfiles
def merge_subfiles(dataset,subfile_name,template_length):
	#Inform the user
	print("\t\t - Merging metadata subfiles for: "+subfile_name+" "+str(datetime.datetime.now()))
	sys.stdout.flush()
	#Open the metadata
	metadata=open("./temporary_directory/"+subfile_name+"_metadata","r")
	#Read through all the comments of the metadata
	line_metadata=metadata.readline()
	while(line_metadata[0]=="@" or line_metadata[0]=="#"):
		line_metadata=metadata.readline()
	#If the dataset is alignemnt, there are two subfiles 
	if(dataset=="Alignment"):
		#Open the updated and discarded files
		updated_file=open("./temporary_directory/"+subfile_name,"a")	#Open the updated file in append mode since it already contains the comments
		discarded_file=open("./temporary_directory/"+subfile_name+".discarded","w")
		#First thing required is to prepare the barcodes in the subfiles: eliminate repeated barcodes
		prepare_barcode(infile="./temporary_directory/sorted_updated_"+subfile_name+"_A",outfile="./temporary_directory/BarcodeReady_"+subfile_name+"_A",dataset=dataset)
		prepare_barcode(infile="./temporary_directory/sorted_updated_"+subfile_name+"_B",outfile="./temporary_directory/BarcodeReady_"+subfile_name+"_B",dataset=dataset)
		#Open the files
		file_A=open("./temporary_directory/BarcodeReady_"+subfile_name+"_A","r")
		file_B=open("./temporary_directory/BarcodeReady_"+subfile_name+"_B","r")
		#Read the first line of metadata and the files A & B. Split them by the tabs
		line_metadata=line_metadata.split("\t")
		line_A=file_A.readline().rstrip().split("\t")
		line_B=file_B.readline().rstrip().split("\t")
		#Start looping through the metadata file
		while(len(line_metadata)>1):
			#If the one of the reads is already empty, discard the rest of the metadata
			if(len(line_A)!=3 or len(line_B)!=3):
				discarded_file.write("\t".join(line_metadata[1:]))
				line_metadata=metadata.readline().split("\t")
			#Otherwise updata the metadata
			else:
				#Obtain the barcodes
				barcode_metadata=line_metadata[0]
				barcode_A=line_A[2]
				barcode_B=line_B[2]
				#Check if the barcodes match. If they do, write the new line with the new TLENGHT
				if(barcode_metadata==barcode_A and barcode_metadata==barcode_B):
					#Read A new coordiantes and newID
					line_metadata[3]=line_A[0]
					line_metadata[4]=line_A[1]
					#New coordinate for read B
					line_metadata[8]=line_B[1]
					#New seqID for read B: if the seqID is identical in both reads, place a = and determine the TLENGTH
					if(line_B[0]==line_A[0]):
						line_metadata[7]="="
						line_metadata[9]=str(int(line_B[1])-int(line_A[1]))
						#Write the line in the output file,omit the barcode.
						updated_file.write("\t".join(line_metadata[1:]))
						#Discard if it is higher than the threshold or lower than zero
						#if(int(line_metadata[9])>template_length or int(line_metadata[9])<0): CURRENTLY ONLY NOT CHECKING TLENGTH VALUE
							#Tlength is too high/low, discard the entry
							#discarded_file.write("\t".join(line_metadata[1:]))
						#Otherwise it is good to go into the updated file
						#else:
							#Write the line in the output file,omit the barcode.
							#updated_file.write("\t".join(line_metadata[1:]))
					#Otherwise they map different reads, place a 0 as TLENGTH
					else:
						line_metadata[7]=line_B[0]
						line_metadata[9]="0"
						#Write the line in the output file,omit the barcode.
						updated_file.write("\t".join(line_metadata[1:]))
					#Read the next line of the three files and spplit them using the tabs
					line_metadata=metadata.readline().split("\t")
					line_A=file_A.readline().rstrip().split("\t")
					line_B=file_B.readline().rstrip().split("\t")
				#If the barcodes dont match, one of the reads was discarded at some point. If the barcode of both reads match, then both reads were discarded. Discard the 
				# entire entry and loop to the next line of the metadata
				elif(barcode_A==barcode_B):
					#Discard the entry
					discarded_file.write("\t".join(line_metadata[1:]))
					#Read the next line
					line_metadata=metadata.readline().split("\t")
				#Otherwise only one of the reads was discarded. Need to discard the entire line and re-adjust both reads
				else:
					#Discard the entry
					discarded_file.write("\t".join(line_metadata[1:]))
					#Determine which read was discarded
					if(int(barcode_A)>int(barcode_B)):
						#If line A was discarded, read the next line in the file B and split it by the tabs
						line_B=file_B.readline().rstrip().split("\t")
						#Read as well the next one in the metadata
						line_metadata=metadata.readline().split("\t")
					else:
						#If line B was discarded, read the next line in the file A and split it by the tabs
						line_A=file_A.readline().rstrip().split("\t")
						#Read as well the next one in the metadata
						line_metadata=metadata.readline().split("\t")
		#Close the files
		file_A.close()
		file_B.close()
		updated_file.close()
		discarded_file.close()
	#If it is annotation
	elif(dataset=="Annotation"):
		#Open the updated and discarded files
		updated_file=open("./temporary_directory/"+subfile_name,"a")	#Open the updated file in append mode since it already contains the comments
		discarded_file=open("./temporary_directory/"+subfile_name+".discarded","w")
		#First thing required is to prepare the barcodes in the subfiles: eliminate repeated barcodes
		prepare_barcode(infile="./temporary_directory/sorted_updated_"+subfile_name+"_A",outfile="./temporary_directory/BarcodeReady_"+subfile_name+"_A",dataset=dataset)
		prepare_barcode(infile="./temporary_directory/sorted_updated_"+subfile_name+"_B",outfile="./temporary_directory/BarcodeReady_"+subfile_name+"_B",dataset=dataset)
		#Open the files
		file_A=open("./temporary_directory/BarcodeReady_"+subfile_name+"_A","r")
		file_B=open("./temporary_directory/BarcodeReady_"+subfile_name+"_B","r")
		#Read the first line of metadata and the files A & B. Split them by the tabs
		line_metadata=line_metadata.split("\t")
		line_A=file_A.readline().rstrip().split("\t")
		line_B=file_B.readline().rstrip().split("\t")
		#Start looping through the metadata file
		while(len(line_metadata)>1):
			#If the one of the reads is already empty, discard the rest of the metadata
			if(len(line_A)!=3 or len(line_B)!=3):
				discarded_file.write("\t".join(line_metadata[1:]))
				line_metadata=metadata.readline().split("\t")
			#Otherwise updata the metadata
			else:
				#Obtain the barcodes
				barcode_metadata=line_metadata[0]
				barcode_A=line_A[2]
				barcode_B=line_B[2]
				#Check if the barcodes match. If they do, write the new line with the new coordinates and seqIDs
				if(barcode_metadata==barcode_A and barcode_metadata==barcode_B):
					#Only if both parts have been mapped to the same region, add them to the updated file
					if(line_A[0]==line_B[0]):
						#Part A new coordiantes and newID
						line_metadata[1]=line_A[0]
						line_metadata[4]=line_A[1]
						#New coordinate for part B
						line_metadata[5]=line_B[1]
						#Write the entry
						updated_file.write("\t".join(line_metadata[1:]))
					#Otherwise discard the read
					else:
						discarded_file.write("\t".join(line_metadata[1:]))
					#Read the next line of the three files and spplit them using the tabs
					line_metadata=metadata.readline().split("\t")
					line_A=file_A.readline().rstrip().split("\t")
					line_B=file_B.readline().rstrip().split("\t")
				#If the barcodes dont match, one of the reads was discarded at some point. If the barcode of both reads match, then both reads were discarded. Discard the 
				# entire entry and loop to the next line of the metadata
				elif(barcode_A==barcode_B):
					#Discard the entry
					discarded_file.write("\t".join(line_metadata[1:]))
					#Read the next line
					line_metadata=metadata.readline().split("\t")
				#Otherwise only one of the reads was discarded. Need to discard the entire line and re-adjust both reads
				else:
					#Discard the entry
					discarded_file.write("\t".join(line_metadata[1:]))
					#Determine which read was discarded
					if(int(barcode_A)>int(barcode_B)):
						#If line A was discarded, read the next line in the file B and split it by the tabs
						line_B=file_B.readline().rstrip().split("\t")
						#Read as well the next one in the metadata
						line_metadata=metadata.readline().split("\t")
					else:
						#If line B was discarded, read the next line in the file A and split it by the tabs
						line_A=file_A.readline().rstrip().split("\t")
						#Read as well the next one in the metadata
						line_metadata=metadata.readline().split("\t")
		#Close the files
		file_A.close()
		file_B.close()
		updated_file.close()
		discarded_file.close()
	#Otherwise it is variants: only one subfile
	else:
		#Open the updated and discarded files
		updated_file=open("./temporary_directory/"+subfile_name,"a")	#Open the updated file in append mode since it already contains the comments
		discarded_file=open("./temporary_directory/"+subfile_name+".discarded","w")
		#First thing required is to prepare the barcodes in the subfiles: eliminate repeated barcodes
		prepare_barcode(infile="./temporary_directory/sorted_updated_"+subfile_name+"_A",outfile="./temporary_directory/BarcodeReady_"+subfile_name+"_A",dataset=dataset)
		#Open the files
		file_A=open("./temporary_directory/BarcodeReady_"+subfile_name+"_A","r")
		#Read the first line of metadata and A. Split them by the tabs
		line_metadata=line_metadata.split("\t")
		line_A=file_A.readline().rstrip().split("\t")
		#Start looping through the metadata file
		while(len(line_metadata)>1):
			#If the one of the reads is already empty, discard the rest of the metadata
			if(len(line_A)!=4):
				discarded_file.write("\t".join(line_metadata[1:]))
				line_metadata=metadata.readline().split("\t")
			#Otherwise update the metadata
			else:
				#Obtain the barcodes
				barcode_metadata=line_metadata[0]
				barcode_A=line_A[3]
				#Check if the barcodes match. If they do, write the new line with the new TLENGHT
				if(barcode_metadata==barcode_A):
					#Read A new coordiantes and newID
					line_metadata[1]=line_A[0]
					line_metadata[2]=line_A[1]
					#Write the line in the output file,omit the barcode.
					updated_file.write("\t".join(line_metadata[1:]))
					#Read the next line
					line_metadata=metadata.readline().split("\t")
					line_A=file_A.readline().rstrip().split("\t")
				#Otherwise the variant was discarded. Need to discard the entire line and re-adjust
				else:
					#Discard the entry
					discarded_file.write("\t".join(line_metadata[1:]))
					#Re-adjust by reading the next line in the metadata and split it by the tabs
					line_metadata=metadata.readline().split("\t")
		#Close the files
		file_A.close()
		updated_file.close()
		discarded_file.close()

	#Close metadata
	metadata.close()

#Create a function to start a given function and its arguments in a subprocess in the shell. It will need to update the number of processes being processed and the number already processed
def call_batch(function,argument):
	#Load variables
	global processing
	global processed
	#Create a process and start it
	p = Process(target=function, args=(argument,))
	p.start()
	p.join()
	#Update variables
	processing-=1
	processed+=1

#Create a funciton to analyse the alignment of a given sequence and append the updated/discarded entries into the global output dictionaries
def update_sequence(query):
	#Load global variables
	global discarded
	global ToUpdate
	global OldNewID_Dict
	global tabix_queries
	#Execute the query
	ShellCommand=Popen(query,shell=True).wait()
	#If the query corresponds with an omitted region or an identical sequence, it is only necessary to perform the query. Exit the function.
	if(tabix_queries[query]=="identical"):
		#Exit the function
		return
	#Check if the sequence was reversed
	elif(tabix_queries[query][0]=="reversed"):
		#Open the query output file, updated and there is no need for discarded file (identical reverse sequence)
		query_outfile=open(tabix_queries[query][5],"r")
		updated_file=open(tabix_queries[query][6],"w")
		#Determine the lenght of the sequence
		length=int(tabix_queries[query][8])
		#Loop through the lines of the outfile
		for entry in subset_file:
			#Split the entry by the tabs and modify the entries so that they correspond with the reverse index
			entry=entry.split("\t")
			entry[1]=str(length-int(entry[1]))
			#Modify the oldID with the newID
			entry[0]=OldNewID_Dict[entry[0]]
			#Join the entry and append it into the overwrite dict
			updated_file.write("\t".join(entry))
		#Close the files
		query_outfile.close()
		updated_file.close()
		#Delete the subsets
		#os.remove("./temporary_directory/"+tabix_queries[query][3]+"_"+tabix_queries[query][4]+":0"+tabix_queries[query][1]) THIS IS FOR TESTING
	#Otherwise sequence alignment must be evaluated.
	else:
		#Determine the displacement factor
		displacement_factor=int(tabix_queries[query][8][2])-int(tabix_queries[query][8][0])
		#Open the query output file, updated and the discarded file
		query_outfile=open(tabix_queries[query][5],"r")
		updated_file=open(tabix_queries[query][6],"w")
		#Loop through the entries resulting out of the tabix query
		for entry in query_outfile:
			#Split the entry by the tabs
			entry=entry.split("\t")
			#Modify the oldID with the newID
			entry[0]=OldNewID_Dict[entry[0]]
			#Determine the entry index
			entry_index=int(entry[1])
			#Evaluate the entry inly if it is inside the block. NO NEED FOR THIS, IT WILL ALWAYS BE IN THE BLOCK THANKS TO TABIX
			#if(entry_index>=old_block_start) and (entry_index<=old_block_stop):
			#Determine if it is necessary to update the entry (maybe there is no change)
			if(displacement_factor==0) and (tabix_queries[query][8][1]==tabix_queries[query][8][3]) and (len(tabix_queries[query][9])==0):
				#Add the unchanged entry in the updated list.
				updated_file.write("\t".join(entry))
			#Otherwise there was a change
			else:
				#If there were no modifications, just add the displacement factor
				if(len(tabix_queries[query][9])==0):
					entry[1]=str(entry_index+displacement_factor)
					updated_file.write("\t".join(entry))
				#If there are some modifications in the alignment, both factors need to be considered
				else:
					entry_index=entry_index+displacement_factor
					#Update the entry index with the displacement related to the snps
					for snp in tabix_queries[query][9]:
						#Store the snp information
						snp=snp.split()
						snp[3]=int(snp[3])
						#If the snp is before the entry position, it affects the entry.
						# CHECK IT WITH THE NEW INDEX, NOT THE OLD!
						if(entry_index>=snp[3]):
							#If insertion, add one to the index
							if(snp[1]=="."):
								#Add one to the entry index
								entry_index+=1
							#It is a deletion
							elif(snp[2]=="."):
								#Remove one to the entry indexes
								entry_index-=1
							#Otherwise it is a substitution: If there is a substitution at the entry_index, modify the ref base (only if this is variants).
							elif(entry_index==snp[3] and tabix_queries[query][2]=="Variants"):
								entry[2]=str(snp[2])
					#Update the list
					entry[1]=str(entry_index)
					updated_file.write("\t".join(entry))
		#Close the files
		query_outfile.close()
		updated_file.close()
		
#Create a function to interpret the alignment of a tabix_queries
def interpret_alignment(queries,oldnew,threads,ToUpdate,tlength,filecrack):
	#Create global variables for the threads
	global tabix_queries
	tabix_queries=queries
	global OldNewID_Dict
	OldNewID_Dict=oldnew
	global number_threads
	number_threads=threads
	global processed
	processed=0
	global processing
	processing=0
	#Add the comments of the original files into the begining of the updated ones. Loop through the datasets and act accordingly.
	for dataset in ToUpdate.keys():	#{dataset:[[filename.extension,directory,size],[]...]}
		#No comments in the genome dataset
		if(dataset!="Genome"):
			#Loop through the files of the dataset
			for subfile in ToUpdate[dataset]:
				#Add the comments of the original file into the updated one
				ShellCommand=Popen("cp ./"+subfile[1]+"/Comments.txt ./temporary_directory/"+subfile[0],shell=True).wait()
	#From now on this part can be threaded, one thread analysing one tabix query at a time. The dictionaries must be made global otherwise they will be overwriten.
	query_count=0
	all_queries=tabix_queries.keys()
	prev_time=time.time()
	#While there are sequence alignments to be processed
	print("\n\t\t - Now processing tabix queries "+str(datetime.datetime.now()))
	while(processed!=len(all_queries)):
		#If the number of currently processing threads is lower than the thread number and there are commands to start execution, start their execution
		if (processing<number_threads) and (query_count<len(all_queries)):
			t=Thread(target=call_batch, args=(update_sequence,all_queries[query_count],))
			t.daemon = True
			t.start()
			query_count+=1
			processing+=1
		#Inform the user of the progress. Do this only every 30 seconds.
		if(time.time()-prev_time>30):
			print("\t\t\t "+str(processed)+" queries processed out of "+str(len(all_queries))+" "+str(datetime.datetime.now()))
			prev_time=time.time()
	#When the threads are done, merge the files in the filecrack. Inform the user.
	print("\n\t\t - Now concatenating updated subfiles resulting from multi-threaded mode "+subfile[0]+" "+str(datetime.datetime.now()))
	sys.stdout.flush()
	weld_cracks(filecrack)
	#Finally, when the files are created it is required to sort them. Loop through the datasets.
	for dataset in ToUpdate.keys():
		if(dataset!="Genome"):
			#Loop through the files of the dataset
			for subfile in ToUpdate[dataset]:
				#Inform the user
				print("\n\t\t - Now parsing updated file "+subfile[0]+ " into the repository structure "+str(datetime.datetime.now()))
				sys.stdout.flush()
				#If the dataset is alignment
				if(dataset=="Alignment"):
					#Append the missing unmaped reads intot he updated files (these reads wont be included otherwise)
					ShellCommand=Popen("tabix ./temporary_directory/"+subfile[0]+"_A.gz *:0 >> ./temporary_directory/updated_"+subfile[0]+"_A",shell=True).wait()
					ShellCommand=Popen("tabix ./temporary_directory/"+subfile[0]+"_B.gz *:0 >> ./temporary_directory/updated_"+subfile[0]+"_B",shell=True).wait()
					#Sort the reads in both A and B files using the barcode
					ShellCommand=Popen("sort --numeric-sort -k 3 ./temporary_directory/updated_"+subfile[0]+"_A > ./temporary_directory/sorted_updated_"+subfile[0]+"_A",shell=True).wait()
					ShellCommand=Popen("sort --numeric-sort -k 3 ./temporary_directory/updated_"+subfile[0]+"_B > ./temporary_directory/sorted_updated_"+subfile[0]+"_B",shell=True).wait()
				elif(dataset=="Annotation"):
					#Sort the reads in both A and B files using the barcode
					ShellCommand=Popen("sort --numeric-sort -k 3 ./temporary_directory/updated_"+subfile[0]+"_A > ./temporary_directory/sorted_updated_"+subfile[0]+"_A",shell=True).wait()
					ShellCommand=Popen("sort --numeric-sort -k 3 ./temporary_directory/updated_"+subfile[0]+"_B > ./temporary_directory/sorted_updated_"+subfile[0]+"_B",shell=True).wait()
				#Otherwise it is variants
				else:
					#Sort the reads in both A and B files using the barcode
					ShellCommand=Popen("sort --numeric-sort -k 4 ./temporary_directory/updated_"+subfile[0]+"_A > ./temporary_directory/sorted_updated_"+subfile[0]+"_A",shell=True).wait()
				#Merge all the files
				merge_subfiles(dataset=dataset,subfile_name=subfile[0],template_length=int(tlength))
				#Parse the updated file into the repository (not necessary to git add, that is done in the genomegit main wrapper). First, delete the old directory.
				shutil.rmtree("./"+dataset+"/"+subfile[0])
				#Create a new directory a move inside
				os.mkdir("./"+dataset+"/"+subfile[0])
				os.chdir("./"+dataset+"/"+subfile[0])
				#Parse the file
				parse_dataset(dataset=dataset,input_path="../../temporary_directory/"+subfile[0],size=str(os.path.getsize("../../temporary_directory/"+subfile[0]))[:-1],update="1")
				#Go back to the base directory
				os.chdir("../../")
				#Add the discarded entries into the file directory
				shutil.copyfile("./temporary_directory/"+subfile[0]+".discarded","./"+dataset+"/"+subfile[0]+"/Discarded")
