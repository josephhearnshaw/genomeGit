#!/usr/bin/env python

# PYTHON FUNCTIONS REQUIRED IN THE REPOSITORY AUTO-UPDATES PART I: ALIGNMENT OBTENTION.

#Make the imports
import hashlib
import pickle

#Create a store_variables function to store the provided variables into a pickle
def store_variables(variables,alignment_pickle):
	#Open a pickle
	with open(alignment_pickle+"/pickle", "wb") as pickle_file:
		#Store the variable list
		pickle.dump(variables, pickle_file,protocol=1)
	#Close the pickle
	pickle_file.close()

#Create a load_variables function to load variabled from a pickle file
def load_variables(alignment_pickle):
	#Open the pickle
	with open(alignment_pickle, "rb") as pickle_file:
		#Load the alignment pickle.
		variables=pickle.load(pickle_file)
	#Close the pickle
	pickle_file.close()
	#Return the data
	return variables

#Create a obtain_alignment_pickle to obtain the name of the pickle of two given assemblies
def obtain_alignment_pickle(old_assembly,new_assembly):
	#Inititate the content list
	content=[]
	#Open the new assembly and loop through the lines
	with open(new_assembly,"r") as new_assembly_file:
		for line in new_assembly_file:
			#Append the line to the content list
			line=line.rstrip()
			content.append(line)
	#Close the file
	new_assembly_file.close()
	#Obtain the hash
	new_assembly_hash=str(hashlib.sha1("".join(content)).hexdigest())
	#Reinititate the content list
	content=[]
	#Open the new assembly and loop through the lines
	with open(old_assembly,"r") as old_assembly_file:
		for line in old_assembly_file:
			#Append the line to the content list
			line=line.rstrip()
			content.append(line)
	#Close the file
	old_assembly_file.close()
	#Obtain the hash
	old_assembly_hash=str(hashlib.sha1("".join(content)).hexdigest())
	#Return the combination of both hashes
	return str("./Delta/"+new_assembly_hash+"_"+old_assembly_hash)



