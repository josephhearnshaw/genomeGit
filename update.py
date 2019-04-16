#!/usr/bin/env python

# INSTRUCTIONS: The script will update the coordinates of the variants and annotation files of a given
# genome assembly to match those of a new one. Arguments:
# update_dependent_datasets.py -new_file- -threads- -file_size- -template_length-

# Make the imports
import os
import sys
import datetime
from subprocess import Popen
from update_functions import detect_updates
from reconstruct_functions import reconstruct_dataset
from ObtainAlignment_functions import obtain_alignment
from StoreAlignment import store_variables, load_variables, obtain_alignment_pickle
from InterpretAlignment_functions import interpret_alignment
import shutil

# Load arguments: new_assembly, threads size, tlength, which aligner to use, and associated mashmap args
new_assembly = str(sys.argv[1])
number_threads = int(sys.argv[2])
file_size = str(sys.argv[3])
template_length = int(sys.argv[4])
aligner_switch = int(sys.argv[5])
aligner = ""
percent_identity = int(sys.argv[6])
kmer = int(sys.argv[7])
segLength = int(sys.argv[8])


def aligner_status(ToUpdate):
    """Informs the user of the aligner selected
    Requires ToUpdate and global variables defined at the start of Update.py """

    if (ToUpdate == "empty"):
        print("No datasets to be updated were detected.\n")
        sys.exit()
    else:
        if aligner_switch == 1:
            aligner = "MashMap + Nucmer"
            print("\n\n**Aligner chosen is MashMap + Nucmer (GenomeGit 2.1) with the following "
                  "parameters for MashMap:**\n\n\t -percent identity: {}"
                  "\n\n\t -Kmer: {}\n\n\t -Segment length: {}\n\n"
                  .format(percent_identity, kmer, segLength))
        else:
            aligner = "Nucmer"
            print("\n\n**Aligner chosen is Nucmer only (GenomeGit 2.0)**\n\n")
    return aligner


def update_inform_user(ToUpdate):
    """ Simply checks if there are any files to be updated and reports the file information, including
    file size and filename.

    Requires: ToUpdate"""

    aligner = aligner_status(ToUpdate)

    print("Number of threads selected: {}\nMaximum insert size: {}\nThe following files that "
          "need to be updated have been detected:\n".format(str(number_threads), str(template_length)))
    for dataset in ToUpdate.keys():
        # If the dataset is not a genome dataset
        if(dataset != "Genome"):
            print("\t-Files detected in the " + dataset + " dataset\n")
            # If there are files (len of dataset > 0), print the files out
            if(len(ToUpdate[dataset]) != 0):
                for subfile in ToUpdate[dataset]:
                    # If the size of the file is lower than 1MB then print the below message
                    if(subfile[2] == "1"):
                        print("\t\t--{} (<1 MB)\n".format(subfile[0]))
                    # Alternatively, print the below if >1MB
                    else:
                        print(
                            "\t\t--{} ({} MB)\n".format(subfile[0], subfile[2]))
            # If there are no files in this dataset, inform the user
            else:
                print("\t\t--No files detected in this {} dataset.\n".format(dataset))

    # Inform the user that the update will start now
    # If the file size of the genome data is 1 Mb perform the following:
    if(ToUpdate["Genome"][0][2] == "1"):
        if(file_size == "1"):
            print("\n*** NOW STARTING UPDATE OF THE REPOSITORY: {} (<1 MB) ---> {} (<1 MB)***\n"
                  .format(ToUpdate["Genome"][0][0], os.path.basename(new_assembly)))
        else:
            print("\n*** NOW STARTING UPDATE OF THE REPOSITORY: {} (<1 MB) ---> {} ({} MB)***\n"
                  .format(ToUpdate["Genome"][0][0], os.path.basename(new_assembly), file_size))
    else:
        if(file_size == "1"):
            print("\n*** NOW STARTING UPDATE OF THE REPOSITORY: {} ({} MB) --> {} (<1 MB)***\n"
                  .format(ToUpdate["Genome"][0][0], ToUpdate["Genome"][0][2],
                          os.path.basename(new_assembly)))
        else:
            print("\n*** NOW STARTING UPDATE OF THE REPOSITORY: {} ({} MB) ---> {} ({} MB)***\n"
                  .format(ToUpdate["Genome"][0][0], ToUpdate["Genome"][0][2],
                          os.path.basename(new_assembly), file_size))


def reconstruct_annotation_variants(ToUpdate):
    """ Reconstructs the variant and annotation dataset in the temporary directory.

    Requires: ToUpdate ({dataset:[[filename.extension,directory,size],[],...]}) which contains INFORMATION
    of what requires updating. """

    for dataset in ToUpdate.keys():
        # If the data are variants and annotation datasets, and there is data present, perform the following block
        if(dataset != "Genome" and len(ToUpdate[dataset]) != 0):
            # Loop through the files of the dataset and reconstruct the file in the temporary
            # directory with the same original name
            for subfile in ToUpdate[dataset]:
                reconstruct_dataset(
                    size=1, directory=subfile[1], output_file="./temporary_directory/{}".format(subfile[0]),
                    mode=dataset, update=True, seqID="0", region="0")
                # Create an appropiate tabix library for the file. If it is annotation or
                # alignment there are two pseudo files
                if(dataset == "Annotation" or dataset == "Alignment"):
                    # Get all the headers (starts wiht #) for the alignment files, and all non header data,
                    # sort column 1 and column 2 (numerically).
                    # Then compress and create indexes and use tabix to create the library
                    Popen("(grep ""^#"" ./temporary_directory/{}_A; grep -v ""^#"" ./temporary_directory/{}_A | "
                          "sort -V -k1,1 -k2,2n) | bgzip > ./temporary_directory/{}_A.gz; "
                          "tabix -b 2 -e 2 ./temporary_directory/{}_A.gz;"
                          .format(subfile[0], subfile[0], subfile[0], subfile[0]), shell=True).wait()
                    # Do this for dataset B too
                    Popen("(grep ""^#"" ./temporary_directory/{}_B; grep -v ""^#"" ./temporary_directory/{}_B | "
                          "sort -V -k1,1 -k2,2n) | bgzip > ./temporary_directory/{}_B.gz; "
                          "tabix -b 2 -e 2 ./temporary_directory/{}_B.gz;"
                          .format(subfile[0], subfile[0], subfile[0], subfile[0]), shell=True).wait()
                elif(dataset == "Variants"):
                    # Same as above but with variant data instead
                    Popen("(grep ""^#"" ./temporary_directory/{}_A; grep -v ""^#"" ./temporary_directory/{}_A | "
                          "sort -V -k1,1 -k2,2n) | bgzip > ./temporary_directory/{}_A.gz; "
                          "tabix -b 2 -e 2 ./temporary_directory/{}_A.gz;"
                          .format(subfile[0], subfile[0], subfile[0], subfile[0]), shell=True).wait()
                # Otherwise there was an error
                else:
                    print(
                        "***INTERNAL ERROR*** DATASET NOT RECOGNIZED: {}".format(dataset))


def obtain_variables(alignment_pickle):
    """ This function will obtain the alignment if present, otherwise it'll call to perform a new alignment if not present.
    It will return the variables obtained, including the new alignment if performed.

    Requires: alignment_pickle"""

    # If there is an alignment already stored in the repository
    if(os.path.isdir(alignment_pickle)):
        # Inform the user
        print("\n\t*PART II. OBTAINING GENOME ALIGNMENT: STORED ALIGNMENT DETECTED, NOW LOADING SAVED DATA.*")
        print("\t" + str(datetime.datetime.now()))
        # Load the variables stored in the pickle [tabix_queries,alignment_pickle,summary_Dict,file_crack]
        variables = load_variables(alignment_pickle + "/pickle")
    # Otherwise it is necessary to obtain the alignment
    else:
        # Inform the user
        print("\n\t*PART II. OBTAINING GENOME ALIGNMENT: NO STORED ALIGNMENT DETECTED, NOW CREATING NEW ALIGNMENT*")
        print("\t" + str(datetime.datetime.now()))
        # Obtain the variables
        variables = obtain_alignment(old_assembly="./temporary_directory/genome_old.fa",
                                     new_assembly=new_assembly, directory="./temporary_directory",
                                     threads=number_threads, ToUpdate=ToUpdate, alignment_pickle=alignment_pickle,
                                     aligner_switch=aligner_switch, percent_identity=percent_identity,
                                     kmer=kmer, segLength=segLength)
    return variables

####################################################
# PART 0.DETERMINE FILES THAT NEED TO BE UPDATED #
####################################################


# Determine the files to be updated. ToUpdate = {dataset:[[filename.extension,directory,size],[],...]}
ToUpdate = detect_updates("./RepoMap.txt")
# If the there are no files to update, inform the user
update_inform_user(ToUpdate)


###############################################
# PART 1. RECONSTRUCTION OF REPOSITORY DATA #
###############################################

# Inform the user
print("\n\t*PART I. RECONSTRUCTION OF REPOSITORY DATA.*")
print("\t {}".format(str(datetime.datetime.now())))
# Create temporary directory
os.mkdir("./temporary_directory")
# Reconstruct the necessary datasets. First the old genome.
reconstruct_dataset(size=60, directory="./Genome", output_file="./temporary_directory/genome_old.fa",
                    mode="Genome", seqID="0", region="0")
# Reconstruct all the files related to the variants and annotation datasets
reconstruct_annotation_variants(ToUpdate)

##############################################
# PART 2. OBTAIN THE ALIGNMENT INFORMATION #
##############################################

# Obtain the information of the alignment between both assemblies.
# tabix queries:	tabix ./temporary_directory/filename.ext_AB.gz chr:x-y > ./temporary_directory/filename.ext_chr:x_y_AB
#
# {tabix_query:["compare","_A/_B//","Annotation/Variants/Alignment",filename.ext,oldID,query_outfile,sub_updated,sub_discarded,[old_block_start,old_block_stop,new_block_start,new_block_stop,block_modifications],[snps...]]}
# {tabix_query:["reversed","_A/_B//","Annotation/Variants/Alignment",filename.ext,oldID,query_outfile,sub_updated,sub_discarded,length]}
# {tabix_query:"omitted"}
# {tabix_query:"identical"}
#
# file_crack {finalfilename:[sub_updated/discarded...]}
#
# finalfilename:	updated_filename.ext	/	discarded_filename.ext
# sub_updated/discarded:	updated_filename.ext_chr:x_y_AB	/	discarded_filename.ext_chr:x_y_AB
#
# summary_Dict {oldID:[status,length]}

# Determine the alignment pickle
alignment_pickle = obtain_alignment_pickle("./temporary_directory/genome_old.fa",
                                           new_assembly)
variables = obtain_variables(alignment_pickle)
store_variables(variables=variables, alignment_pickle=alignment_pickle)
queries, alignment_pickle, summary_Dict, file_crack, oldSeqs, newSeqs = variables

###################################################################
# PART 3. START THE INTERPRETATION OF THE ALIGNMENT INFORMATION #
###################################################################

# Inform the user
print("\n\t*PART III. INTERPRETATION OF THE ALIGNMENT INFORMATION AND CREATION OF UPDATED FILES")
print("\t{}".format(str(datetime.datetime.now())))
# Interpret the information contained in the delta_dict and obtain the updated files.
interpret_alignment(queries=queries, threads=number_threads, ToUpdate=ToUpdate,
                    tlength=template_length, filecrack=file_crack)

# Inform the user the update is completed
print("\n\t***UPDATE COMPLETED: NOW PARSING THE GENOME DATASET***")
