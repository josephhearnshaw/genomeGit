#!/usr/bin/env python

# PYTHON FUNCTIONS USED FOR FILE RECOSTRUCTION

# Make the imports
import re
import os
import sys
import subprocess
from subprocess import Popen


def MessageToHash(message):
    """
    Function to return the commit hash for a given commit message
    """

    # Isolate the first commit hash of a given message and return it
    commit_hash = subprocess.check_output(
        'git log --grep="' + message + '" --format=%H -n 1')
    return str(commit_hash).rstrip()


def get_files(Dir, Seq, RegionStart, RegionStop):
    """
    Function to create a list with the subfiles included in a region of interest
    """

    # Create an empty list
    ToReconstruct = []
    # If the user is interested in a region in particular
    if(RegionStart != "none"):
        # Determine the 1MB regions of interest
        RegionStart = RegionStart / 1000000
        RegionStop = RegionStop / 1000000
        # If the regions are identical
        if(RegionStart == RegionStop):
            all_regions = [RegionStart]
        # Otherwise append them
        else:
            all_regions = range(RegionStart, RegionStop + 1)
        # Append the regions into a list, only if the file exists
        for region in all_regions:
            if(os.path.isfile(Dir + "/Region" + str(region))):
                ToReconstruct.append(Dir + "/Region" + str(region))
        # Return the list
        return ToReconstruct
    # Otherwise the user wants the entire sequence
    else:
        # Open the seqID map and
        with open(Dir + "/SeqIDs_Map.txt", "r") as IDfile:
            # Read the first line
            line = IDfile.readline()
            # Loop through the lines
            while(line):
                # If the line is the sequence
                if(line[0] == ">" and line.split("\t")[0][1:] == Seq):
                    # Read the next line
                    line = IDfile.readline()
                    # Append all the lines until the next sequence into the list
                    while(line and line[0] != ">"):
                        ToReconstruct.append(line.rstrip())
                        # Read the next line
                        line = IDfile.readline()
                    # Close the file
                    IDfile.close()
                    # Return the list
                    return ToReconstruct
                # Read the next line
                line = IDfile.readline()


def check_ID(Seq, Dir):
    """
    Function to determine if a sequence exists in a file
    """

    # Open the SeqIDs_Map.txt file
    with open(Dir + "/SeqIDs_Map.txt", "r") as IDfile:
        # Loop through the lines
        for line in IDfile:
            # Check if the line corresponds with a seqID
            if(line[0] == ">"):
                # Split the line
                line = line.split("\t")
                # Check if the seqID is the one of interest
                if(line[0][1:] == Seq):
                    return line[len(line) - 1].rstrip()
    # Close the file
    IDfile.close()
    # SeqID not found
    return "none"


def reconstruct_dataset(size, directory, output_file, mode, seqID="0", region="0",
                        update=False, filemode="w", convertToBam="0"):
    """
    Function to reconstruct files back from their folders
    """

    # GENOME RECONSTRUCTION

    # If the user wants to reconstruct a genome dataset
    if(mode == "Genome"):
        # Open a reconstruct file which will be the reconstructed fasta
        reconstruct = open(output_file, filemode)

        # PART 1. STORE THE SEQUENCE IDs

        # Open the seqID file and store each of the IDs (one per line) in a list
        seq_list = []
        seq_file = open(directory + "/SeqIDs.txt", "r")
        for line in seq_file:
            seq_list.append(line)
        # Close the seq ID file
        seq_file.close()

        # PART 2. APPEND THE CONTIG SEQUENCES TO THE RECONSTRUCTED FASTA FILE

        # Open the input map file and read the first line
        map = open(directory + "/Map.txt", "r")
        # Initiate required variables
        seq_count = -1
        current_sequence = ""
        # Loop through the lines of the map file
        for line in map:
            # If the line starts with ">", it corresponds with a new sequence
            if(line[0] == ">"):
                seq_count += 1
                # Unless this is the first sequence, write the currently stored sequence in the reconstructed file.
                if(seq_count != 0):
                    # Write the seq ID
                    reconstruct.write(seq_list[seq_count - 1])
                    # Write the sequence of nucleotides using the line size specified
                    reconstruct.write('\n'.join(
                        current_sequence[y:y + size] for y in range(0, len(current_sequence), size)) + "\n")
                    # Empty current sequence
                    current_sequence = ""
            # If the line does not start with an > or a N, it corresponds with the path to a sequence file
            elif(line[0] != "N"):
                # Split the line using the \t. The first field should correspond with the full path to the contig file.
                line = re.split("\t", line)
                input_file = open(directory + line[0][1:], "r")
                # Open the file, store all its lines and remove all whitespaces and \n characters
                file_content = input_file.read()
                file_content = re.sub('\s+', '', file_content)
                # Store all the sequence of the contig file
                current_sequence = current_sequence + file_content
            # Otherwise the line corresponds to a gap
            else:
                # Split the line using the \t and then with the \s. The start is
                # the first number and the stop the second one
                line = re.split("\t", line)
                line = line[1]
                line = re.split(" ", line)
                # Add the gap into the current sequence. The size of the gap is the difference
                # between the stop and the start plus one
                current_sequence = current_sequence + "N" * \
                    int(int(line[1]) - int(line[0]) + 1)
        # Close the file
        map.close()
        # Add the last sequence left out of the loop
        reconstruct.write(seq_list[seq_count])
        reconstruct.write('\n'.join(
            current_sequence[y:y + size] for y in range(0, len(current_sequence), size)) + "\n")
        # Close the reconstruct ile
        reconstruct.close()

    # UPDATE RECONSTRUCTION

    # If this forms part of an update
    elif(update):
        # If the dataset to reconstruct is aligment
        if(mode == "Alignment"):
            # Open a second reconstruct file
            reconstruct_A = open(output_file + "_A", "w")
            reconstruct_B = open(output_file + "_B", "w")
            reconstruct_metadata = open(output_file + "_metadata", "w")
            # Open the comment file and append it to the reconstructed metadata file
            comments = open(directory + "/Comments.txt", "r")
            for line in comments:
                reconstruct_metadata.write(line)
            # Open the seqID map
            seqID_map = open(directory + "/SeqIDs_Map.txt", "r")
            # Loop through the map
            line_count = 0
            for line in seqID_map:
                # If the line starts with >, it is a new seqID. Ignore it
                if(line[0] != ">"):
                    # Open the file specified in the line
                    comments.close()
                    comments = open(directory + line.rstrip(), "r")
                    # Append all the lines of the file to the reconstructed files
                    for l in comments:
                        line_count += 1
                        # Add the metadata to the metadata file
                        reconstruct_metadata.write(str(line_count) + "\t" + l)
                        # Split the line
                        l = l.split("\t")
                        # Add the first pair to the reconstruct A
                        reconstruct_A.write(
                            l[2] + "\t" + l[3] + "\t" + str(line_count) + "\n")
                        # If the left read has an = symbol, it maps the same seqID
                        if(l[6] == "="):
                            reconstruct_B.write(
                                l[2] + "\t" + l[7] + "\t" + str(line_count) + "\n")
                        # Otherwise it is a different seqID
                        else:
                            reconstruct_B.write(
                                l[6] + "\t" + l[7] + "\t" + str(line_count) + "\n")
            # Close the files
            reconstruct_metadata.close()
            reconstruct_A.close()
            reconstruct_B.close()
            comments.close()
            seqID_map.close()
        # If the dataset is variants
        elif(mode == "Variants"):
            # Open a second reconstruct file
            reconstruct = open(output_file + "_A", "w")
            reconstruct_metadata = open(output_file + "_metadata", "w")
            # Open the comment file and append it to the reconstructed metadata file
            comments = open(directory + "/Comments.txt", "r")
            for line in comments:
                reconstruct_metadata.write(line)
            # Open the seqID map
            seqID_map = open(directory + "/SeqIDs_Map.txt", "r")
            # Loop through the map
            line_count = 0
            for line in seqID_map:
                # If the line starts with >, it is a new seqID. Ignore it
                if(line[0] != ">"):
                    # Open the file specified in the line
                    comments.close()
                    comments = open(directory + line.rstrip(), "r")
                    # Append all the lines of the file to the reconstructed files
                    for l in comments:
                        line_count += 1
                        # Add the metadata to the metadata file
                        reconstruct_metadata.write(str(line_count) + "\t" + l)
                        # Split the line
                        l = l.split("\t")
                        # Add the seqID, the index and the reference nucleotide into the reconstruct
                        reconstruct.write(
                            # l[0] + "\t" + l[1] + "\t" + l[3] + "\t" + str(line_count) + "\n")
                            l[0] + "\t" + l[1] + "\t" + l[3] + "\t" + l[4] + "\t" + str(line_count) + "\n")
            # Close the files
            reconstruct_metadata.close()
            reconstruct.close()
            comments.close()
            seqID_map.close()
        # If the dataset to reconstruct is annotation, create two subfiles (one for each coordinate)
        elif(mode == "Annotation"):
            # Open a second reconstruct file
            reconstruct_A = open(output_file + "_A", "w")
            reconstruct_B = open(output_file + "_B", "w")
            reconstruct_metadata = open(output_file + "_metadata", "w")
            # Open the comment file and append it to the reconstructed metadata file
            comments = open(directory + "/Comments.txt", "r")
            for line in comments:
                reconstruct_metadata.write(line)
            # Open the seqID map
            seqID_map = open(directory + "/SeqIDs_Map.txt", "r")
            # Loop through the map
            line_count = 0
            for line in seqID_map:
                # If the line starts with >, it is a new seqID. Ignore it
                if(line[0] != ">"):
                    # Open the file specified in the line
                    comments.close()
                    comments = open(directory + line.rstrip(), "r")
                    # Append all the lines of the file to the reconstructed files
                    for l in comments:
                        line_count += 1
                        # Add the metadata to the metadata file
                        reconstruct_metadata.write(str(line_count) + "\t" + l)
                        # Split the line
                        l = l.split("\t")
                        # Add the first pair to the reconstruct A
                        reconstruct_A.write(
                            l[0] + "\t" + l[3] + "\t" + str(line_count) + "\n")
                        reconstruct_B.write(
                            l[0] + "\t" + l[4] + "\t" + str(line_count) + "\n")
            # Close the files
            reconstruct_metadata.close()
            reconstruct_A.close()
            reconstruct_B.close()
            comments.close()
            seqID_map.close()
        # There was an error somwhere
        else:
            print("*INTERNAL ERROR. RECONSTRUCTION MODE NOT RECOGNISED: " + mode)

    # REGION OF INTEREST RECONSTRUCTION

    # If the user is interested in sequence of interest
    elif(seqID != "0"):
        # Get the direcetory of the sequence of interest
        SeqDir = check_ID(Seq=seqID, Dir=directory)
        # Check that the seqID exists in the file, if this is not a dataset extraction exit the program if it is not present
        if(filemode != "a" and SeqDir == "none"):
            # Inform the user and exit if the ID is not in the file
            print(
                "*ERROR: The provided sequence to be extracted is not present in the file.")
            sys.exit()
        # If this is a dataset extraction and the sequence is not in the current file, return the function to jump to the next file
        elif(SeqDir == "none"):
            return
        # Check if the user is interested in a region in particular
        elif(region != "0"):
            # Check if the regions are integers
            try:
                # Determine the region start / stop
                region_start = int(region.split("-")[0])
                region_stop = int(region.split("-")[1])
                # Open a recontruct file
                reconstruct = open(output_file, filemode)
                # Only append comments if filemode is different than append (no comments if the user reconstruct a region)
                if(filemode != "a"):
                    comment_file = open(directory + "/Comments.txt", "r")
                    for line in comment_file:
                        reconstruct.write(line)
                    comment_file.close()
                # Determine the files of interest (those files containing the region of interest) and create a list [filename,filename...]
                ToReconstruct = get_files(
                    Dir=directory + SeqDir, Seq=seqID, RegionStart=region_start, RegionStop=region_stop)
                # If variants
                if(mode == "Variants"):
                    # Loop through the list of files and append thos entries inside the region
                    for filename in ToReconstruct:
                        # Open the file
                        with open(filename, "r") as subfile:
                            # Loop through the lines
                            for line in subfile:
                                # Split the line into its fields
                                line = line.split("\t")
                                # Check that the line index is within the region, otherwise stop
                                if(int(line[1]) < region_stop and int(line[1]) > region_start):
                                    reconstruct.write("\t".join(line))
                        # Close the file
                        subfile.close()
                # If alignment
                elif(mode == "Alignment"):
                    # Loop through the list of files and append thos entries inside the region
                    for filename in ToReconstruct:
                        # Open the file
                        with open(filename, "r") as subfile:
                            # Loop through the lines
                            for line in subfile:
                                # Split the line into its fields
                                line = line.split("\t")
                                # Check that the line index is within the region, otherwise stop
                                if(int(line[3]) < region_stop and int(line[3]) > region_start):
                                    reconstruct.write("\t".join(line))

                        # Check if user wants to convert to BAM
                        if (convertToBam == "1"):
                            convertSAMtoBAM(subfile.name)

                # Otherwise it is annotation
                else:
                    # Loop through the list of files and append thos entries inside the region
                    for filename in ToReconstruct:
                        # Open the file
                        with open(filename, "r") as subfile:
                            # Loop through the lines
                            for line in subfile:
                                # Split the line into its fields
                                line = line.split("\t")
                                # Check that the line index is within the region, otherwise stop
                                if(int(line[3]) < region_stop and int(line[3]) > region_start or int(line[4]) < region_stop and int(line[4]) > region_start):
                                    reconstruct.write("\t".join(line))
                        # Close the file
                        subfile.close()
            # If the user does not provide an integer
            except ValueError:
                # Inform the user and exit
                print(
                    '*ERROR: You must provide two integers separated with a "-" as a region to be extracted.')
        # Otherwise users is only interested in the entire sequence
        else:
            # Open a recontruct file
            reconstruct = open(output_file, filemode)
            # Only append comments if filemode is different than append (no comments if the user reconstruct a region)
            if(filemode != "a"):
                comment_file = open(directory + "/Comments.txt", "r")
                for line in comment_file:
                    reconstruct.write(line)
                comment_file.close()
            # Determine the files of interest (those files containing the region of interest) and create a list [filename,filename...]
            ToReconstruct = get_files(
                Dir=directory, Seq=seqID, RegionStart="none", RegionStop="none")
            # Loop through the files
            for filename in ToReconstruct:
                # Open the file
                with open(directory + filename, "r") as subfile:
                    # Loop through the lines
                    for line in subfile:
                        # Write the line into the reconstruct file
                        reconstruct.write(line)
                # Close the file
                subfile.close()
            # Close the reconstruct file
            reconstruct.close()

            # Check if user wants to convert to bam
            if (convertToBam == "1"):
                convertSAMtoBAM(reconstruct.name)

    # REGULAR RECONSTRUCTION

    # Otherwise this is a normal reconstruction of annotation/variants/alignment dataset
    else:
        # Open a reconstruct file which will be the reconstructed file
        reconstruct = open(output_file, filemode)
        # Open the comment file and append it to the reconstructed file
        comments = open(directory + "/Comments.txt", "r")
        # Only append comments if filemode is different than append (no comments if the user reconstruct a region)

        # Add headers if it's a bam or sam file
        if(output_file.endswith('.bam') or output_file.endswith('.sam') or filemode != 'a'):
            for line in comments:
                reconstruct.write(line)
        # Open the seqID map
        seqID_map = open(directory + "/SeqIDs_Map.txt", "r")
        # Loop through the map
        for line in seqID_map:
            # If the line starts with >, it is a new seqID. Ignore it
            if(line[0] != ">"):
                # Open the file specified in the line
                comments.close()
                comments = open(directory + line.rstrip(), "r")
                # Append all the lines of the file to the reconstruct
                for l in comments:
                    reconstruct.write(l)

        # Close the files
        reconstruct.close()
        seqID_map.close()
        comments.close()

        # Check if the user wants to convert to BAM
        if (convertToBam == "1"):
            convertSAMtoBAM(reconstruct.name)


def convertSAMtoBAM(filename):
    """
    Convert from SAM to BAM
    """

    print("Now converting SAM to BAM")
    output_name = filename.replace('.sam', '.bam')
    Popen("samtools view -S -b -h -o " + output_name +
          " " + filename, shell=True).wait()
    print("Conversion done")
    os.remove(filename)


def obtain_file(target, mode):
    """
    Obtain the dataset and line size of a given file
    """

    # If the provided target is a filename
    if(mode == "filename"):
        # Open the repo map in read mode
        repomap = open("./RepoMap.txt", "r")
        # Loop through the lines of the map
        for line in repomap:
            line = line.split("\t")
            # If the file corresponds with the target return its features
            if(line[0] == target):
                # Close the file and exit the function
                repomap.close()
                return line[1], line[2], line[3], line[4]
        # If there is no match for the file of interest, is not stored in the repo
        repomap.close()
        return "error", "error", "error", "error"
    # Otherwise the target provided is a dataset
    else:
        # Create an empty file list to be returned as a result [[filename,dataset,directory,linesize],[...]]
        file_list = []
        # Open the repo map in read mode
        with open("./RepoMap.txt", "r") as repomap:
            # Loop through the lines of the map
            for line in repomap:
                line = line.split("\t")
                # If the dataset corresponds with the target append its features to the list
                if(line[1] == target):
                    file_list.append([line[0], line[1], line[2], line[3]])
        # Close the file
        repomap.close()
        # Return the list
        return file_list


def get_extension(dataset):
    """
    Function to return the extension of a given dataset
    """

    # Return the extension corresponding with each dataset
    if(dataset == "Genome"):
        return ".fa"
    elif(dataset == "Variants"):
        return ".vcf"
    elif(dataset == "Annotation"):
        return ".gff"
    else:
        return ".sam"


def extract_dataset(dataset, seqID, region, convertToBam):
    """
    Function to extract an entire dataset (all the files contained in that dataset)
    """

    # Extract a list of files for the dataset of interest [[filename,dataset,directory,linesize],[...]]
    file_list = obtain_file(target=dataset, mode="dataset")
    # Get the file extension
    extension = get_extension(dataset)
    # Check if the dataset exists in the repository.
    if(not os.path.isdir("./" + dataset) or len(file_list) == 0):
        print("***WARNING: No stored data was found for the dataset of interest: " +
              dataset + "\nNow aborting.")
    # Check if the user provided a seqID
    elif(seqID != "0"):
        # Inform the user
        if(region != "0"):
            print("Now reconstructing data contained in the " + dataset +
                  " dataset for sequence " + seqID + " in the region " + region)
        else:
            print("Now reconstructing data contained in the " +
                  dataset + " dataset for sequence " + seqID)
        # Delete the file if it already exists
        if(os.path.isfile("../Extracted_" + dataset + extension)):
            os.remove("../Extracted_" + dataset + extension)
        # Loop through the files contained in the dataset
        for subfile in file_list:
            # Reconstruct the subfile
            reconstruct_dataset(size=int(subfile[3]), directory=subfile[2], output_file="../Extracted_" +
                                dataset + extension, mode=subfile[1], seqID=seqID, region=region, filemode="a", convertToBam=convertToBam)
        # Inform the user
        print("Reconstruction completed.")
    # Otherwise the user wants to reconstruct the entire dataset
    else:
        # Inform the user
        if(region != "0"):
            print("Now reconstructing data contained in the " + dataset +
                  " dataset for sequence " + seqID + " in the region " + region)
        else:
            print("Now reconstructing data contained in the " + dataset +
                  " dataset for sequence " + seqID)  # Delete the file if it already exists
        if(os.path.isfile("../Extracted_" + dataset + extension)):
            os.remove("../Extracted_" + dataset + extension)
        # Loop through the files contained in the dataset
        for subfile in file_list:
            # Reconstruct the subfile
            reconstruct_dataset(size=int(subfile[3]), directory=subfile[2], output_file="../Extracted_" +
                                dataset + extension, mode=subfile[1], seqID="0", region="0", filemode="a", convertToBam=convertToBam)
        # Inform the user
        print("Reconstruction completed.")
