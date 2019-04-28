#!/usr/bin/env python

# PYTHON FUNCTIONS USED FOR FILE PARSING


# Make the imports
import os
import re
import sys
from subprocess import Popen


def parse_cntg(seq, number, map_file, folder, seq_index):
    """
    Define the parse contig function. This function will take a nucleotide sequence,
    the number of the sequence count, the map file, the current folder and the
    sequence index as input parameters. It will create a file for each contig,
    but if the 250.000 nt threshold is reached it will create a folder contaiing equally
    sized files. Additionally it has to return the updated index.
    """

    # Save the input arguments
    current_nucleotides = seq
    i = number
    map = map_file
    current_dir = folder
    index = seq_index
    # If the file is empty, this means the sequence started with Ns (this may occur during scaffold parsing).
    # In this case, do not write a contig file (it would be empty)
    if(len(current_nucleotides) == 0):
        pass
    # If the length of the contig is smaller than 250.000 nt, write a single file
    elif (len(current_nucleotides) < 250000):
        # Write the file
        output_file = open("{}/Sequence{}".format(current_dir, str(i)), "w")
        output_file.write(current_nucleotides)
        # output_file.close()
        # Write in the map file
        map.write("{}/Sequence{}\t{} {}\n".format(current_dir, str(i), str(index),
                                                  str(index + len(current_nucleotides) - 1)))
        index = index + len(current_nucleotides)
    # If this contig is bigger than 250.000 nt, split it into files of 100.000 nt each
    else:
        # Divide the current contig into chunks of 100.000 nt
        chunks = [current_nucleotides[y:y + 100000]
                  for y in range(0, len(current_nucleotides), 100000)]
        # Create a new folder to store those chunks
        current_dir = "{}/Sequence{}".format(current_dir, str(i))
        os.mkdir(current_dir)
        # Create chunk_count
        chunk_count = 0
        # Write each chunk in a separate file
        for chunk in chunks:
            chunk_count = chunk_count + 1
            # Write the file
            output_file = open(
                "{}/Sequence{}".format(current_dir, str(chunk_count)), "w")
            output_file.write(chunk)
            # output_file.close()
            # Write in the map file
            map.write("{}/Sequence{}\t{} {}\n".format(current_dir, str(chunk_count),
                                                      str(index), str(index + len(chunk) - 1)))
            index = index + len(chunk)
    # Return the updated index
    return index



def parse_scf(seq, folder, map_file, seq_index):
    """
    Define the parse scaffold function. This function will take a nucleotide sequence, the map file,
    the current folder and the sequence index as input parameters.
    It will create a folder for each scaffold, and pass each contig to the parse_cntg function.
    It has to return an updated index.
    """

    # Save the arguments
    current_nucleotides = seq
    current_dir = folder
    map = map_file
    # Store all the contigs in one list
    contig_list = re.split("N+", current_nucleotides)
    # Store all the gaps in one list and remove any empty element
    gap_list = re.split("[^N]+", current_nucleotides)
    gap_list = test = filter(None, gap_list)
    # Initiate the index, the count of subsequences and the current subsequence
    subseq_count = 0
    index = seq_index
    current_subsequence = ""
    # Loop through all the contigs in the list
    for x in range(0, len(contig_list)):
        # Add the current contig into the current subsequence
        current_subsequence = current_subsequence + contig_list[x]
        # If this is the last contig, write the current subsequence into a file
        if(x == len(contig_list) - 1):
            subseq_count = subseq_count + 1
            index = parse_cntg(seq=current_subsequence, number=subseq_count,
                               map_file=map, folder=current_dir, seq_index=index)
        # If the next gap is smaller than 25, iterate to the next contig. DONT FORGET TO ADD THOSE Ns!!!
        elif(len(gap_list[x]) < 25):
            current_subsequence = current_subsequence + gap_list[x]
            continue
        # Otherwise store the current subsequence in a new file
        else:
            subseq_count = subseq_count + 1
            index = parse_cntg(seq=current_subsequence, number=subseq_count,
                               map_file=map, folder=current_dir, seq_index=index)
            map.write("N\t" + str(index) + " "
                      + str(index + len(gap_list[x]) - 1) + "\n")
            index = index + len(gap_list[x])
            current_subsequence = ""
    # Return the updated index
    return index


def parse_super(seq_list, gap_list, folder, map_file):
    """
    Create a method to parse the super subsequences. This function will take a nucleotide sequence,
    the map file, the current folder and the sequence index as input parameters.
    It will create a folder for each supersubseq, and pass each scaffold or contig to their
    respective parsing function. It has to return an updated index.
    """

    # Save the arguments
    super_subseq_list = seq_list
    super_gaps_list = gap_list
    current_dir = folder
    map = map_file
    # Subsequence count is 0 and index is 1
    subseq_count = 0
    index = 1
    # For each of the super subsequences, determine their class and act accordingly
    for x in range(0, len(super_subseq_list)):
        subseq_count = subseq_count + 1
        # If the super subsequence is empty (this will occur when the sequence starts or ends with Ns) dont do anything
        # only update the map if necessary
        if(len(super_subseq_list[x]) == 0):
            # Unless this empty supergap is at the end (this would mean the sequence ends with Ns), update the map file
            if(x != len(super_gaps_list)):
                map.write("N\t{} {}\n".format(str(index), str(
                    index + len(super_gaps_list[x]) - 1)))
                index = index + len(super_gaps_list[x])
        # If this sequence is a contig, store its content in a new file and write the roadmap
        elif(IsContig(super_subseq_list[x])):
            # Parse the contig
            index = parse_cntg(
                seq=super_subseq_list[x], number=subseq_count, map_file=map, folder=current_dir, seq_index=index)
            # Only if this is not the last super subsequence, write in the map file the Ns
            if(x != len(super_gaps_list)):
                # Write in the map file
                map.write("N\t{} {}\n".format(str(index), str(
                    index + len(super_gaps_list[x]) - 1)))
                index = index + len(super_gaps_list[x])
        # If this sequence is a scaffold, store its contigs in a new folder
        else:
            # Create a new folder to allocate the scaffold
            current_subfolder = "{}/Sequence{}".format(
                current_dir, str(subseq_count))
            os.mkdir(current_subfolder)
            # Parse the scaffold
            index = parse_scf(
                seq=super_subseq_list[x], folder=current_subfolder, map_file=map, seq_index=index)
            # Only if this is not the last super subsequence, write in the map file the N gap
            if(x != len(super_gaps_list)):
                # Write the map file
                map.write("N\t{} {}\n".format(str(index), str(
                    index + len(super_gaps_list[x]) - 1)))
                index = index + len(super_gaps_list[x])

def parse_nucleotides(seq, number, map_file):
    """
    Create a parse nucleotides function. This function parses the nucleotides of a sequence using parse_cntg,
    parse_scf and parse_super.
    """

    # Store the arguments
    current_nucleotides = seq
    i = number
    map = map_file
    # If this sequence is a contig, store its content in a new file and write the roadmap
    if(IsContig(current_nucleotides)):
        # Write the info about the current sequence in the map file
        map.write(">Sequence" + str(i) + "\t"
                  + str(len(current_nucleotides)) + "\n")
        # Parse the contig
        parse_cntg(seq=current_nucleotides, number=i,
                   map_file=map, folder=".", seq_index=1)
    # If this sequence is a scaffold, store its contigs in a new folder
    elif(IsScaffold(current_nucleotides)):
        # Create a directory and write the map
        current_dir = "./Sequence" + str(i)
        os.mkdir(current_dir)
        map.write(">Sequence" + str(i) + "\t"
                  + str(len(current_nucleotides)) + "\n")
        # Parse the scaffold
        parse_scf(seq=current_nucleotides, folder=current_dir,
                  map_file=map, seq_index=1)
    # If it is not a contig or a scaffold, it is a chromosome
    else:
        # Create a directory and write the map
        current_dir = "./Sequence" + str(i)
        os.mkdir(current_dir)
        map.write(">Sequence" + str(i) + "\t"
                  + str(len(current_nucleotides)) + "\n")
        # Store all the subsequences separated by N-10000+ gaps in a list
        pattern = "N" * 10000
        super_subseq_list = re.split(pattern + "+", current_nucleotides)
        # Store all the N-10000+ gaps in one list, remove empty elements and those which are lower than 10.000
        super_gaps_list = re.split("[^N]+", current_nucleotides)
        super_gaps_list = filter(None, super_gaps_list)
        super_gaps_list = [g for g in super_gaps_list if len(g) >= 10000]
        # If there are no super gaps, the whole sequence is like a large scaffold
        if(len(super_gaps_list) == 0):
            # All the sequence is a large scaffold
            parse_scf(seq=current_nucleotides, folder=current_dir,
                      map_file=map, seq_index=1)
        # Otherwise parse each super subsequence
        else:
            # Parse each of the super subsequences
            parse_super(seq_list=super_subseq_list,
                        gap_list=super_gaps_list, folder=current_dir, map_file=map)

def IsContig(seq):
    """
    Define a the contig and the scaffold classification function.
    They will return true when one of these sequences is detected
    """

    # Check if it has gaps bigger than N-25, in which case this is a contig
    if(seq.find("NNNNNNNNNNNNNNNNNNNNNNNNN") == -1):
        return True
    else:
        return False


def IsScaffold(seq):
    # If it has N-25 gaps and its length is lower than 1.000.000, this is a scaffold
    if(len(seq) < 1000000):
        return True
    else:
        return False


def parse_dependent_dict(DependentDict):
    """
    Create a function parse dependent file to parse a the information of a dependent file contained
    in a dependent dictionary into the repository
    """

    # Initiate required variables
    seq_count = 0
    # Create a SeqID_Map file. This map file stores information about the seqID, the location of files and the count of sequences
    seqID_map = open("./SeqIDs_Map.txt", "w")
    # Loop through the sequences in the dependent file
    for key in DependentDict.keys():
        seq_count += 1
        # Create subfiles for 1MB regions
        # Write in the seqID map
        seqID_map.write(">" + key + "\t/Sequence" + str(seq_count) + "\n")
        # Create a folder for the sequence and open the first subfile
        os.mkdir("./Sequence" + str(seq_count))
        # Loop through the regions stored in the dictionary {seqID:{1M:[line[1:],],2M:[],}}
        for region in DependentDict[key].keys():
            # Write in the map
            seqID_map.write("/Sequence" + str(seq_count)
                            + "/Region" + str(region) + "\n")
            # Open a new file for the region
            output_file = open("./Sequence" + str(seq_count)
                               + "/Region" + str(region), "w")
            # Write the entries of the region into the file
            for line in DependentDict[key][region]:
                output_file.write(line)
            # Close the file
            output_file.close()
    # Close seqID map
    seqID_map.close()


def format_bytes(size):
    # 2**10 = 1024
    power = 2**10
    n = 0
    power_labels = {0 : '', 1: 'K', 2: 'M', 3: 'G', 4: 'T'}
    while(size > power):
        size /= power
        n += 1
    return str(size), power_labels[n]+'B'


def parse_dataset(dataset, input_path, size, update):
    """
    Create a parse dataset function to parse a given dataset in a file into the repository structure
    """

    file_size = os.path.getsize(input_path)
    size, label = format_bytes(file_size)
    size = str(size)

    # PARSE A GENOME DATASET

    # If the file to be parsed is a genome file.
    if(dataset == "Genome"):
        # Load the input file
        input_file = open(input_path, "r")
        # Initiate necessary variables
        i = 0
        current_nucleotides = ""
        mapfile = open("Map.txt", "w")
        invalid_input = True
        # Loop through the lines of the file
        for line in input_file:
            # Only read the line if it does not start with a # (description lines can start with this)
            if(line[0] != "#"):
                # Remove all new lines (sequence stored as a big line)
                line = line.rstrip()
                # If the line starts with that symbol, this is a new sequence
                if (line[0] == ">"):
                    # The input file is not invalid anymore
                    invalid_input = False
                    # Unless this is the first sequence, write a file with its nucleotides
                    if(current_nucleotides != ""):
                        parse_nucleotides(
                            seq=current_nucleotides, number=i, map_file=mapfile)
                    # Add one to the sequence count
                    i = i + 1
                    # Empty the current nucleotides
                    current_nucleotides = ""
                # Otherwise store the current sequence in the current_nucleotides variable
                else:
                    if(current_nucleotides == ""):
                        line_size = str(len(line))
                    current_nucleotides = current_nucleotides + line

        # If when the parsing is finished, the input file has not been validated, raise an exception
        if(invalid_input):
            # Inform the user
            print(
                "\n\n\t***INPUT ERROR: Please make sure to input a valid FASTA format file***\n")
            # Delete the invalid Genome folder
            os.system("rm -r ../Genome")
            # Exit the system
            sys.exit()

        # Parse the last sequence left out of the loop
        parse_nucleotides(seq=current_nucleotides, number=i, map_file=mapfile)
        # Close the map file and the input file
        input_file.close()
        mapfile.close()
        # If it is an update it is required to update the repomap
        if(update):
            # Open the repository map file in read mode and append all the lines to a list
            repo_list = []
            with open("../RepoMap.txt", "r") as repomap:
                # Loop through the lines
                for line in repomap:
                    # Split the line
                    line = line.split("\t")
                    # If the dataset of the current file is genome, substitute the fields of the name, line size and size of the file
                    if("Genome" == line[1]):
                        line[0] = os.path.basename(input_path)
                        line[3] = line_size
                        line[4] = "{} {}\n".format(size, label)
                    # Store the line with or without modification
                    repo_list.append("\t".join(line))
            # Close the file
            repomap.close()
            # Open the file with write permit and write all the elements of the repo list
            with open("../RepoMap.txt", "w") as repomap:
                for line in repo_list:
                    repomap.write(line)
            # Close the file
            repomap.close()
        # Otherwise simply append a new line to the repomap
        else:
            # Open the repository map file in append mode
            with open("../RepoMap.txt", "a") as repomap:
                # Write the information corresponding with the added file
                repomap.write(os.path.basename(input_path) + "\t" + dataset
                              + "\t" + "./" + dataset + "\t" + line_size + "\t" + "{} {}\n".format(size, label))
            # Close the file
            repomap.close()

    # PARSE A ANNOTATION DATASET

    # If it is an annotation or variants dataset
    elif(dataset == "Annotation"):
        # Open the dependent file and read the first line
        annotation_file = open(input_path, "r")
        # Create a comment file to store the comments at the begining of the file
        comment_file = open("./Comments.txt", "w")
        # Initiate required variables
        annotation_dict = {}  # {seqID:{1M:[line[1:],],2M:[],}}
        line_number = 0
        # Loop through the annotation file lines
        for line in annotation_file:
            line_number += 1
            # If the line starts with # or @; it corresponds with a comment, store it in the comment file
            if(line[0] == "#"):
                comment_file.write(line)
            # Otherwise the line corresponds with an entry
            else:
                # Split the line using \t
                line = line.split("\t")
                # Annotation should have at least 8 fields and the fields 3 and 4 should be a number
                if(len(line) < 8 or not line[3].isdigit() or not line[4].isdigit()):
                    # Inform the user
                    print(
                        "\n\n\t***INPUT ERROR: Invalid input file. Please make sure to provide a valid GFF file for your Annotation dataset.***")
                    print("\n\t The problem was raised in line number "
                          + str(line_number) + " :\n\n\t" + "\t".join(line))
                    os.system("rm -r ../$(basename " + input_path + ")")
                    # Close the files and exit
                    comment_file.close()
                    annotation_file.close()
                    sys.exit()

                # If the seqId is already a key in the dictionary, append the new line to the list (but without the seqID)
                elif (line[0] in annotation_dict.keys()):
                    # Determine the region of the current line using the start index
                    region = (int(line[3]) / 1000000)
                    # Check if the region already has a list, if it does, append the new line to the list
                    if(region in annotation_dict[line[0]].keys()):
                        annotation_dict[line[0]][region].append(
                            "\t".join(line))
                    # Otherwise create a new list with the line
                    else:
                        annotation_dict[line[0]][region] = ["\t".join(line)]

                # If this seqId is not in the dictionary, create a new list and add the current line, and update the key count
                # This key count will help to know the order of the sequences when reconstructing the file (dictionaries are not sorted, so
                # its entries will be shown in random order).
                else:
                    annotation_dict[line[0]] = {}
                    # Determine the region of the current line using the start index
                    region = (int(line[3]) / 1000000)
                    annotation_dict[line[0]][region] = ["\t".join(line)]
        # Close the files
        comment_file.close()
        annotation_file.close()
        # Use the parse dependent dict function to parse the information contained in the resulting dictionary into the repo
        parse_dependent_dict(DependentDict=annotation_dict)
        # If it is an update it is required to update the repomap
        if(update):
            # Open the repository map file in read mode and append all the lines to a list
            repo_list = []
            with open("../../RepoMap.txt", "r") as repomap:
                # Loop through the lines
                for line in repomap:
                    # Split the line
                    line = line.split("\t")
                    # If the current file is has the same name, substitute the field of the size of the file with the new one
                    if(os.path.basename(input_path) == line[0]):
                        line[4] = "{} {}\n".format(size, label)
                    # Store the line with or without modification
                    repo_list.append("\t".join(line))
            # Close the file
            repomap.close()
            # Open the file with write permit and write all the elements of the repo list
            with open("../../RepoMap.txt", "w") as repomap:
                for line in repo_list:
                    repomap.write(line)
            # Close the file
            repomap.close()
        # Otherwise simply append a new line to the repomap
        else:
            # Open the repository map file in append mode
            with open("../../RepoMap.txt", "a") as repomap:
                # Write the information corresponding with the added file
                repomap.write(os.path.basename(input_path) + "\t" + dataset + "\t" + "./"
                              + dataset + "/" + os.path.basename(input_path) + "\t1\t" + "{} {}\n".format(size, label))
            # Close the file
            repomap.close()

    # PARSE A VARIANTS DATASET

    elif(dataset == "Variants"):
        # Open the dependent file and read the first line
        variants_file = open(input_path, "r")
        # Create a comment file to store the comments at the begining of the file
        comment_file = open("./Comments.txt", "w")
        # Initiate required variables
        variants_dict = {}  # {seqID:{1M:[line[1:],],2M:[],}}
        line_number = 0
        # Loop through the file lines
        for line in variants_file:
            line_number += 1
            # If the line starts with # or @; it corresponds with a comment, store it in the comment file
            if(line[0] == "#"):
                comment_file.write(line)
            # Otherwise the line corresponds with an entry
            else:
                # Split the line using \t
                line = line.split("\t")
                # Variants should have at least 8 fields and the field 1 should be a number
                if(len(line) < 8 or not line[1].isdigit()):
                    # Inform the user
                    print(
                        "\n\n\t***INPUT ERROR: Invalid input file. Please make sure to provide a valid VCF file for your Variants dataset.***")
                    print("\n\t The problem was raised in line number "
                          + str(line_number) + " :\n\n\t" + "\t".join(line))
                    os.system("rm -r ../$(basename " + input_path + ")")
                    # Close the files and exit
                    variants_file.close()
                    comment_file.close()
                    sys.exit()
                # If the seqId is already a key in the dictionary, append the new line to the list (but without the seqID)
                elif (line[0] in variants_dict.keys()):
                    # Determine the region of the current line using the start index
                    region = (int(line[1]) / 1000000)
                    # Check if the region already has a list, if it does, append the new line to the list
                    if(region in variants_dict[line[0]].keys()):
                        variants_dict[line[0]][region].append("\t".join(line))
                    # Otherwise create a new list with the line
                    else:
                        variants_dict[line[0]][region] = ["\t".join(line)]

                # If this seqId is not in the dictionary, create a new list and add the current line, and update the key count
                # This key count will help to know the order of the sequences when reconstructing the file
                # (dictionaries are not sorted, so its entries will be shown in random order).
                else:
                    variants_dict[line[0]] = {}
                    # Determine the region of the current line using the start index
                    region = (int(line[1]) / 1000000)
                    # Add the line to its corresponding region
                    variants_dict[line[0]][region] = ["\t".join(line)]
        # Close the files
        comment_file.close()
        variants_file.close()
        # Use the parse dependent dict function to parse the information
        # contained in the resulting dictionary into the repo
        parse_dependent_dict(DependentDict=variants_dict)
        # If this is an update it is required to update the repomap
        if(update):
            # Open the repository map file in read mode and append all the lines to a list
            repo_list = []
            with open("../../RepoMap.txt", "r") as repomap:
                # Loop through the lines
                for line in repomap:
                    # Split the line
                    line = line.split("\t")
                    # If the current file is has the same name,
                    # substitute the field of the size of the file with the new one
                    if(os.path.basename(input_path) == line[0]):
                        line[4] = "{} {}\n".format(size, label)
                    # Store the line with or without modification
                    repo_list.append("\t".join(line))
            # Close the file
            repomap.close()
            # Open the file with write permit and write all the elements of the repo list
            with open("../../RepoMap.txt", "w") as repomap:
                for line in repo_list:
                    repomap.write(line)
            # Close the file
            repomap.close()
        # Otherwise simply append a new line to the repomap
        else:
            # Open the repository map file in append mode
            with open("../../RepoMap.txt", "a") as repomap:
                # Write the information corresponding with the added file
                repomap.write(os.path.basename(input_path) + "\t" + dataset + "\t" + "./"
                              + dataset + "/" + os.path.basename(input_path) + "\t1\t" + "{} {}\n".format(size, label))
            # Close the file
            repomap.close()

    #	PARSE AN ALIGNMENT DATASET

    elif(dataset == "Alignment"):
        # If this is an update it is required to update the repomap

        # Convert BAM to SAM if the input file is BAM
        if ('.bam' in input_path):
            samOutput_path = './' + os.path.basename(input_path)[:-3] + 'sam'
            print("Converting BAM to SAM.")
            Popen("samtools view -h -o {} {}".format(samOutput_path,
                                                     input_path), shell=True).wait()
            print("Conversion done.")
            # New path for this file is .sam instead of .bam
            input_path = samOutput_path

            # Change the size of the file to SAM's size
            size = str(os.path.getsize(input_path) / (1024 * 1024))

        if(update):
            # Open the repository map file in read mode and append all the lines to a list
            repo_list = []
            repomap = open("../../RepoMap.txt", "r")
            # Loop through the lines
            for line in repomap:
                # Split the line
                line = line.split("\t")
                # If the current file is has the same name,
                # substitute the field of the size of the file with the new one
                if(os.path.basename(input_path) == line[0]):
                    line[4] = "{} {}\n".format(size, label)
                # Store the line with or without modification
                repo_list.append("\t".join(line))
            # Close the file
            repomap.close()
            # Open the file with write permit and write all the elements of the repo list
            repomap = open("../../RepoMap.txt", "w")
            for line in repo_list:
                repomap.write(line)
        # Otherwise simply append a new line to the repomap
        else:
            # Open the repository map file in append mode
            repomap = open("../../RepoMap.txt", "a")
            # Write the information corresponding with the added file
            repomap.write(os.path.basename(input_path) + "\t" + dataset + "\t" + "./"
                          + dataset + "/" + os.path.basename(input_path) + "\t1\t" + "{} {}\n".format(size, label))
        # Create a comment file to store the comments at the begining of the file
        comment_file = open("./Comments.txt", "w")
        # Open the sam file
        sam_file = open(input_path, "r")
        # Inititate varibles
        seqID = ""
        seqID_dict = {}  # {seqID:directory}
        seq_count = 0
        # Loop through the lines of the file and store the comments
        for line in sam_file:
            # If the line is a comment, write in the comment file
            if(line[0] == "@"):
                comment_file.write(line)
            # Otherwise it is a read, store it
            else:
                # Split the line
                line = line.split("\t")
                # If the seqID is new (there is no directory already in place for the sequence)
                if(line[2] not in seqID_dict.keys()):
                    # Add one to the seq_count
                    seq_count += 1
                    # Add it to the dictionary, with its count
                    seqID_dict[line[2]] = "Sequence" + str(seq_count)
                    # Create a new directory for the sequence
                    os.mkdir("./Sequence" + str(seq_count))
                    # Determine the region of the line
                    region = (int(line[3]) / 1000000)
                    # Append the line to the file of the corresponding region
                    repomap.close()
                    repomap = open("./Sequence" + str(seq_count)
                                   + "/Region" + str(region), "a")
                    repomap.write("\t".join(line))
                # If the sequence ID is in the dictionary, no need to create a new directory
                else:
                    # Determine the region of the line
                    region = (int(line[3]) / 1000000)
                    # If the region and seqID has not changed since the last line, do not close and open a new file
                    if(repomap.name == seqID_dict[line[2]] + "/Region" + str(region)):
                        repomap.write("\t".join(line))
                    # Otherwise it is required to open a new file
                    else:
                        # Append the line to the file of the corresponding region
                        repomap.close()
                        repomap = open(
                            seqID_dict[line[2]] + "/Region" + str(region), "a")
                        repomap.write("\t".join(line))
        # Close the files
        repomap.close()
        sam_file.close()
        comment_file.close()
        # Create a map file to inform the structure of the dataset
        mapfile = open("./SeqIDs_Map.txt", "w")
        # Include all the sequences with their subfiles
        for seqID in seqID_dict.keys():
            # Write the seqID and its directory
            mapfile.write(">" + seqID + "\t/" + seqID_dict[seqID] + "\n")
            # Create a list of subfiles and loop through them. Write them in the map file.
            # Sort the files according to their region.
            region_list = sorted(os.listdir(
                "./" + seqID_dict[seqID]), key=lambda x: x[6:])
            for subfile in region_list:
                mapfile.write("/" + seqID_dict[seqID] + "/" + subfile + "\n")
        # Close the map file
        mapfile.close()

        # Remove the SAM file that is a result of conversion
        if(os.path.isfile("./converted.txt")):
            os.remove(input_path)

    # RAISE AN ERROR

    # Otherwise there was a problem
    else:
        print("\n***INTERNAL ERROR: DATASET NOT RECOGNISED: " + dataset + "\n")

    # Exit the function
    return



def parse_delta_file(file_path):
    """
    Create function parse delta file to store the information contained in the mummer delta file
    into a list with an entry for every alignment [oldID, newID, oldLength, newLength]
    """

    delta_list = []
    # Loop through the lines of the delta file
    with open(file_path, "r") as delta:
        for line in delta:
            if(line[0] == ">"):
                line = line[1:].split()
                delta_list.append(line[:4])
    # Return the list
    return delta_list


def parse_snp_file(file_path):
    """
    Create a parse_snp_file function to parse a snp file into a dictionary containing
    the number of insertions, deletions and substitutions for each alignment
    """

    # Initiate an empty dictionary {(oldID, newID):[insertions,deletions,substitutions]...}
    snp_dict = {}
    # Open the snp file and loop through its lines
    with open(file_path, "r") as snp:
        for line in snp:
            # Split the line into its fields
            line = line.rstrip().split("\t")
            key = oldID, newID = line[10], line[11]
            # If the oldID is already in the dictionary, add the current snp into the count
            if(key in snp_dict.keys()):
                # If it is an insertion, add one to the count
                if(line[1] == "."):
                    snp_dict[key][0] += 1
                # If it is a deletion, add one to the count
                elif(line[2] == "."):
                    snp_dict[key][1] += 1
                # Otherwise it is a substitution, add one to the count
                else:
                    snp_dict[key][2] += 1
            # If the current oldID is not a key in the dictionary, create a new entry
            else:
                # If it is an insertion, add one to the count
                if(line[1] == "."):
                    snp_dict[key] = [1, 0, 0]
                # If it is a deletion, add one to the count
                elif(line[2] == "."):
                    snp_dict[key] = [0, 1, 0]
                # Otherwise it is a substitution, add one to the count
                else:
                    snp_dict[key] = [0, 0, 1]

    # Return the dictionary
    return snp_dict


def parse_coords_file(file_path):
    """
    Create a parse_coords_file function to return a dictionary
    with the information contained in the summary coords file
    """

    # Create the dictionary {oldID:[lines]}
    coords_dict = {}
    # Open the file and loop through the lines
    with open(file_path, "r") as coords:
        for line in coords:
            line = line.rstrip().split("\t")
            key = (line[9], line[10])
            # If the oldId is already in the dictionary, add the line to the current list
            if(key in coords_dict.keys()):
                coords_dict[key].append(line[:9])
            # Otherwise create a new entry in the dictionary
            else:
                coords_dict[key] = [line[:9]]

    # Return the dictionary
    return coords_dict

# # Create lowest_highest function that will be used for parsing annotation and variants datasets.
# # It returns the highest and the lowest indexes in the output file.
# def lowest_highest(output_file, mode):
#     # When the file is completely written, is necessary to write in the seqID map file the information
#     # Initiate lowest and highest values (indicate the lowest and the highest index inside the file)
#     lowest = 99999999999999999999999999999999
#     highest = 0
#     # If it is a variants file, determine if its index is higher than highest or lower than lowest, and substitute it if so.
#     if(mode == "Variants"):
#         # Loop through the lines of the file
#         for line in output_file:
#             if(int(line.split("\t")[0]) > highest):
#                 highest = int(line.split("\t")[0])
#             if(int(line.split("\t")[0]) < lowest):
#                 lowest = int(line.split("\t")[0])
#     # If it is an alignment file, the index is in the field 3
#     if(mode == "Alignment"):
#         # Loop through the lines of the file
#         for line in output_file:
#             if(int(line.split("\t")[3]) > highest):
#                 highest = int(line.split("\t")[3])
#             if(int(line.split("\t")[3]) < lowest):
#                 lowest = int(line.split("\t")[3])
#     # Same for annotation
#     else:
#         # Loop through the lines of the file
#         for line in output_file:
#             if(int(line.split("\t")[3]) > highest):
#                 highest = int(line.split("\t")[3])
#             if(int(line.split("\t")[2]) < lowest):
#                 lowest = int(line.split("\t")[2])
#     # Return results
#     return [str(lowest), str(highest)]