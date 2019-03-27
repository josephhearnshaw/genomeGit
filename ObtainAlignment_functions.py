#!/usr/bin/env python

# PYTHON FUNCTIONS REQUIRED IN THE REPOSITORY AUTO-UPDATES PART I: ALIGNMENT OBTENTION.
from pyfaidx import Fasta
from subprocess import Popen, call
import multiprocessing
import shutil
import os
import math as ma
import datetime
import signal
import re
import hashlib
import sys
from threading import Thread
from multiprocessing import Process
from StoreAlignment import obtain_alignment_pickle

# Create global variables. They will be used for the threads
tabix_queries = {}
alignment_pickle = ""


def obtain_SHA1(input_string):
    """ This function calculates the SHA-1 hash of the input sequneces
    It requires the sequence obtained from the assembly dictionary """
    # Calculate the SHA-1 hash of both input sequences
    return str(hashlib.sha1(input_string).hexdigest())


def parse_assembly(assembly_file):
    """ This function parses the assembly charactersitcis into a dictionary, requiring
    the assembly files. It'll produce a dictionary with caharacteris as follows:
    {seqKey:[oldID,[line,line...],hash,hash-1]}  """

    # Create output dictionary
    assembly_dictionary = {}
    # Open the asembly file
    with open(assembly_file, "r") as input_file:
        # Loop through the lines
        sequence_count = 0
        for line in input_file:
            line = line.rstrip()
            # If the line starts with >, it is a new sequence
            if(line[0] == ">"):
                # The key is the seqID - need to split and only get the seq ID to avoid errors
                current_key = line.split()[0][1:]
                # Add the first two values (the seqID and an empty string for the sequence)
                assembly_dictionary[current_key] = [line[1:], []]
            # Otherwise the line corresponds with the nucleotides of the sequence
            else:
                assembly_dictionary[current_key][1].append(line)
        # After looping through all the lines of the file, determine the SHA-1 of the sequences in the dictionary
        for key in assembly_dictionary:
            # Forward SHA1
            assembly_dictionary[key].append(obtain_SHA1("".join(assembly_dictionary[key][1])))
            # Reverse all the SHA1
            assembly_dictionary[key].append(obtain_SHA1("".join(assembly_dictionary[key][1])[::-1]))
    # Return the assembly dictionary
    return assembly_dictionary


def parse_snp_file(file_path):
    """ Parse all the SNP data into a dictionary as {SeqID: "oldIndex oldChar newChar newIndex", "...",} """
    # Initate variables
    snp_dict = {}
    # Open the file and loop through the lines
    with open(file_path, "r") as snp_file:
        for line in snp_file:
            # Split the line fields
            line = line.split("\t")
            # If the seqID is already a key in the dictionary, append the new line to the list
            oldIndex, oldChar, newChar, newIndex = line[:4]
            snp_seqID = line[10]
            if(snp_seqID in snp_dict.keys()):
                snp_dict[snp_seqID].append(" ".join([oldIndex, oldChar, newChar, newIndex, snp_seqID]))
            # Otherwise create a new list
            else:
                snp_dict[snp_seqID] = [" ".join([oldIndex, oldChar, newChar, newIndex])]

    return snp_dict


def calculate_MashMap_score(Blocks, lengthA, lengthB):
    """ This function calculates the score, as described further in documentation.
    It requires each block of the output, and the lengths of the query and reference blocks.
    A resultant score is produced, which is used later to filter data to find which alignments
    should be used for further analysis. """

    # Initiate the coherence and block count
    coherence = 0
    block_number = 1
    # Calculate length factor
    length_factor = 2 - (float(lengthA) / float(lengthB) +
                         float(lengthB) / float(lengthA))
    # Save the first block, named as prevblock
    prevblock = [int(x) for x in Blocks[0]]
    # Declare all elements of the prevblock
    prev_query_start, prev_query_end, prev_ref_start, prev_ref_end = prevblock

    for current_block in Blocks[1:]:
        # initate the current block of nucleotides
        current_block = [int(x) for x in current_block]
        # Declare all elements of the current block
        query_start, query_end, ref_start, ref_end = current_block
        # The current query starts are sorted in ascending order; no need to check their order
        # If current query start pos is greater than previous start pos, add one.
        # Otherwise deduct 1 from coherence
        if(query_start > prev_query_start and query_end > prev_query_end):
            coherence += 1
        elif(query_start > prev_query_start and query_end < prev_query_end):
            coherence -= 1
        # Make prevblock equal to the previous block
        prevblock = current_block
    # Caclulate the score - see documentation for further info on score
    score = float(coherence) / block_number + length_factor
    # Return the score
    return score


def parse_mashMap_line(line):
    """ Split each line and parse each index into its respective elements i.e. query = line[0] """
    line = line.split()
    query, query_length, query_start, query_end, ref, ref_length, ref_start, ref_end = \
        line[0], line[1], line[2], line[3], line[5], line[6], line[7], line[8]
    # Return as a nested list
    return((query, query_length), (ref, ref_length), [query_start, query_end,
                                                      ref_start, ref_end])


def parse_mashMap_output_to_Dict(input_file):
    """ This function parses the mashmap output into a dictionary. It uses a query key and a refernece key;
these keys correspodn to the query and reference names and their associated lengths, respectively. Each
    key will contain its respective block. The function will check to see if there is a new alignment. If
    there is, it'll add the new keys."""

    # Declare the unfiltered dictionary
    Unfiltered_MashDict = {}  # {[blocks and modifications]}
    # Read the first line
    line = input_file.readline().rstrip()
    # Initiate query and ref key, and blocks
    query_key, ref_key, block = parse_mashMap_line(line)
    # Initiate the unfiltered MashMap dictionary
    Unfiltered_MashDict[query_key] = {ref_key: [block]}
    # Store the first query and ref key as previous query and ref keys
    prev_query_key, prev_ref_key = query_key, ref_key
    # Read through the input file
    for line in input_file:
        # Strip the line of whitespace
        line = line.rstrip()
        # Add the respective elements in from the file
        query_key, ref_key, block = parse_mashMap_line(line)
        # Check if the previous query and ref are equal to eachother, if not it's a new block
        if query_key == prev_query_key:
            if ref_key == prev_ref_key:
                Unfiltered_MashDict[query_key][ref_key].append(block)
            else:
                Unfiltered_MashDict[query_key][ref_key] = [block]
        else:
            Unfiltered_MashDict[query_key] = {ref_key: [block]}
        # Reassign the previous keys so it lags by one
        prev_query_key, prev_ref_key = query_key, ref_key

    return Unfiltered_MashDict


def find_best_mashmap_alignment(query_length, alignment_dict):
    """This function takes the query length and the alignment dictionary. It will then
find the best alignments and return the reference key and score associated with those
alignments

This is for MashMap"""

    # Obtain the keys
    refs = list(alignment_dict.keys())
    max_ref_key = refs[0]
    blocklist = alignment_dict[max_ref_key]
    ref_length = max_ref_key[1]
    # Obtain the score
    max_score = calculate_MashMap_score(blocklist, query_length, ref_length)
    # Iterate through the keys
    for ref_key in refs[1:]:
        # Obtain the blocklist associated with each key
        blocklist = alignment_dict[ref_key]
        ref_length = ref_key[1]
        score = calculate_MashMap_score(blocklist, query_length, ref_length)
        # If the score is greater than the current max score, then reassign the max store
        # as the current score and the associated ref key as the max ref key - return the best
        if score > max_score:
            max_score = score
            max_ref_key = ref_key
    return max_ref_key, max_score


def retrieve_fasta_sequence_and_make_fifo(pyfaidx_obj, ID, fifo_path):
    """ This method pipes the sequences as fifos, which are thread safe
and memory efficient; improving the speed of reading sequences and using less memory.
    No writing is performed to long-term memory with this method."""

    # Obtain the sequence using Pyfadix (used as opposed to samtools due to being faster)
    seq = ">{}\n{}".format(ID, pyfaidx_obj[ID][:].seq)
    # Make the temporary fifo file-path
    os.mkfifo(fifo_path)
    # Create the fork processes (original with PID of child and forked with PID of 0)
    pid = os.fork()
    # If the process ID is forked it'll equal 0, so get the PID for that PID
    if pid == 0:
        pid = os.getpid()
        # Write the sequence to the memory and call it 'fifo'
        with open(fifo_path, 'w') as fifo:
            fifo.write(seq)
        # Remove the temporary directory
        os.remove(fifo_path)
        # Terminate all Processes
        os.kill(pid, signal.SIGTERM)


def filter_mashMap(mash_file, directory, alignment_pickle, out_file, threads, inversion_flag):
    """ This function Produces Fasta files for the best alignments and passes them
to nucmer for alignment."""
    # Obtain Reference and query fasta files
    reference_directory = "{}/Compare_OldAssembly.fa".format(directory)
    query_directory = "{}/Compare_NewAssembly.fa".format(directory)
    fasta_query = Fasta(query_directory)
    fasta_ref = Fasta(reference_directory)
# Create a mashMap dictionary out of the contents of the unfiltered MashMap file
    with open(mash_file, "r") as input_file:
        Unfiltered_MashDict = parse_mashMap_output_to_Dict(input_file)

    # Set the initital number of process wrokers equal to none
    num_of_workers = None
    # Obtain the thread pool
    tp = multiprocessing.Pool(num_of_workers)
    # Set count to be 0
    count = 0
    # Create temporary directories
    os.mkdir('tmp_delta')
    # Get the number of alignments
    key_length = len(Unfiltered_MashDict)
    print("\n\n\t\t\t\t\t***MASHMAP FINISHED, NOW RUNNING NUCMER***\n\n")
    # Get the number of zeros required for naming delta files in order (see below Zfill for clarification)
    num_of_zero = ma.floor(ma.log10(key_length)) + 1
    # Loop through the alignments in the dictionary and delete those which do not meet
    # the score as calculated in score function.
    for query_key in Unfiltered_MashDict.keys()[:]:
        query_ID = query_key[0]
        query_length = query_key[1]
        alignments = Unfiltered_MashDict[query_key]

        ref_key, score = find_best_mashmap_alignment(query_length, alignments)
        ref_ID = ref_key[0]
        # Query

        query_ID, ref_ID = ref_ID, query_ID

        query_file_name = 'tmp_delta/{}_{}_query'.format(query_ID, ref_ID)
        retrieve_fasta_sequence_and_make_fifo(fasta_query, query_ID, query_file_name)

        ref_file_name = 'tmp_delta/{}_{}_ref'.format(query_ID, ref_ID)
        retrieve_fasta_sequence_and_make_fifo(fasta_ref, ref_ID, ref_file_name)

        # Create the name for each delta file - Each file must start as 01 ... 010 ... 0100 etc
        # Important for subsequent commands using cat && awk; ensures files are sorted in order
        zfill_count = str(count).zfill(int(num_of_zero))
        name_out = "tmp_delta/{}_{}_{}".format(zfill_count, query_key[0], ref_key[0])
        print("\n\t\t -Running Nucmer now for reference {} and query {} at {}"
              .format(ref_ID, query_ID, str(datetime.datetime.now())))
        # print('Running nucmer for {} and {} at {}\n'.format(query_ID, ref_ID, str(datetime.datetime.now())))
        tp.apply_async(get_sequences_mashmap, (ref_file_name, query_file_name,
                                               name_out, threads, score, inversion_flag))
        # Increase the count by one
        count += 1
    # Close and join the thread pool and delete the temporary directories once done
    tp.close()
    tp.join()
    # Remove old file directories, add the real directories and add the NUCMER using sed && awk
    call("cat tmp_delta/* | awk '/[/]/ || /^NUCMER/ {{ next }} 2' > {} ".format(out_file), shell=True)
    # Add file directory paths for the query and reference to the delta file
    call("sed -i '1i {} {}' {} ".
         format(reference_directory, query_directory, out_file), shell=True)
    # Add NUCMER to the following line
    sed_command = ['sed', '-i', '2iNUCMER', out_file]
    Popen(sed_command).wait()
    # Remove the temporary directory
    shutil.rmtree('tmp_delta')


def get_sequences_mashmap(ref_file_name, query_file_name, name_out, threads, score, inversion_flag):
    """ This function runs the best alignments through nucmer """

    if(inversion_flag == 1):
        command = ['nucmer',  '--mum', ref_file_name,
                   query_file_name, '-p', name_out, '-t', threads]
    # Add score to each delta file at the end of the alignment header line
        Popen(command).wait()
    else:
        command = ['nucmer', str(inversion_flag),  '--mum', ref_file_name,
                   query_file_name, '-p', name_out, '-t', threads]
    # Add score to each delta file at the end of the alignment header line
        Popen(command).wait()

    Popen("sed 's/^>.*/& {}/' {} -i".format(score, name_out + ".delta"), shell=True).wait()


def calculate_nucmer_score(Blocks, length_A, length_B):
    """ This function calculates the score as described in documentation.
    Requires each block and the respective alignment lengths"""

    coherence = 0
    block_number = 1
    # Calculate length factor
    length_factor = 2 - (float(length_A) / float(length_B) +
                         float(length_B) / float(length_A))
    # Get the previous block and assign to start and end positions
    prevblock = [int(x) for x in Blocks[0].split()]
    prev_start, prev_end = prevblock[0], prevblock[2]
    # Loop through the blocks and modifications
    for block in Blocks[1:]:
        # Check if it's an alignment block; not a mod index
        if(" " in block):
            # Get the current block and assign to start and end positions
            current_block = [int(x) for x in block.split()]
            start, end = current_block[0], current_block[2]
            # Add one to the count of fragments
            block_number += 1
            if(start > prev_start and end > prev_end):
                coherence += 1
            elif(start > prev_start and end < prev_end):
                coherence -= 1
            # Reassign the previous block to equal  the current block
            prev_start, prev_end = start, end
    # Determine the score of this alignment
    score = float(coherence) / float(block_number) + length_factor
    # Return the obtained score
    return score


def prase_nucmer_delta_file(input_file):
    """ Parses the delta file into an 'unfiltered delta dict'.
    Requires the input delta file """
    Unfiltered_DeltaDict = {}  # {>alignment:[blocks and modifications]}
    for line in input_file:
        line = line.rstrip()
        # Identify the new alignments and store them as keys
        if (line[0] == ">"):
            current_key = line
            Unfiltered_DeltaDict[current_key] = []
            # Add the corresponding block for this allignment, including the mod index
            # See Mummer Manual regarding mod index.
        else:
            Unfiltered_DeltaDict[current_key].append(line)
    return Unfiltered_DeltaDict


def find_best_nucmer_alignment(Unfiltered_DeltaDict):
    """This function finds the best alignments, based on score, and adds it to the
    filtered dictionary. Requires the unfiltered dictionary """
    Filtered_DeltaDict = {}  # {old_seqID:[>alignment,[blocks and modifications],score]}
    for alignment in Unfiltered_DeltaDict.keys():
        # Obtain all alignments
        all_alignments = alignment

        # Obtain blocks
        Blocks = Unfiltered_DeltaDict[alignment]
        alignment = alignment.split()
        # Obtain the lengths
        length_A, length_B = alignment[2:4]
        ref, query = alignment[0:2]
        ref = ref[1:]
        score = calculate_nucmer_score(Blocks, length_A, length_B)
        # If the score is greater than the current score, then add to the new dictionary
        try:
            old_score = Filtered_DeltaDict[query][2]
            if(score > old_score):
                Filtered_DeltaDict[query] = [all_alignments,
                                             Blocks, score]
        # If it is the first alignment of this old seqID, add it to the filtered dictionary
        except KeyError:
            Filtered_DeltaDict[query] = [all_alignments,
                                         Blocks, score]
    return Filtered_DeltaDict


def write_nucmer_delta(Filtered_DeltaDict, file_directory_line, nucmer_line, out_file):
    """ This function writes the Filtered Dictionary out for further use """
    # Write the contents of Filtered_delta into a new Filtered delta file
    with open(out_file, "w") as output_file:
        # Write the first and second lines (file path to fasta files and NUCMER line)
        output_file.write("{}\n".format(file_directory_line))
        output_file.write("{}\n".format(nucmer_line))
        # Loop through the keys in the Filtered_DeltaDict
        for old_seqID in Filtered_DeltaDict.keys():
            # Write the alignment header. Add an extra field at the end indicating the score.
            output_file.write("{} {} \n".format(Filtered_DeltaDict[old_seqID][0],
                                                str(Filtered_DeltaDict[old_seqID][2])))
            # Write the blocks and modifications
            for BlockModifications in Filtered_DeltaDict[old_seqID][1]:
                output_file.write("{}\n".format(BlockModifications))


# Create a function filter_delta to filter the results in the multi-delta so that it only contains the best alignment of each query sequence
def filter_nucmer_delta(delta_file, out_file):
    """This function filters the delta file according to score, and writes it out.
    Requires a delta-file and out directory"""
    Filtered_DeltaDict = {}  # {old_seqID:[>alignment,[blocks and modifications],score]}
    with open(delta_file, "r") as input_file:
        # Obtaining the file directory and NUCMER line
        file_directory_line = input_file.readline().rstrip()
        nucmer_line = input_file.readline().rstrip()
        Unfiltered_DeltaDict = prase_nucmer_delta_file(input_file)
    Filtered_DeltaDict = find_best_nucmer_alignment(Unfiltered_DeltaDict)
    write_nucmer_delta(Filtered_DeltaDict, file_directory_line, nucmer_line, out_file)


def multiFasta_construct(directory, OldAssembly_Dict, NewAssembly_Dict, OldNewID_Dict):
    """ creates a multifasta file for each assembly
            where seqIDs were not reversed or identical.
            Requires the directories and dictionaries """

    # Old assembly multi-fastafile
    with open("{}/Compare_OldAssembly.fa".format(directory), "w") as output_file:
        for oldID in OldAssembly_Dict.keys():
            if not oldID in OldNewID_Dict.keys():
                output_file.write(">{}\n{}\n".format(oldID, "".join(OldAssembly_Dict[oldID][1])))
        # New assembly multi-fastafile
    output_file.close()
    with open("{}/Compare_NewAssembly.fa".format(directory), "w") as output_file:
        for newID in NewAssembly_Dict.keys():
            if not newID in OldNewID_Dict.values():
                output_file.write(">{}\n{}\n".format(newID, "".join(NewAssembly_Dict[newID][1])))
    output_file.close()


def file_cracker(ToUpdate, file_crack):
    """ Creates entries within the file crack for the subfiles in ToUpdate.
    Requires ToUpdate and file crack """

    for dataset in ToUpdate.keys():
        if(dataset == "Variants"):
            # Loop through the subfiles of the dataset and create a new key for each subfile with an empty list as a value.
            for subfile in ToUpdate[dataset]:
                file_crack["./temporary_directory/updated_{}_A".format(subfile[0])] = []
        # If the dataset is onyl variants (non-genome data) then do this for B too
        elif(dataset != "Genome"):
            for subfile in ToUpdate[dataset]:
                file_crack["./temporary_directory/updated_{}_A".format(subfile[0])] = []
                file_crack["./temporary_directory/updated_{}_B".format(subfile[0])] = []


def inversion_setting(inversions):
    if(inversions == 0):
        inversion_flag = '--forward'
    else:
        inversion_flag = 1
    return inversion_flag


def aligner_caller(aligner_switch, threads, nucmer_directory, reference_directory, query_directory, mashmap_directory,
                   alignment_pickle, directory, delta_directory, percent_identity, kmer, segLength, inversion_flag):
    """
    This function calls the correct aligner, dependent on what the user selects, and executes
    the respective score algorithm for the selected aligner.

    """
    if aligner_switch == 2:
        # Run nucmer with the resulting multifastas. Use threads.
        print("\n\t\t - Running nucmer {}".format(str(datetime.datetime.now())))
        if(inversion_flag == 1):
            Shellcommand_nucmer = ['nucmer', '--mum', '-t', str(threads), reference_directory,
                                   query_directory, '-p', nucmer_directory]
            Popen(Shellcommand_nucmer).wait()
        else:
            Shellcommand_nucmer = ['nucmer', str(inversion_flag), '--mum', '-t', str(threads), reference_directory,
                                   query_directory, '-p', nucmer_directory]
            Popen(Shellcommand_nucmer).wait()
        # Filter the delta file
        print("\n\t\t - Filtering resulting alignment to obtain equivalent "
              "sequences across versions {}".format(str(datetime.datetime.now())))
        filter_nucmer_delta(delta_file="{}.delta".format(nucmer_directory), out_file=delta_directory)
    elif aligner_switch == 1:
        print("\t\t\t  - Running MashMap " + str(datetime.datetime.now()))
        ShellCommand_mashMap = ['/usr/local/bin/mashmap', '-s', str(segLength), '-k', str(kmer), '--pi',
                                str(percent_identity), '-t', str(threads), '-r', query_directory, '-q', reference_directory,
                                '-o', mashmap_directory]
        Popen(ShellCommand_mashMap).wait()
        print("\n\t\t - Filtering and realigning resulting alignment(s) to obtain equivalent "
              "sequences across versions {}".format(str(datetime.datetime.now())))
        # Filter the MashMap output
        filter_mashMap(mash_file=mashmap_directory, directory=directory, alignment_pickle=alignment_pickle,
                       out_file=delta_directory,  threads=str(threads), inversion_flag=inversion_flag)


################
## Main method #
################

def obtain_alignment(old_assembly, new_assembly, directory, threads, ToUpdate, alignment_pickle, aligner_switch,
                     percent_identity, kmer, segLength, inversions):
    """
    ### Function ###
    Performs alignment for the given assemblies.

    ### Requires ###
    Requires the the directories of the assemblies, user specified threads, dictionary of what needs to be
    updated, the alignment pickle, user specified aligner switch, and additional args for MashMap (if used)

    ### returns ###
    Returns a list of the tabix queries, ID dictionary, modified alignment pickle, a summary dictionary, and the
    file crack.

    """
    # Create the pickle directory
    os.makedirs(alignment_pickle)
    # Initiate variables:
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
    inversion_flag = inversion_setting(inversions)

    # summary_Dict {oldID:[status,length]}
    tabix_queries = {}
    OldNewID_Dict = {}
    summary_Dict = {}
    file_crack = {}
    # Define the directories
    query_directory = '{}/Compare_NewAssembly.fa'.format(directory)
    reference_directory = '{}/Compare_OldAssembly.fa'.format(directory)
    mashmap_directory = '{}/genomeGitMash.out'.format(directory)
    nucmer_directory = '{}/ToFilter'.format(directory)
    delta_directory = '{}/Filtered.delta'.format(alignment_pickle)
    snp_directory = '{}/Filtered.snp'.format(alignment_pickle)
    snp_coords_dir = '{}/summary.coords'.format(alignment_pickle)

    for dataset in ToUpdate.keys():
        # If the dataset is variants
        if(dataset == "Variants"):
            # Loop through the subfiles of the dataset and create a new key for each subfile with an empty list as a value.
            for subfile in ToUpdate[dataset]:
                file_crack["./temporary_directory/updated_{}_A".format(subfile[0])] = []
        # Only if the dataset is annotation or variants, do it for B as well
        elif(dataset != "Genome"):
            for subfile in ToUpdate[dataset]:
                file_crack["./temporary_directory/updated_{}_A".format(subfile[0])] = []
                file_crack["./temporary_directory/updated_{}_B".format(subfile[0])] = []

    # Parse the genome assemblies into dictionaries with the following characteristics: {seqKey:[oldID,[line,line...],hash,hash-1]}
    OldAssembly_Dict = parse_assembly(old_assembly)
    NewAssembly_Dict = parse_assembly(new_assembly)

    ##########################################################
    # Create tabix dictionary of sequences for each dataset  #
    ##########################################################

    for oldID in OldAssembly_Dict.keys():
        # Loop through the keys in NewAssembly_Dict
        for newID in NewAssembly_Dict.keys():
            # Check if the SHA1 is identical in the NewAssembly_Dict
            if(OldAssembly_Dict[oldID][2] == NewAssembly_Dict[newID][2]):
                # Add the info to the summary_Dict
                summary_Dict[oldID] = ["identical", str(len("".join(OldAssembly_Dict[oldID][1])))]
                # Add it to the tabix queries dict as an identical sequence. Loop through the files of each dataset.
                for dataset in ToUpdate.keys():
                    # If the dataset is not the genome dataset, create the tabix query
                    if(dataset != "Genome"):
                        # Loop through the files of the dataset  {dataset:[[filename.extension,directory,size],[...]]}
                        for subfile in ToUpdate[dataset]:
                            # Create the query for subfile A
                            query = "tabix ./temporary_directory/{}_A.gz {}:0 > ./temporary_directory/{}_{}:0_A".format(
                                subfile[0], oldID, subfile[0], oldID)
                            tabix_queries[query] = "identical"
                            # Add the information into the file crack. No need to add a discarded file since the sequence is identical
                            file_crack["./temporary_directory/updated_{}_A".format(subfile[0])].append(
                                "./temporary_directory/{}_{}:0_A".format(subfile[0], oldID))
                            # Only Alignment and annotation files need subfile B as well
                            if(dataset == "Alignment" or dataset == "Annotation"):
                                # Repeat with B
                                query = "tabix ./temporary_directory/{}_B.gz {}:0 > ./temporary_directory/{}_{}:0_B".format(
                                    subfile[0], oldID, subfile[0], oldID)
                                tabix_queries[query] = "identical"
                                file_crack["./temporary_directory/updated_{}_B".format(subfile[0])].append(
                                    "./temporary_directory/{}_{}:0_B".format(subfile[0], oldID))
                # It is identical but it their IDs might be different, so add it to the translate dict
                OldNewID_Dict[oldID] = newID
                # No need to keep looking for the equivalent sequence, break
                break
            # Check if the old reverse SHA1 is identical to the new forward SHA1
            elif(NewAssembly_Dict[newID][2] == OldAssembly_Dict[oldID][3]):
                # Add the info to the summary_Dict
                summary_Dict[oldID] = ["reversed", str(len("".join(OldAssembly_Dict[oldID][1])))]
                # Add it to the tabix queries dict as an identical sequence. Loop through the files of each dataset.
                for dataset in ToUpdate:
                    # If the dataset is not the genome dataset, create the tabix query
                    if(dataset != "Genome"):
                        # Loop through the files of the dataset
                        for subfile in ToUpdate[dataset]:
                            # Create the query for subfile A
                            query = "tabix ./temporary_directory/{}_A.gz {}:0 > ./temporary_directory/{}_{}:0_A".format(subfile[0],
                                                                                                                        oldID, subfile[0], oldID)
                            tabix_queries[query] = ["reversed", "_A", dataset, subfile[0], oldID,
                                                    "./temporary_directory/{}_{}:0_A".format(subfile[0], oldID),
                                                    "./temporary_directory/updated_{}_{}:0_A".format(subfile[0], oldID),
                                                    "./temporary_directory/discarded_{}:0_A".format(subfile[0]),
                                                    str(len("".join(OldAssembly_Dict[oldID][1])))]
                            # Add the information into the file crack. No need to add a discarded file since the sequence is identical reversed.
                            file_crack["./temporary_directory/updated_{}_A".format(subfile[0])].append(
                                "./temporary_directory/updated_{}_{}:0_A".format(subfile[0], oldID))
                            # Only Alignment and annotation files need subfile B as well
                            if(dataset == "Alignment" or dataset == "Annotation"):
                                # Do the same with B
                                query = "tabix ./temporary_directory/{}_B.gz {}:0 > ./temporary_directory/{}_{}:0_B".format(subfile[0],
                                                                                                                            oldID, subfile[0], oldID)
                                tabix_queries[query] = ["reversed", "_B", dataset, subfile[0], oldID,
                                                        "./temporary_directory/{}_{}:0_B".format(subfile[0], oldID),
                                                        "./temporary_directory/updated_{}_{}:0_B".format(
                                                            subfile[0], oldID),
                                                        "./temporary_directory/discarded_{}:0_B".format(subfile[0]),
                                                        str(len("".join(OldAssembly_Dict[oldID][1])))]
                                # Add the information into the file crack. No need to add a discarded file since the sequence is identical reversed.
                                file_crack["./temporary_directory/updated_{}_B".format(subfile[0])].append(
                                    "./temporary_directory/updated_{}_{}:0_B".format(subfile[0], oldID))
                # Add Ids to the transalte dict
                OldNewID_Dict[oldID] = newID
                # No need to keep looking for the equivalent sequence, break
                break
    #########################################################
    ## multi FASTA construction from old and new assemblies #
    #########################################################

    # Construct the multifasta files from the old and new assembly
    multiFasta_construct(directory, OldAssembly_Dict, NewAssembly_Dict, OldNewID_Dict)

    ###############################
    ## Alignment and showing snps #
    ###############################

    aligner_caller(aligner_switch, threads, nucmer_directory, reference_directory,
                   query_directory, mashmap_directory, alignment_pickle, directory, delta_directory,
                   percent_identity, kmer, segLength, inversion_flag)

# Run show-snps with the resulting filtered delta file
    print("\n\t\t - Identification of SNPs in sequences {}".format(str(datetime.datetime.now())))
    Popen("show-snps -H -T {} > {}".format(delta_directory, snp_directory), shell=True).wait()
    # Create a dictionary of snps in the filtered snp file {SeqID:"oldIndex oldChar newChar newIndex","...",}
    snp_dict = parse_snp_file("{}/Filtered.snp".format(alignment_pickle))
    # Loop through the lines of the filtered delta file
    print("\n\t\t - Alignment finished. Now creating tabix queries {}".format(str(datetime.datetime.now())))

    ###############################
    ## Create the Tabix queries   #
    ###############################

    with open(delta_directory, "r") as delta_file:
            # Ommit the first line
        delta_file.readline()
        # Initate variables
        prev_snp = 0
        prev_end = "1"
        # Loop through the lines
        for line in delta_file:
                # If the line starts with >, it is a new alignment
            if(line[0] == ">"):
                new_alignment = line.split()
                current_oldID = new_alignment[0][1:]
                # Add the new and old Id to the OldNewID_Dict
                OldNewID_Dict[current_oldID] = new_alignment[1]
                # Re-start the previous end index and prev_snp
                prev_end = "1"
                prev_snp = 0
                # If the seqID has snps, load the new snp list
                if(current_oldID in snp_dict.keys()):
                    snp_list = snp_dict[current_oldID]
                # Otherwise create an empty list
                else:
                    snp_list = []
            # Otherwise add the block (if it has whitespaces it is a block) and its snps into the delta_dict
            elif(" " in line):
                # Split the line
                block = line.rstrip().split()
                # Determine the current snps
                number_snp = int(block[4])
                current_snp = snp_list[int(prev_snp):int(number_snp+prev_snp)]
                prev_snp = number_snp + prev_snp
                # Add it to the tabix queries dict as a to be compared sequence. Loop through the files of each dataset.
                for dataset in ToUpdate:
                    # If the dataset is not the genome dataset, create the tabix query
                    if(dataset != "Genome"):
                            # Loop through the files of the dataset
                        for subfile in ToUpdate[dataset]:
                                # Add the query A
                            query = "tabix ./temporary_directory/{}_A.gz {}:{}-{} > ./temporary_directory/{}_{}:{}_{}_A".format(subfile[0],
                                                                                                                                current_oldID, block[0], block[1], subfile[0], current_oldID, block[0], block[1])
                            tabix_queries[query] = ["compare", "_A", dataset, subfile[0], current_oldID,
                                                    "./temporary_directory/{}_{}:{}_{}_A".format(
                                                        subfile[0], current_oldID, block[0], block[1]),
                                                    "./temporary_directory/updated_{}_{}:{}_{}_A".format(
                                                        subfile[0], current_oldID, block[0], block[1]),
                                                    "./temporary_directory/discarded_{}_{}:{}_{}_A".format(subfile[0], current_oldID, block[0], block[1]), block[0:4], current_snp]
                            # Add the to the file crack the updated and discarded files
                            file_crack["./temporary_directory/updated_{}_A".format(subfile[0])].append(
                                "./temporary_directory/updated_{}_{}:{}_{}_A".format(subfile[0], current_oldID, block[0], block[1]))
                            # Create a query B for alignment and annotation.
                            if(dataset != "Variants"):
                                # Do the same with B
                                query = "tabix ./temporary_directory/{}_B.gz {}:{}-{} > ./temporary_directory/{}_{}:{}_{}_B".format(subfile[0],
                                                                                                                                    current_oldID, block[0], block[1], subfile[0], current_oldID, block[0], block[1])
                                tabix_queries[query] = ["compare", "_B", dataset, subfile[0], current_oldID,
                                                        "./temporary_directory/{}_{}:{}_{}_B".format(
                                                            subfile[0], current_oldID, block[0], block[1]),
                                                        "./temporary_directory/updated_{}_{}:{}_{}_B".format(
                                                            subfile[0], current_oldID, block[0], block[1]),
                                                        "./temporary_directory/discarded_{}_{}:{}_{}_B".format(subfile[0], current_oldID, block[0], block[1]), block[0:4], current_snp]
                                # Add the to the file crack the updated and discarded files
                                file_crack["./temporary_directory/updated_{}_A".format(subfile[0])].append(
                                    "./temporary_directory/updated_{}_{}:{}_{}_A".format(subfile[0], current_oldID, block[0], block[1]))
                # The previous end is now equal to the current region end
                prev_end = str(block[1])

        # If the oldID is not in the OldNew_Dict, it must be absent from the new assembly
        for old_key in OldAssembly_Dict.keys():
            if not(old_key in OldNewID_Dict.keys()):
                # Add the info to the summary_Dict
                summary_Dict[OldAssembly_Dict[old_key][0]] = [
                    "deleted", str(len("".join(OldAssembly_Dict[old_key][1])))]
        # Run show coords
        Popen("show-coords -c -T -H {}/Filtered.delta > {}/summary.coords".format(alignment_pickle,
                                                                                  alignment_pickle), shell=True).wait()

        return [tabix_queries, OldNewID_Dict, alignment_pickle, summary_Dict, file_crack]
