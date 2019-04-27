#!/usr/bin/env python

# PYTHON FUNCTIONS REQUIRED IN THE REPOSITORY AUTO-UPDATES PART II: INTERPRETATION OF ALIGNEMNT

# Make the imports
import os
import sys
import datetime
import shutil
from subprocess import Popen
from parse_functions import parse_dataset
import multiprocessing
from ObtainAlignment_functions import reverse_complement
import tabix

# Create a lock for update_sequence
lock = multiprocessing.Lock()


# Create a function to eliminate repeated barcodes in sub SAM files (this can mess arround with the merge process)
def prepare_barcode(infile, outfile):
    # Initiate the index and open the input file
    old_index = "0"
    with open(infile, "r") as input_file:
        # Open the output file and loop through the lines of the input file
        with open(outfile, 'w') as output_file:
            for line in input_file:
                line = line.split("\t")
                new_index = line[-1]
                # If the new index is equal to the old one, skip the line
                if (new_index == old_index):
                    continue
                # Otherwise include the line in the output file
                else:
                    output_file.write("\t".join(line))
                old_index = new_index


# Create a function to merge the two SAM/GFF/VCF subfiles
def merge_subfiles(dataset, subfile_name, template_length):
    # Inform the user
    print("\t\t - Merging metadata subfiles for: {} at {}".format(subfile_name,
                                                                  str(datetime.datetime.now())))
    sys.stdout.flush()
    # Open the metadata
    metadata = open(
        "./temporary_directory/{}_metadata".format(subfile_name), "r")
    # Read through all the comments of the metadata
    line_metadata = metadata.readline()
    while (line_metadata[0] == "@" or line_metadata[0] == "#"):
        line_metadata = metadata.readline()
    # If the dataset is alignment, there are two subfiles
    if (dataset == "Alignment"):
        # Open the updated and discarded files
        # Open the updated file in append mode since it already contains the comments
        updated_file = open(
            "./temporary_directory/{}".format(subfile_name), "a")
        discarded_file = open(
            "./temporary_directory/{}.discarded".format(subfile_name), "w")
        # First thing required is to prepare the barcodes in the subfiles: eliminate repeated barcodes
        prepare_barcode(infile="./temporary_directory/sorted_updated_{}_A".format(subfile_name),
                        outfile="./temporary_directory/BarcodeReady_{}_A".format(subfile_name))
        prepare_barcode(infile="./temporary_directory/sorted_updated_{}_B".format(subfile_name),
                        outfile="./temporary_directory/BarcodeReady_{}_B".format(subfile_name))
        # Open the files
        file_A = open(
            "./temporary_directory/BarcodeReady_{}_A".format(subfile_name), "r")
        file_B = open(
            "./temporary_directory/BarcodeReady_{}_B".format(subfile_name), "r")
        # Read the first line of metadata and the files A & B. Split them by the tabs
        line_metadata = line_metadata.split("\t")
        line_A = file_A.readline().rstrip().split("\t")
        line_B = file_B.readline().rstrip().split("\t")
        # Start looping through the metadata file
        while (len(line_metadata) > 1):
            # If the one of the reads is already empty, discard the rest of the metadata
            if (len(line_A) != 3 or len(line_B) != 3):
                discarded_file.write("\t".join(line_metadata[1:]))
                line_metadata = metadata.readline().split("\t")
            # Otherwise updata the metadata
            else:
                # Obtain the barcodes
                barcode_metadata = line_metadata[0]
                barcode_A = line_A[2]
                barcode_B = line_B[2]
                # Check if the barcodes match. If they do, write the new line with the new TLENGHT
                if (barcode_metadata == barcode_A and barcode_metadata == barcode_B):
                    # Read A new coordiantes and newID
                    line_metadata[3] = line_A[0]
                    line_metadata[4] = line_A[1]
                    # New coordinate for read B
                    line_metadata[8] = line_B[1]
                    # New seqID for read B: if the seqID is identical in both reads, place a = and determine the TLENGTH
                    if (line_B[0] == line_A[0]):
                        line_metadata[7] = "="
                        line_metadata[9] = str(int(line_B[1]) - int(line_A[1]))
                        # Write the line in the output file,omit the barcode.
                        updated_file.write("\t".join(line_metadata[1:]))
                        # Discard if it is higher than the threshold or lower than zero
                        # if(int(line_metadata[9])>template_length or int(line_metadata[9])<0):
                        # CURRENTLY ONLY NOT CHECKING TLENGTH VALUE
                        # Tlength is too high/low, discard the entry
                        # discarded_file.write("\t".join(line_metadata[1:]))
                        # Otherwise it is good to go into the updated file
                        # else:
                        # Write the line in the output file,omit the barcode.
                        # updated_file.write("\t".join(line_metadata[1:]))
                    # Otherwise they map different reads, place a 0 as TLENGTH
                    else:
                        line_metadata[7] = line_B[0]
                        line_metadata[9] = "0"
                        # Write the line in the output file,omit the barcode.
                        updated_file.write("\t".join(line_metadata[1:]))
                    # Read the next line of the three files and spplit them using the tabs
                    line_metadata = metadata.readline().split("\t")
                    line_A = file_A.readline().rstrip().split("\t")
                    line_B = file_B.readline().rstrip().split("\t")
                # If the barcodes dont match, one of the reads was discarded at some point.
                # If the barcode of both reads match, then both reads were discarded. Discard the
                # entire entry and loop to the next line of the metadata
                elif (barcode_A == barcode_B):
                    # Discard the entry
                    discarded_file.write("\t".join(line_metadata[1:]))
                    # Read the next line
                    line_metadata = metadata.readline().split("\t")
                # Otherwise only one of the reads was discarded.
                # Need to discard the entire line and re-adjust both reads
                else:
                    # Discard the entry
                    discarded_file.write("\t".join(line_metadata[1:]))
                    # Determine which read was discarded
                    if (int(barcode_A) > int(barcode_B)):
                        # If line A was discarded, read the next line in the file B and split it by the tabs
                        line_B = file_B.readline().rstrip().split("\t")
                        # Read as well the next one in the metadata
                        line_metadata = metadata.readline().split("\t")
                    else:
                        # If line B was discarded, read the next line in the file A and split it by the tabs
                        line_A = file_A.readline().rstrip().split("\t")
                        # Read as well the next one in the metadata
                        line_metadata = metadata.readline().split("\t")
        # Close the files
        file_A.close()
        file_B.close()
        updated_file.close()
        discarded_file.close()
    # If it is annotation
    elif (dataset == "Annotation"):
        # Open the updated and discarded files
        # Open the updated file in append mode since it already contains the comments
        updated_file = open(
            "./temporary_directory/{}".format(subfile_name), "a")
        discarded_file = open(
            "./temporary_directory/{}.discarded".format(subfile_name), "w")
        # First thing required is to prepare the barcodes in the subfiles: eliminate repeated barcodes
        prepare_barcode(infile="./temporary_directory/sorted_updated_{}_A".format(subfile_name),
                        outfile="./temporary_directory/BarcodeReady_{}_A".format(subfile_name))
        prepare_barcode(infile="./temporary_directory/sorted_updated_{}_B".format(subfile_name),
                        outfile="./temporary_directory/BarcodeReady_{}_B".format(subfile_name))
        # Open the files
        file_A = open(
            "./temporary_directory/BarcodeReady_{}_A".format(subfile_name), "r")
        file_B = open(
            "./temporary_directory/BarcodeReady_{}_B".format(subfile_name), "r")
        # Read the first line of metadata and the files A & B. Split them by the tabs
        line_metadata = line_metadata.split("\t")
        line_A = file_A.readline().rstrip().split("\t")
        line_B = file_B.readline().rstrip().split("\t")
        # Start looping through the metadata file
        while (len(line_metadata) > 1):
            # If the one of the reads is already empty, discard the rest of the metadata
            if (len(line_A) != 3 or len(line_B) != 3):
                discarded_file.write("\t".join(line_metadata[1:]))
                line_metadata = metadata.readline().split("\t")
            # Otherwise updata the metadata
            else:
                # Obtain the barcodes
                barcode_metadata = line_metadata[0]
                barcode_A = line_A[2]
                barcode_B = line_B[2]
                # Check if the barcodes match. If they do, write the new line with the new coordinates and seqIDs
                if (barcode_metadata == barcode_A and barcode_metadata == barcode_B):
                    # Only if both parts have been mapped to the same region, add them to the updated file
                    if (line_A[0] == line_B[0]):
                        # Part A new coordiantes and newID
                        line_metadata[1] = line_A[0]
                        line_metadata[4] = line_A[1]
                        # New coordinate for part B
                        line_metadata[5] = line_B[1]
                        # Write the entry
                        updated_file.write("\t".join(line_metadata[1:]))
                    # Otherwise discard the read
                    else:
                        discarded_file.write("\t".join(line_metadata[1:]))
                    # Read the next line of the three files and spplit them using the tabs
                    line_metadata = metadata.readline().split("\t")
                    line_A = file_A.readline().rstrip().split("\t")
                    line_B = file_B.readline().rstrip().split("\t")
                # If the barcodes dont match, one of the reads was discarded at some point.
                # If the barcode of both reads match, then both reads were discarded. Discard the
                # entire entry and loop to the next line of the metadata
                elif (barcode_A == barcode_B):
                    # Discard the entry
                    discarded_file.write("\t".join(line_metadata[1:]))
                    # Read the next line
                    line_metadata = metadata.readline().split("\t")
                # Otherwise only one of the reads was discarded.
                # Need to discard the entire line and re-adjust both reads
                else:
                    # Discard the entry
                    discarded_file.write("\t".join(line_metadata[1:]))
                    # Determine which read was discarded
                    if (int(barcode_A) > int(barcode_B)):
                        # If line A was discarded, read the next line in the file B and split it by the tabs
                        line_B = file_B.readline().rstrip().split("\t")
                        # Read as well the next one in the metadata
                        line_metadata = metadata.readline().split("\t")
                    else:
                        # If line B was discarded, read the next line in the file A and split it by the tabs
                        line_A = file_A.readline().rstrip().split("\t")
                        # Read as well the next one in the metadata
                        line_metadata = metadata.readline().split("\t")
        # Close the files
        file_A.close()
        file_B.close()
        updated_file.close()
        discarded_file.close()
    # Otherwise it is variants: only one subfile
    else:
        # Open the updated and discarded files
        # Open the updated file in append mode since it already contains the comments
        updated_file = open(
            "./temporary_directory/{}".format(subfile_name), "a")
        discarded_file = open(
            "./temporary_directory/{}.discarded".format(subfile_name), "w")
        # First thing required is to prepare the barcodes in the subfiles: eliminate repeated barcodes
        prepare_barcode(infile="./temporary_directory/sorted_updated_{}_A".format(subfile_name),
                        outfile="./temporary_directory/BarcodeReady_{}_A".format(
                            subfile_name))
        # Open the files
        file_A = open(
            "./temporary_directory/BarcodeReady_{}_A".format(subfile_name), "r")
        # Read the first line of metadata and A. Split them by the tabs
        line_metadata = line_metadata.split("\t")
        line_A = file_A.readline().rstrip().split("\t")
        # Start looping through the metadata file
        while (len(line_metadata) > 1):
            # If the one of the reads is already empty, discard the rest of the metadata
            if (len(line_A) != 5):
                discarded_file.write("\t".join(line_metadata[1:]))
                line_metadata = metadata.readline().split("\t")
            # Otherwise update the metadata
            else:
                # print(line_metadata)
                # print(line_A)
                # Obtain the barcodes
                barcode_metadata = line_metadata[0]
                barcode_A = line_A[4]
                # Check if the barcodes match. If they do, write the new line with the new TLENGHT
                if (barcode_metadata == barcode_A):
                    # Read A new coordiantes, newID and updated bases
                    line_metadata[1] = line_A[0]
                    line_metadata[2] = line_A[1]
                    line_metadata[4] = line_A[2]
                    line_metadata[5] = line_A[3]
                    # Write the line in the output file,omit the barcode.
                    updated_file.write("\t".join(line_metadata[1:]))
                    # Read the next line
                    line_metadata = metadata.readline().split("\t")
                    line_A = file_A.readline().rstrip().split("\t")
                # Otherwise the variant was discarded. Need to discard the entire line and re-adjust
                else:
                    # Discard the entry
                    discarded_file.write("\t".join(line_metadata[1:]))
                    # Re-adjust by reading the next line in the metadata and split it by the tabs
                    line_metadata = metadata.readline().split("\t")
        # Close the files
        file_A.close()
        updated_file.close()
        discarded_file.close()

    # Close metadata
    metadata.close()



def update_sequence(query_obj, query_count, number_of_queries):
    """
    Create a function to analyse the alignment of a given sequence and append the updated/discarded
    entries into the global output files
    only change seqIDs if sequence is identical
    """

    # Change seqIDs
    if query_obj.type == "identical":
        process_identical_query(query_obj)

    # Change seqIDs and coordinates (plus bases for variants), if the sequence has been reversed
    elif query_obj.type == "reversed":
        process_reversed_query(query_obj)

    # Otherwise sequence alignment must be evaluated.
    elif query_obj.type == 'compare':
        process_compare_query(query_obj)

    if (query_count + 1) % 1000 == 0:
        print('\t - Approximately. {} out of {} tabix queries processed'.format(
            query_count + 1, number_of_queries))


def process_identical_query(query_obj):
    """
    Loops through the output of a 'identical' query and processes all entries. The function used to update
    the entries is defined before the loop to reduce the if-clauses evaluated in every iteration (slight
    performance increase)
    """

    # Execute the tabix query, skip query if empty result
    try:
        tb = tabix.open(query_obj.originalFile)
        records = tb.query(query_obj.oldSeqID, 0, query_obj.oldSeqLength)
    except tabix.TabixError:
        return

    # get frequently used attributes for faster lookup
    newID = query_obj.newSeqID

    # open files
    # lock is acquired then released after writing is done
    try:
        lock.acquire()
        with open(query_obj.dependentFile, 'a') as updated_file:
            # Loop through tabix output
            for entry in records:
                # Modify the oldID with the newID
                entry[0] = newID
                # check if the updated coordinate is negative. if yes, the entry is discarded
                if int(entry[1]) < 0:
                    continue
                updated_file.write("\t".join(entry))
                updated_file.write("\n")
    finally:
        lock.release()


def process_reversed_query(query_obj):
    """
    Loops through the output of a 'reversed' query and processes all entries. The function used to update
    the entries is defined before the loop to reduce the if-clauses evaluated in every iteration (slight
    performance increase)
    """

    # Execute the tabix query, skip query if empty result
    try:
        tb = tabix.open(query_obj.originalFile)
        records = tb.query(query_obj.oldSeqID, 0, query_obj.oldSeqLength)
    except tabix.TabixError:
        return
    # get frequently used attributes for faster lookup
    newID = query_obj.newSeqID
    length = query_obj.newSeqLength
    # determine which update function to use
    update_func = find_update_entry_function(
        inversed=True, dataset=query_obj.dataset)

    # initialize entry list
    entryList = []

    # Loop through tabix output
    for entry in records:
        # update and write the entry
        entry = update_func(entry, newID, length + 1)
        # check if the updated coordinate is negative. if yes, the entry is discarded
        if int(entry[1]) < 0:
            # print('negative:\t{}'.format(entry))
            continue
        entryList.append('\t'.join(entry))

    # lock is acquired then released after writing is done
    try:
        lock.acquire()
        # open the files
        with open(query_obj.dependentFile, 'a') as updated_file:
            for entry in entryList:
                updated_file.write(entry)
                updated_file.write("\n")
    finally:
        lock.release()


def process_compare_query(query_obj):
    """
    Loops through the output of a 'compare' query and processes all entries. The function used to update
    the entries is defined before the loop to reduce the if-clauses evaluated in every iteration (slight
    performance increase).
    Uses the SNPs list (as generator) to loop through entries and SNPs concomitantly.
    """

    # Execute the tabix query, skip query if empty result
    try:
        tb = tabix.open(query_obj.originalFile)
        start = int(query_obj.block[0])
        end = int(query_obj.block[1])
        records = tb.query(query_obj.oldSeqID, start, end)
    except tabix.TabixError:
        return

    # put some attributes into variables for faster lookup (slightly better performance on many calls)
    newID = query_obj.newSeqID
    ref_start, ref_end, qry_start, qry_end = (int(x) for x in query_obj.block)
    # determine if the aligned black is inversed
    inversed = True if qry_start > qry_end else False
    # the displacement_factor needed for updating coordinates downstream is calculated differently
    # depending on whether the block is inversed or not
    displacement_factor = qry_start + ref_start if inversed else qry_start - ref_start
    # determine which update function to use
    update_func = find_update_entry_function(inversed, query_obj.dataset)

    entryList = []
    # if there are SNPs, process them
    if (query_obj.SNPs):
        # make SNPs list iterator from the list so that it can be iterated over with next()
        SNPs = iter(query_obj.SNPs)
        curr_snp = next(SNPs)
        snp_coord = int(curr_snp[0])
        # loop through the entries, process SNPs and update coordinates and bases (if inversed)
        for entry in records:
            entry_coord = int(entry[1])
            while (curr_snp != 'done' and snp_coord <= entry_coord):
                entry, displacement_factor = process_snp(entry, curr_snp, entry_coord, snp_coord,
                                                         query_obj.dataset, displacement_factor,
                                                         inversed)
                try:
                    curr_snp = next(SNPs)
                    snp_coord = int(curr_snp[0])
                except StopIteration:
                    # the iterator is exhausted --> set curr_snp to done so that the
                    # while loop won't be entered any more.
                    curr_snp = 'done'
            if entry[-1] == 'discarded':
                continue
            entry = update_func(entry, newID, displacement_factor)
            # check if updated coordinate is negative and discard entry in this case
            if int(entry[1]) < 0:
                continue
            # join updated entry and print to file
            entryList.append('\t'.join(entry))

    # iterate over the entries returned by the query
    else:
        for entry in records:
            entry = update_func(entry, newID, displacement_factor)
            # join updated entry and print to file
            entryList.append('\t'.join(entry))

    # lock is acquired then released after writing is done
    try:
        lock.acquire()
        # open the files
        with open(query_obj.dependentFile, 'a') as updated_file:
            for entry in entryList:
                updated_file.write(entry)
                updated_file.write("\n")
    finally:
        lock.release()


def find_update_entry_function(inversed, dataset):
    if dataset == 'Variants':
        if inversed:
            return update_inversed_variance_entry
        else:
            return update_variance_entry
    else:
        if inversed:
            return update_inversed_non_variance_entry
        else:
            return update_non_variance_entry


def update_inversed_non_variance_entry(entry, newID, displacement_factor):
    new_entry_coord = str(displacement_factor - int(entry[1]))
    entry[0], entry[1] = newID, new_entry_coord
    return entry


def update_non_variance_entry(entry, newID, displacement_factor):
    new_entry_coord = str(displacement_factor + int(entry[1]))
    entry[0], entry[1] = newID, new_entry_coord
    return entry


def update_inversed_variance_entry(entry, newID, displacement_factor):
    base1 = reverse_complement(entry[2])
    base2 = reverse_complement(entry[3])
    new_entry_coord = str(displacement_factor - int(entry[1]))
    entry = (newID, new_entry_coord, base1, base2, entry[4])
    return entry


def update_variance_entry(entry, newID, displacement_factor):
    base1 = entry[2]
    base2 = entry[3]
    new_entry_coord = str(displacement_factor + int(entry[1]))
    entry = (newID, new_entry_coord, base1, base2, entry[4])
    return entry


def process_snp(entry, curr_snp, entry_coord, snp_coord, dataset, displacement_factor, inversed=False):
    displacement_change = 0
    if (entry_coord == snp_coord and dataset == "Variants"):
        # if the variant has been deleted, replace the barcode with 'discarded' so that it is discarded downstream
        if curr_snp[2] == '.':
            entry[-1] = 'discarded'
        entry[2] = curr_snp[2]
    if (curr_snp[1] == "."):
        # Add one to the entry index
        displacement_change = 1
    # It is a deletion
    elif (curr_snp[2] == "."):
        # Remove one to the entry indexes
        displacement_change = -1
    if inversed:
        displacement_change *= -1
    return entry, displacement_factor + displacement_change


def interpret_alignment(queries, threads, ToUpdate, tlength, new_assembly):
    """
    Add the comments of the original files into the beginning of the updated ones.
    Loop through the datasets and act accordingly.
    {dataset:[[filename.extension,directory,size],[]...]}
    """

    number_of_queries = len(queries)

    for dataset in ToUpdate.keys():
        # No comments in the genome dataset
        if (dataset != "Genome"):
            # Loop through the files of the dataset
            for subfile in ToUpdate[dataset]:
                # Add the comments of the original file into the updated one
                # Ignore alignment, new alignment headers will be generated later
                if (dataset != "Alignment"):
                    ShellCommand = Popen(
                        "cp ./" + subfile[1] + "/Comments.txt ./temporary_directory/" + subfile[0], shell=True).wait()
    # From now on this part can be threaded, one thread analysing one tabix query at a time.
    # The dictionaries must be made global otherwise they will be overwriten.

    # While there are sequence alignments to be processed
    print("\n\t\t - Now processing {} tabix queries at {}".format(number_of_queries,
                                                                  str(datetime.datetime.now())))

    # Create _A and _B updated files before processing tabix queries
    # This is to avoid multiple processes creating files at the same time
    updatedFileList = []
    for dataset in ToUpdate.keys():
        # If the dataset is variants
        if (dataset == "Variants"):
            # Loop through the subfiles of the dataset and create a new key
            # for each subfile with an empty list as a value.
            for subfile in ToUpdate[dataset]:
                updatedFileList.append(
                    open("./temporary_directory/updated_{}_A".format(subfile[0]), 'w'))
        # Only if the dataset is annotation or variants, do it for B as well
        elif (dataset != "Genome"):
            for subfile in ToUpdate[dataset]:
                updatedFileList.append(
                    open("./temporary_directory/updated_{}_A".format(subfile[0]), 'w'))
                updatedFileList.append(
                    open("./temporary_directory/updated_{}_B".format(subfile[0]), 'w'))

    # Process tabix queries with multi-processing
    pool = multiprocessing.Pool(threads)
    for i, query in enumerate(queries):
        pool.apply_async(update_sequence, args=(query, i, number_of_queries))
    pool.close()
    pool.join()

    # Close all _A and _B files
    for file in updatedFileList:
        file.close()

    # Finally, when the files are created it is required to sort them. Loop through the datasets.
    for dataset in ToUpdate.keys():
        if (dataset != "Genome"):
            # Loop through the files of the dataset
            for subfile in ToUpdate[dataset]:
                # Inform the user
                print("\n\t\t - Now parsing updated file {} into the repository structure at {}".format(
                    subfile[0], str(datetime.datetime.now())))
                sys.stdout.flush()
                # If the dataset is alignment
                if (dataset == "Alignment"):
                    # Append the missing unmapped reads into the updated files (these reads wont be included otherwise)
                    ShellCommand = Popen("tabix ./temporary_directory/{}_A.gz *:0 >> "
                                         "./temporary_directory/updated_{}_A"
                                         .format(subfile[0], subfile[0]), shell=True).wait()
                    ShellCommand = Popen("tabix ./temporary_directory/{}_B.gz *:0 >> "
                                         "./temporary_directory/updated_{}_B"
                                         .format(subfile[0], subfile[0]), shell=True).wait()
                    # Sort the reads in both A and B files using the barcode
                    ShellCommand = Popen("sort --numeric-sort -k 3 ./temporary_directory/updated_{}_A > "
                                         "./temporary_directory/sorted_updated_{}_A"
                                         .format(subfile[0], subfile[0]), shell=True).wait()
                    ShellCommand = Popen("sort --numeric-sort -k 3 ./temporary_directory/updated_{}_B > "
                                         "./temporary_directory/sorted_updated_{}_B"
                                         .format(subfile[0], subfile[0]), shell=True).wait()
                elif (dataset == "Annotation"):
                    # Sort the reads in both A and B files using the barcode
                    ShellCommand = Popen("sort --numeric-sort -k 3 ./temporary_directory/updated_{}_A > "
                                         "./temporary_directory/sorted_updated_{}_A"
                                         .format(subfile[0], subfile[0]), shell=True).wait()
                    ShellCommand = Popen("sort --numeric-sort -k 3 ./temporary_directory/updated_{}_B > "
                                         "./temporary_directory/sorted_updated_{}_B"
                                         .format(subfile[0], subfile[0]), shell=True).wait()
                # Otherwise it is variants
                else:
                    # Sort the reads by barcode
                    ShellCommand = Popen("sort --numeric-sort -k 5 ./temporary_directory/updated_{}_A > "
                                         "./temporary_directory/sorted_updated_{}_A"
                                         .format(subfile[0], subfile[0]), shell=True).wait()
                # Merge all the files
                merge_subfiles(
                    dataset=dataset, subfile_name=subfile[0], template_length=int(tlength))

                # Add new headers to the new alignment file
                if (dataset == "Alignment"):
                    ShellCommand = Popen("samtools faidx " + new_assembly, shell=True).wait()
                    ShellCommand = Popen("samtools view -ht " + new_assembly + ".fai " +
                                         "./temporary_directory/{}".format(subfile[0]) +
                                         " > ./temporary_directory/{}".format(subfile[0] + "_commented"),
                                         shell=True).wait()
                    os.rename("./temporary_directory/{}".format(subfile[0] + "_commented"),
                              "./temporary_directory/{}".format(subfile[0]))
                    os.remove(new_assembly + ".fai")

                # make the directory
                git_directory = "./{}/{}".format(dataset, subfile[0])
                # Parse the updated file into the repository (not necessary to git add, that is
                # done in the genomegit main wrapper). First, delete the old directory.
                shutil.rmtree(git_directory)
                # Create a new directory a move inside
                os.mkdir(git_directory)
                os.chdir(git_directory)
                # Parse the file
                parse_dataset(dataset=dataset, input_path="../../temporary_directory/{}".format(subfile[0]),
                              size=str(os.path.getsize("../../temporary_directory/" + subfile[0]))[:-1], update="1")
                # Go back to the base directory
                os.chdir("../../")
                # Add the discarded entries into the file directory

                shutil.copyfile("./temporary_directory/{}.discarded".format(
                    subfile[0]), "./{}/{}/Discarded".format(dataset, subfile[0]))
