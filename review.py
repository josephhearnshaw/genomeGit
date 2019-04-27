#!/usr/bin/env python

# Arguments: commitA commitB threads merge

# Make the imports
import sys
import os
import shutil
import subprocess
from subprocess import Popen
from StoreAlignment import load_variables, obtain_alignment_pickle, store_variables
from reconstruct_functions import reconstruct_dataset
from parse_functions import parse_delta_file, parse_snp_file, parse_coords_file
from ObtainAlignment_functions import obtain_alignment
from update_functions import detect_updates
from collections import OrderedDict
import gc


# Load arguments: commit_A commit_B threads report
commit_A = str(sys.argv[1])
commit_B = str(sys.argv[2])
number_threads = str(sys.argv[3])
merge = int(sys.argv[4])

# Create a temporary directory
if(os.path.isdir("./.git/info/temporary_directory")):
    shutil.rmtree("./.git/info/temporary_directory")
os.mkdir("./.git/info/temporary_directory")
# Obtain the commit messages
print("\nNow producing a summary review between commits " + commit_A + " and " + commit_B +
      ". Warning: If you have added any data into the repository but did not commit, these changes will be lost.\n")
message_A = subprocess.check_output(
    "git log --format=%B -n 1 " + commit_A, shell=True)
message_B = subprocess.check_output(
    "git log --format=%B -n 1 " + commit_B, shell=True)
# Reconstruct the first version of the assembly
ShellCommand = Popen(
    "git checkout " + str(sys.argv[1]) + " 2> /dev/null", shell=True).wait()
reconstruct_dataset(size=60, directory="./Genome",
                    output_file="./.git/info/temporary_directory/assembly_1.fa", mode="Genome")
# Reconstruct the second version of the assembly
ShellCommand = Popen(
    "git checkout " + str(sys.argv[2]) + " 2> /dev/null", shell=True).wait()
reconstruct_dataset(size=60, directory="./Genome",
                    output_file="./.git/info/temporary_directory/assembly_2.fa", mode="Genome")
# Checkout back to the master
ShellCommand = Popen("git checkout master 2> /dev/null", shell=True).wait()

# Obtain alignment pickle of both assemblies
alignment_pickle1 = obtain_alignment_pickle(
    "./.git/info/temporary_directory/assembly_1.fa", "./.git/info/temporary_directory/assembly_2.fa")
alignment_pickle2 = obtain_alignment_pickle(
    "./.git/info/temporary_directory/assembly_2.fa", "./.git/info/temporary_directory/assembly_1.fa")

# Check If the alignemnt containing the information between these two assemblies already exists
alignment_pickle = ""
if(os.path.isdir(alignment_pickle1)):
    alignment_pickle = alignment_pickle1
elif(os.path.isdir(alignment_pickle2)):
    alignment_pickle = alignment_pickle2

# Parse the information in the alignment files if exists
if(alignment_pickle != ""):
    # Inform the user
    print("A stored alignment has been found. Now loading data.\n")
    # Load the OldNewID dictionary
    queries, alignment_pickle, summary_Dict, \
        oldSeqs, newSeqs = \
        load_variables("{}/pickle".format(alignment_pickle))

# Otherwise it is required to perform the alignment between both assemblies
else:
    # Inform the user
    print("No stored alignment was found. Now creating genome alignment. This will take some time.\n")
    # Detect the updates
    ToUpdate = detect_updates("./RepoMap.txt")
    # Obtain the variables
    variables = obtain_alignment(old_assembly="./.git/info/temporary_directory/assembly_1.fa",
                                 new_assembly="./.git/info/temporary_directory/assembly_2.fa",
                                 directory="./.git/info/temporary_directory",
                                 threads=number_threads, ToUpdate=ToUpdate,
                                 alignment_pickle=alignment_pickle, aligner_switch=2,
                                 percent_identity=95, kmer=15, segLength=5000, c_flag=2000, b_flag=1, ms_flag=1)

    # Store the pickle [tabix_queries,OldNewID_Dict,alignment_pickle,summary_Dict]
    store_variables(variables=variables, alignment_pickle=alignment_pickle)
    # Load the variables
    queries, alignment_pickle, summary_Dict, oldSeqs, newSeqs = variables

# Parse the coords file {oldID:[contigs]...}
coords_dict = parse_coords_file(
    "{}/summary.coords".format(alignment_pickle))

# Parse the delta file {oldID:[length_A,length_B,score]...}
delta_list = parse_delta_file("{}/Filtered.delta".format(alignment_pickle))


# initiate empty containers
identical_list = []
reversed_list = []
del_list = []
merge_dict = OrderedDict()
split_dict = OrderedDict()

# iterate through the summary_Dict to get identical, reversed and deleted sequences
for oldID, (tag, oldLength, newID) in summary_Dict.items():
    if tag == 'identical':
        identical_list.append("\t------\n\t|\n\t| Sequence: {} ({} nt) --> {} (identical)\n\t|\n\t------".format(
            oldID, oldLength, newID))
    elif tag == 'reversed':
        reversed_list.append("\t------\n\t|\n\t| Sequence: {} ({} nt) --> {} (reversed)\n\t|\n\t------".format(
            oldID, oldLength, newID))
    elif tag == 'deleted':
        del_list.append("\t------\n\t|\n\t| Sequence: {} ({} nt) --> deleted\n\t|\n\t------".format(
            oldID, oldLength, newID))

# iterate through the delta_list and populate the split and merge dicts
for oldID, newID, oldLength, newLength in delta_list:
    if newID in merge_dict.keys():
        merge_dict[newID].append(oldID)
    else:
        merge_dict[newID] = [oldID]
    if oldID in split_dict.keys():
        split_dict[oldID].append(newID)
    else:
        split_dict[oldID] = [newID]

# open the output file and write the header
with open("../GenomeGit_report.txt", "w") as report_file:
    report_file.write('DIFF REPORT FOR:')
    report_file.write("\n\t###\n\t#\tVERSION A: {} ({})\n\t###\n".format(
        message_A.rstrip(), commit_A))
    report_file.write("\n\t###\n\t#\tVERSION B: {} ({})\n\t###\n\n".format(
        message_B.rstrip(), commit_B))

    # write the identical seqs
    report_file.write('\n\n\tIdentical sequences:\n')
    for entry in identical_list:
        report_file.write(entry)

    # write the reversed seqs
    report_file.write('\n\n\tReversed sequences:\n')
    for entry in reversed_list:
        report_file.write(entry)

    # write the deleted seqs
    report_file.write('\n\n\tDeleted sequences:\n')
    for entry in del_list:
        report_file.write(entry)

    # write the merged seqs
    report_file.write('\n\n\tMerges:\n')
    for newID, oldIDs in merge_dict.items():
        if len(oldIDs) > 1:
            report_file.write('\t------\n\t|\n\t| {}\t-->\t{}\n\t|\n\t------'.format(
                ', '.join(oldIDs), newID))

    report_file.write('\n\n\tSplits:\n')
    for oldID, newIDs in split_dict.items():
        if len(newIDs) > 1:
            report_file.write('\t------\n\t|\n\t| {}\t-->\t{}\n\t|\n\t------'.format(
                oldID, ', '.join(newIDs)))

    # call the garbage collector to free up the memory used by the lists and dicts
    gc.collect()

    # iterate through the alignments, perform analysis (count SNPs and indels,
    # get the block coverage) and write to the output file
    report_file.write('\n\n\tAlignments:\n')
    for oldID, newID, oldLength, newLength in delta_list:
        report_file.write(
            "\n\n\t------\n\t|\n\t| Sequence: {} ({} nt) --> {} ({} nt)\n\t|\n".format(
                oldID, oldLength, newID, newLength))
        # Add information regarding contigs. Initiate variables
        i = 0
        coverage_quer = 0
        coverage_ref = 0
        # retrieve the blocks of the current alignment
        for (ref_start, ref_end, query_start, query_end, ref_length, query_length,
                similarity, ref_coverage, query_coverage) in coords_dict[(oldID, newID)]:
            i += 1
            report_file.write("\t|\t-Block {}: {} - {} --> {} - {}\t ({}% sim.)\n".format(
                str(i), ref_start, ref_end, query_start, query_end, similarity))
            # Sum coverage
            coverage_ref = coverage_ref + float(ref_coverage)
            coverage_quer = coverage_quer + float(query_coverage)
        # Append the closing of the table
        report_file.write(
            "\t|\n\t| Total coverage:\n\t|\n\t|\t-Version A: {}%\n\t|\t-Version B: {}%\n".format(
                coverage_ref, coverage_quer))
        report_file.write("\t|\n \t------\n")
    report_file.write(' ------\n\nREPORT END')
