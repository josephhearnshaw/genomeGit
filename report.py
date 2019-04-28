#!/usr/bin/env python

# Make the imports
import sys
import os
import subprocess
from subprocess import Popen
from report_functions import report_repo

# Load the arguments. Obtain commit hash.
commit = str(sys.argv[1])
# If there is no genome dataset, abort
if not (os.path.isdir("./Genome")):
    print("*** THE REPOSITORY DOES NOT CONTAIN A GENOME DATASET AT THE SELECTED COMMIT. "
          "A GENOME DATASET IS REQUIRED FOR THE CREATION OF A REPORT. ***\nNow aborting.")
    sys.exit()

# Create an output file
output_file = open("../GenomeGit_Report.txt", "w")
# Check if the user has selected a commit. If so, load the data.
if (commit != "0"):
    # Checkout the desired version of the assembly
    message = subprocess.check_output("git checkout " + commit, shell=True)
    # Inform the user
    print(
        "\nNow producing a report of the data contained in the repository for the commit: " + message +
        ". Warning: If you have added any data into the repository but did not commit, these changes will be lost.\n")
    output_file.write(
        "\n\nNow producing a report of the data contained in the repository for the commit: " + message +
        ". Warning: If you have added any data into the repository but did not commit, these changes will be lost.\n")
else:
    # Inform the user
    print("\nNow producing a report of the data currently contained in the repository\n")
    output_file.write("\n\nNow producing a report of the data currently contained in the repository\n")

# Create dictionary about the genomic sequences
# {seqKey:[seqID,lenght,contig_count,N_count,[{file:vcf},{file:gff},{file:sam}]]}
genome_report = report_repo()

# Print the information of the dictionary in form of a table
# Loop through the sequences
for seqKey in genome_report.keys():
    # Print the seqID
    output_file.write(
        "\n\t------\n\t|\n\t| Sequence: " + genome_report[seqKey][0] + " (" + genome_report[seqKey][1] + " nt)")
    # Print the number of contigs, gaps and the percentage of Ns
    output_file.write("\n\t| Number of contigs: " + str(genome_report[seqKey][2]) + "\n\t| Number of gaps: " + str(
        genome_report[seqKey][2] - 1) + "\n\t| Percentage of unknown sequence: " + str(
        genome_report[seqKey][3] / (int(genome_report[seqKey][1]) / 100)) + "%")
    # Print VCF features. If there are no VCF files
    if (len(genome_report[seqKey][4][0].keys()) == 0):
        output_file.write("\n\t|\n\t| VCF features: There are no VCF features mapped into this sequence.\n\t|")
    # If there are VCF files
    else:
        output_file.write("\n\t|\n\t| VCF features:\n\t|")
        for vcf_file in genome_report[seqKey][4][0].keys():
            output_file.write("\n\t|\t-" + vcf_file + ": " + str(genome_report[seqKey][4][0][vcf_file]))
    # Print GFF features. If there are no GFF files
    if (len(genome_report[seqKey][4][1].keys()) == 0):
        output_file.write("\n\t|\n\t| GFF features: There are no GFF features mapped into this sequence.\n\t|")
    # If there are GFF files
    else:
        output_file.write("\n\t|\n\t| GFF features:\n\t|")
        for gff_file in genome_report[seqKey][4][1].keys():
            output_file.write("\n\t|\t-" + gff_file + ": " + str(genome_report[seqKey][4][1][gff_file]))
    # Print SAM features. If there are no SAM files
    if (len(genome_report[seqKey][4][2].keys()) == 0):
        output_file.write("\n\t|\n\t| SAM features: There are no SAM features mapped into this sequence.\n\t|")
    # If there are SAM files
    else:
        output_file.write("\n\t|\n\t| SAM features:\n\t|")
        for sam_file in genome_report[seqKey][4][2].keys():
            output_file.write("\n\t|\t-" + sam_file + ": " + str(genome_report[seqKey][4][2][sam_file]))
    # Print the end of the sequence
    output_file.write("\n\t|\n\t------")
# If the user wanted to report a different version of the repo, go back to the current version
if (commit != "0"):
    # Checkout back to the master
    ShellCommand = Popen("git checkout master", shell=True).wait()
# Inform the user
output_file.write("\n\n*End of the report.\n")
# Close the output file
output_file.close()

# Print in the prompt some basic statistics
with open("./RepoMap.txt", "r") as repomap:
    variants_number = 0
    alignment_number = 0
    annotation_number = 0
    for line in repomap:
        line = line.split("\t")
        if (line[1] == "Annotation"):
            annotation_number += 1
        elif (line[1] == "Alignment"):
            alignment_number += 1
        elif (line[1] == "Variants"):
            variants_number += 1
repomap.close()

print("Number of sequences stored in the Genome dataset: " + str(len(genome_report.keys())))
print("Number of files contained in the Annotation dataset: " + str(annotation_number))
print("Number of files contained in the Variants dataset: " + str(variants_number))
print("Number of files contained in the Alignment dataset: " + str(alignment_number) + "\n")
