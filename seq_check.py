#### install PyVCF, pysam, numpy, and argparse with pip ####
import sys
import vcf
import os
from pyfaidx import Fasta
import argparse
import pysam
import re
import numpy as np


def file_paths(vcf_sam_name, fasta_name):
    """
    Obtains the file paths of a VCF and fasta file path

    Parameters
        User input strings (vcf_sam_name & fasta_name)

    Returns
        Absolute file paths to both files
    """

    vcf_sam_name = os.path.abspath(vcf_sam_name)
    fasta_name = os.path.abspath(fasta_name)
    return(vcf_sam_name, fasta_name)


def capatalize(var):
    """
    Capatalizes the variable inputted

    Parameters
        String

    Returns
        Captalized String
    """

    var = var.upper()
    return var


def info_switch(switch):
    """
    Always default the switch as False if nothing is chosen

        Parameters
            switch (user argument --nosum)

        Returns
            summary bool switch
    """

    if switch is None:
        summary = False
    elif(capatalize(str(switch)) == "TRUE"):
        summary = True
    else:
        summary = False
    return summary


def getfileName(file):
    """
    Obtains the filename of the inputted files

    Parameters
        String of absolute path

    Returns
        Filename from absolute path
    """

    # Split by backslash, get the last index
    file = file.rsplit('/', 1)[-1]
    # Remove )" which is found in Pyfadix files
    file = re.sub('[")]', '', file)
    return file


######################
# Checking VCF files #
######################

def check_vcf_fasta(vcf_name, fasta_name, out_name, noSummary):
    """

    Opens the VCF file and check the bases agaisnt the FASTA file.

    Parameters
        Strings; vcf file path, fasta file path, output name, and the summary switch (see info_switch)

    Returns
        Output file containing information of sequence matches and the number of correct matches

    """

    vcf_reader = vcf.Reader(open('{}'.format(vcf_name), 'r'))
    # Initate the counters
    correctCounter = errorCount = seqLenCount = 0
    if not noSummary:
        print("\nPrinting a summary output only.\n")
    print("\nNow comparing your VCF ({}) and FASTA ({}) file... \n\nWriting out to {}.txt".format(
        getfileName(vcf_name), getfileName(str(fasta_name)), out_name))

    original_stdout = sys.stdout
    file_out = open("{}.txt".format(out_name), 'w')
    sys.stdout = file_out
    # Obtain all records in the VCF file
    for record in vcf_reader:
        # Capatlize to ensure that no match doesn't occur due to differences in characters being lower or upper case
        vcf_base = capatalize(record.REF)
        # Obtain all names from Pyfadix
        for chromosomes in fasta_name.keys():
            # N.b. VCF starts at position 1 (N+1) and pyfaidx starts at index 0.
            fasta_base = capatalize(
                fasta_name[record.CHROM][record.POS - 1].seq)
            # If there are matches at the chromosome then the correct count must increase
            if(chromosomes == record.CHROM and vcf_base == fasta_base):
                if noSummary:
                    print("Base {} (chromosome {}, position {}) is correct".format(
                        vcf_base, record.CHROM, record.POS))
                correctCounter += 1
            # If there aren't any matches and the length of the seuqnece is less than 1, it's an error
            if(chromosomes == record.CHROM and vcf_base != fasta_base and len(vcf_base) <= 1):
                if noSummary:
                    print("***Base {} (chromosome {}, position {}) is INCORRECT. Found {} instead.***".format(
                        vcf_base, record.CHROM, record.POS, fasta_base))
                errorCount += 1
            # If there aren't any matches and the length of the seuqnece is greater than 1, it's due to sequence length differences
            if(chromosomes == record.CHROM and vcf_base != fasta_base and len(vcf_base) > 1):
                if noSummary:
                    print("*****Seq {} (chromosome {}, position {}) is too long and will not match. Found base {} in Fasta file.*****".format(
                        vcf_base, record.CHROM, record.POS, fasta_base))
                seqLenCount += 1

    print("\n***** Number of correct reads is {}*****\n*****Number of incorrect reads is {}*****\n**Number of reads "
          "with sequence length > 1 in VCF file is {}**".format(
              correctCounter, errorCount, seqLenCount))

    sys.stdout = original_stdout
    # File writing complete, close stdout.
    file_out.close()
    print("\nAll finished! Your output is in {}.txt".format(out_name))


######################
# Checking SAM files #
######################

def check_sam_file(sam_name, fasta_file, out_name, noSummary):
    """
         Opens the SAM file and check the Sequences agaisnt the FASTA file Sequences.

         Parameters
             Strings; SAM file path, fasta file path, output name, and the summary switch (see info_switch)

         Returns
             Output file containing information of sequence matches and the number of correct matches.

    """

    # Open the samfile using PySam
    samfile = pysam.AlignmentFile(sam_name, "r")
    # Initate the counters and sequence identity list
    correctCounter = errorCount = 0
    seq_ident_list = []

    if not noSummary:
        print("\nPrinting a summary output only.\n")
    print("\nNow comparing your SAM ({}) and FASTA ({}) file... \n\nWriting out to {}.txt".format(
        getfileName(sam_name), getfileName(str(fasta_file)), out_name))

    original_stdout = sys.stdout
    file_out = open("{}.txt".format(out_name), 'w')
    sys.stdout = file_out
    # N.b pos -> reference start, aend -> reference end -
    # - qstart -> query start and qend -> query end

    # Obtain all the chromosome names from the fasta file
    for chromosomes in fasta_file.keys():
        # Obtain all the names from the sam file
        for names in samfile.references:
            # Get every read in the sam file
            for read in samfile.fetch():
                samfile_seq = capatalize(str(read.seq))
                # Sequence matches
                if(chromosomes == names):
                    # Perform as default, or if chosen, otherwise just print the summary output
                    if noSummary:
                        try:
                            # Capatlize to ensure that no match doesn't occur due to differences in characters being lower or upper case
                            fasta_seq = capatalize(
                                str(fasta_name[chromosomes][read.pos:read.aend].seq))
                            # Split the sequences into characters
                            split_fasta = list(fasta_seq)
                            split_sam = list(samfile_seq)
                            if(samfile_seq == fasta_seq):
                                # Seq ID will always be 100% if sequences match
                                print("Correct seq for, {} , at start:, {}, end:, {}, with seq identity of, 100.00%".format(
                                    chromosomes, read.pos, read.aend))
                                correctCounter += 1
                            # Check if reversed sequences are the same
                            elif read.is_reverse:
                                rev_fasta_seq = fasta_seq[::-1]
                                if(samfile_seq == rev_fasta_seq):
                                    print("***REVERSED: seq for {} at start: {} end: {}***".format(
                                        chromosomes, read.pos, read.aend))
                        # Position errors -  check seq identity below:
                        except ValueError:
                            different_bases = [x for x, y in zip(
                                split_sam, split_fasta) if x != y]
                            # Check reversed sequences
                            if read.is_reverse:
                                rev_fasta_seq = list(fasta_seq[::-1])
                                different_bases = [x for x, y in zip(
                                    split_sam, rev_fasta_seq) if x != y]
                                seq_identity = (
                                    float(len(different_bases)) / float(len(samfile_seq))) * 100
                                seq_ident_list.append(seq_identity)
                                print("***Incorrect REVERSED seq for, {}, at start:, {}, end:, {}, with seq identity of, {:.2f}%***".format(
                                    chromosomes, read.pos, read.aend, seq_identity))
                                errorCount += 1
                            else:
                                # Obtain the differences between the fasta and sam sequences
                                # Obtain the sequence identity
                                seq_identity = (
                                    float(len(different_bases)) / float(len(samfile_seq))) * 100
                                seq_ident_list.append(seq_identity)
                                print("***Incorrect seq for, {}, at start: {}, end: {}, with seq identity of, {:.2f}%***".format(
                                    chromosomes, read.pos, read.aend, seq_identity))
                                errorCount += 1
                    # Else print the summary only
                    else:
                        try:
                            fasta_seq = capatalize(
                                str(fasta_name[chromosomes][read.pos:read.aend].seq))
                            split_fasta = list(fasta_seq)
                            split_sam = list(samfile_seq)
                            if(samfile_seq == fasta_seq):
                                correctCounter += 1
                        # Position errors -  check seq identity below:
                        except ValueError:
                            different_bases = [x for x, y in zip(
                                split_sam, split_fasta) if x != y]
                            if read.is_reverse:
                                rev_fasta_seq = list(fasta_seq[::-1])
                                different_bases = [x for x, y in zip(
                                    split_sam, rev_fasta_seq) if x != y]
                                seq_identity = (
                                    float(len(different_bases)) / float(len(samfile_seq))) * 100
                                seq_ident_list.append(seq_identity)
                                errorCount += 1
                            else:
                                # Obtain the differences between the fasta and sam sequences
                                # Obtain the sequence identity
                                seq_identity = (
                                    float(len(different_bases)) / float(len(samfile_seq))) * 100
                                seq_ident_list.append(seq_identity)
                                errorCount += 1

    samfile.close()
    # Obtain the average identity for all incorrect reads and format to 2dp
    identity_avg = np.array(seq_ident_list).mean()
    print("\n***** Number of correct reads is {}*****\n*****Number of incorrect reads is {} with average sequence identity of {:.2f}%*****".format(
        correctCounter, errorCount, identity_avg))
    # File writing complete, close stdout.
    sys.stdout = original_stdout
    file_out.close()
    print("\nAll finished! Your output is in {}.txt".format(out_name))


######################################################################################################################################################
# Initial                                                                       -                                                              START #
######################################################################################################################################################
parser = argparse.ArgumentParser(description="This will check to see if your VCF lift-over matches up with your FASTA file.\n\
                                 -v is your Vcf file. -f is your Fasta file, -o is your output filename.\n\nExample: "
                                 "python base_compare.py -v ./FruitFly/fruit_fly_75_to_1215.vcf -f ./FruitFly/GCA_002310755.1_ASM231075v1_genomic.fa\
                                  -o test_output")

parser.add_argument(
    '-v', '--vcf', help='Use the directory for your VCF file.', required=False)
parser.add_argument(
    '-f', '--FASTA', help='Use the directory for your FASTA file.', required=True)
parser.add_argument(
    '-o', '--output', help='Name your output file.', required=True)
parser.add_argument(
    '-s', '--SAM', help="Directory of your SAM file.", required=False)
parser.add_argument(
    '-nosum', '--SUM', help="Use this flag to show all reads alongside the summary information (-nosum true).", required=False)

# Obtain arguments
args = vars(parser.parse_args())
# Obtain output name for output file
out_name = args['output']
# Default summary switch is True, so always print everything, unless stated otherwise
noSummary = info_switch(args['SUM'])


if args['vcf'] is not None:
    vcf_name, fasta_name = file_paths(args['vcf'], args['FASTA'])
    print("\nIndexing your Fasta file, one moment...\n")
    fasta_name = Fasta(fasta_name)
    check_vcf_fasta(vcf_name, fasta_name, out_name, noSummary)
elif args['SAM'] is not None:
    sam_name, fasta_name = file_paths(args['SAM'], args['FASTA'])
    print("\nIndexing your Fasta file, one moment...\n")
    fasta_name = Fasta(fasta_name)
    check_sam_file(sam_name, fasta_name, out_name, noSummary)
else:
    pass
