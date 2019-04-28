#!/usr/bin/env python

# FUNCTIONS FOR CREATION OF THE REPO REPORT


def count_lines(file_input):
    """
    Function to count lines of a given file
    """

    # Open the file, loop through the lines and add one to the count in each line
    count = 0
    with open(file_input, "r") as infile:
        for line in infile:
            count += 1
    # Return the line count
    return count


def report_repo():
    """
    Function to create a repo list for the repository
    """

    # Create a dictionary about the features of the genomic sequences
    # {seqKey:[seqID,lenght,conitg_count,N_count,[{file:vcf},{file:gff},{file:sam}]]}
    genome_report = {}
    # Create a dictionary showing the relationship between the seqID and the seqKey  {seqID:seqKey}
    ID_key = {}
    # Open the seqID file and add the seqID
    with open("./Genome/SeqIDs.txt", "r") as seqID_file:
        # Loop through the lines
        seq_count = 0
        for line in seqID_file:
            # Update the seqKey
            seq_count += 1
            seqKey = "Sequence" + str(seq_count)
            # Add the line to the dict as a key and start the values
            # {seqKey:[seqID,lenght,conitg_count,N_count,[{file:vcf},{file:gff},{file:sam}]]}
            genome_report[seqKey] = [line[1:].rstrip(), "", 0, 0, [{}, {}, {}]]
            # Add it to the ID_key
            ID_key[(line[1:].rstrip()).split(" ")[0]] = seqKey
    # Close the file
    seqID_file.close()
    # Open the map file
    with open("./Genome/Map.txt", "r") as map_file:
        # Loop through the lines
        for line in map_file:
            # If the line is a new sequence
            if (line[0] == ">"):
                line = line.rstrip().split("\t")
                seqKey = line[0][1:]
                # Add the length to the corresponding sequence in the dictionary
                genome_report[seqKey][1] = line[1]
            # Otherwise it is a sequence feature
            else:
                # Split the line
                line = line.rstrip().split("\t")
                # If this is a gap add it to the gap lenght
                if (line[0] == "N"):
                    genome_report[seqKey][3] += (int(line[1].split()
                                                     [1]) - int(line[1].split()[0]))
                # Otherwise this is a contig, add it to the count
                else:
                    genome_report[seqKey][2] += 1
    # Close the file
    map_file.close()
    # Open the repomap and save the dataset subfile into a dictionary
    # {dataset:[[filename,subfile_dir],[filename,subfile_dir],]}
    subfiles_dict = {}
    subfiles_dict["Alignment"] = []
    subfiles_dict["Variants"] = []
    subfiles_dict["Annotation"] = []
    with open("./RepoMap.txt", "r") as repomap:
        # Loop through the lines
        for line in repomap:
            # Split the line
            line = line.split("\t")
            # If the line is not a genome dataset
            if (line[1] != "Genome"):
                # Add the subfile directory to the dictionary
                subfiles_dict[line[1]].append([line[0], line[2]])
    # Close the file
    repomap.close()
    # Loop through the different datasets  {dataset:[[filename,subfile_dir],[filename,subfile_dir],]}
    for dataset in subfiles_dict.keys():
        # Loop through the subfiles of the dataset
        for subfile in subfiles_dict[dataset]:
            # Open the file map of the subfile
            with open(subfile[1] + "/SeqIDs_Map.txt", "r") as map_file:
                # Start line count
                line_count = 0
                # Loop through the lines of the map
                for line in map_file:
                    # If the line indicates a file, count its lines and add it to the line count
                    if (line[0] != ">"):
                        # result=subprocess.check_output
                        # ("wc -l "+subfile[1]+line.rstrip().split("\t")[0], shell=True).rstrip()
                        line_count += count_lines(subfile[1] +
                                                  line.rstrip().split("\t")[0])
                    # Otherwise it is a new sequence (check that this is not the first sequence)
                    elif (line_count != 0):
                        # Check the seqID is in the genome assembly
                        if (seqID in ID_key.keys()):
                            # Add the amount of lines to the dicionary of features of the genomic sequences
                            # {seqKey:[seqID,lenght,conitg_count,N_count,[{file:vcf},{file:gff},{file:sam}]]}
                            # Use the ID_key dictionary  {seqID:seqKey}
                            # If the dataset is Alignment
                            if (dataset == "Alignment"):
                                genome_report[ID_key[seqID]
                                ][4][2][subfile[0]] = line_count
                            # If it is annotaiton
                            elif (dataset == "Annotation"):
                                genome_report[ID_key[seqID]
                                ][4][1][subfile[0]] = line_count
                            # Otherwise it is variants
                            else:
                                genome_report[ID_key[seqID]
                                ][4][0][subfile[0]] = line_count
                            # Store the new seqID and re-start the line count
                            line_count = 0
                            seqID = line.split("\t")[0][1:]
                    # Otherwise it is the first sequence
                    else:
                        seqID = line.split("\t")[0][1:]
                # Add the last subfile left out of the loop.
                # Check the seqID is in the genome assembly
                if (seqID in ID_key.keys()):
                    if (dataset == "Alignment"):
                        genome_report[ID_key[seqID]][4][2][subfile[0]] = line_count
                    # If it is annotaiton
                    elif (dataset == "Annotation"):
                        genome_report[ID_key[seqID]][4][1][subfile[0]] = line_count
                    # Otherwise it is variants
                    else:
                        genome_report[ID_key[seqID]][4][0][subfile[0]] = line_count
            # Close the map file
            map_file.close()
    # Return the dictionary {seqKey:[seqID,lenght,conitg_count,N_count,[{file:vcf},{file:gff},{file:sam}]]}
    return genome_report
