#!/usr/bin/env python

#	PYTHON SCRIPT TO PRINT A LIST OF FILES CONTAINED IN THE REPOSITORY.

# Make imports
import sys

# Initate dataset_dict dictionary {dataset:[file_name,...]}
dataset_dict = {"Genome": [], "Annotation": [], "Variants": [], "Alignment": []}

# Open the repomap in read mode
with open(str(sys.argv[1]), "r") as repomap:
    for line in repomap:
        # Split the line
        line = line.split("\t")
        # Save the name of the file on its corresponding dataset.
        dataset_dict[line[1]].append("{} ({})".format(
            line[0], line[4].rstrip()))

# Print a table with the files for each dataset. Loop through all the datasets.
for dataset in dataset_dict.keys():
    # Print dataset table header
    print("\n\t------\n\t|\n\t| {} files:\n\t|".format(dataset))
    # If the dataset has no files
    if(len(dataset_dict[dataset]) == 0):
        # Inform the user the dataset is empty
        print("\t|\tThis dataset contains no files.")
    else:
        # Loop through the files contained in the dataset
        for file_name in dataset_dict[dataset]:
            print("\t|\t{}".format(file_name))
    print("\t|\n\t------")
