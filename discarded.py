#!/usr/bin/env python

# USAGE: The first argument must indicate the file of interest.

# Make imports
import os
import shutil
import sys

# Inititate variables
valid = False
# Open the repomap and look for the file of interest
with open("./RepoMap.txt", "r") as repomap:
    # Loop through the lines of the repomap
    for line in repomap:
        # Split the line
        line = line.split("\t")
        # Check if this is the file of interest
        if (line[0] == str(sys.argv[1])):
            # If the file is a genome dataset, it was invalid user input
            if (line[1] == "Genome"):
                print("***WARNING: The provided file corresponds with a genome dataset:" + str(sys.argv[1]) +
                      "\nGenome datasets do not have discarded entries to be retrieved.")
                sys.exit()
            # The file name provided is valid
            valid = True
            # Save the directory and filename
            directory = line[2]
            filename = line[0]
            # No need to keep looking, break
            break
    # If the file was present in the repo, check if it has discarded entries
    if (valid):
        # Check if there is a discarded file
        if (os.path.isfile(directory + "/Discarded")):
            # Copy the discarded outside the repository
            shutil.copyfile(directory + "/Discarded", "../Discarded_" + filename)
        # Otherwise inform the user
        else:
            print("The selected file does not have any discarded entries: " + str(sys.argv[1]))
    else:
        print("***WARNING: The provided file was not found in the repository: " + str(sys.argv[1]))
