#!/usr/bin/env python

"""	USAGE: The first argument must indicate the name of the file to be deleted.
	It will be deleted both from genomegit and git staging areas. """

# Make the imports
import sys
import shutil
from subprocess import Popen

# Check if the provided file name is contained in the Repomap. Open the repomap withth read permit
with open("./RepoMap.txt", "r") as repomap:
    # Initiate valid variable (false if the file is not found in the repo)
    valid = False
    # Initiate empty file list
    file_list = []
    # Loop through the lines of the file
    for line in repomap:
        # Split the line
        line = line.split("\t")
        # If the name of the current field matches the one the user wants to delete, proceed to delete it
        if (line[0] == str(sys.argv[1])):
            # Delete the file both from the git and genomegit staging areas.
            ShellCommand = Popen("git reset " + line[2], shell=True).wait()
            shutil.rmtree(line[2])
            valid = True
        # Otherwise the line must be kept
        else:
            file_list.append("\t".join(line))
# Close the repomap
repomap.close()

# If the file name was valid, proceed to rewrite the repomap without it
if (valid):
    # Open the repomap withth write permit
    with open("./RepoMap.txt", "w") as repomap:
        # Loop through the file list and write it into the repomap
        for line in file_list:
            repomap.write(line)
    # Close the repomap
    repomap.close()

# If the file name was not valid, inform the user
else:
    print("***WARNING: The provided file was not found in the repository. Now aborting.")
