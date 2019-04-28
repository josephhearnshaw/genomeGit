#!/usr/bin/env python

# PYTHON FUNCTIONS USED FOR REPOSITORY UPDATES


def detect_updates(map_path):
    """
    Create detect_updates function to return a dictionary informing of the files to be updated
    {dataset:[[filename,directory,size],[...]]}
    """

    # Initiate empty variable and empty dictionary {dataset:[[filename.extension,directory,size],[...]]}
    empty = True
    ToUpdate = {"Genome": [], "Annotation": [], "Variants": [], "Alignment": []}

    # Open the repomap
    with open(map_path, "r") as repomap:
        # Loop through the lines of the repomap
        for line in repomap:
            # Split the line
            line = line.split("\t")
            # If the file is a correct dataset add the path of the file to the update dictionary and make empty false
            if(line[1] in ToUpdate.keys()):
                empty = False
                ToUpdate[line[1]].append([line[0], line[2], line[4].rstrip()])
            # Otherwise is not a correct dataset
            else:
                print(
                    "***INTERNAL ERROR*** DATASET NOT RECOGNIZED: {}".format(line[1]))
                return "empty"
    # Close the repomap
    repomap.close()
    # If there are no files to update return empty variable
    if(empty):
        return "empty"
    else:
        return ToUpdate
