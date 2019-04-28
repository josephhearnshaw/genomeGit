#!/usr/bin/env python

# PYTHON FUNCTIONS REQUIRED IN THE REPOSITORY AUTO-UPDATES PART I: ALIGNMENT OBTENTION.

# Make the imports
import hashlib
import pickle


def store_variables(variables, alignment_pickle):
    """
    Create a store_variables function to store the provided variables into a pickle
    """

    # Open a pickle
    with open("{}/pickle".format(alignment_pickle), "wb") as pickle_file:
        # Store the variable list
        pickle.dump(variables, pickle_file, protocol=1)
    # Close the pickle
    pickle_file.close()


def load_variables(alignment_pickle):
    """
    Create a load_variables function to load variabled from a pickle file
    """

    # Open the pickle
    with open(alignment_pickle, "rb") as pickle_file:
        # Load the alignment pickle.
        variables = pickle.load(pickle_file)
    # Close the pickle
    pickle_file.close()
    # Return the data
    return variables


def obtain_alignment_pickle(old_assembly, new_assembly):
    """
    Create a obtain_alignment_pickle to obtain the name of the pickle of two given assemblies
    """

    # Initiate the content list
    content = []
    # Open the new assembly and loop through the lines
    with open(new_assembly, "r") as new_assembly_file:
        for line in new_assembly_file:
            # Append the line to the content list
            line = line.rstrip()
            content.append(line)
    # Close the file
    new_assembly_file.close()
    # Obtain the hash
    new_assembly_hash = str(hashlib.sha1("".join(content)).hexdigest())
    # Reinititate the content list
    content = []
    # Open the new assembly and loop through the lines
    with open(old_assembly, "r") as old_assembly_file:
        for line in old_assembly_file:
            # Append the line to the content list
            line = line.rstrip()
            content.append(line)
    # Close the file
    old_assembly_file.close()
    # Obtain the hash
    old_assembly_hash = str(hashlib.sha1("".join(content)).hexdigest())
    # Return the combination of both hashes
    return str("./Delta/{}_{}".format(new_assembly_hash, old_assembly_hash))
