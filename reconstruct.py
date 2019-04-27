#!/usr/bin/env python

# Make the imports
import sys
from reconstruct_functions import reconstruct_dataset, obtain_file, extract_dataset
from subprocess import Popen

# Check if the user wants to extract an entire dataset
if(str(sys.argv[2]) != "0"):
    # Check if the dataset is invalid
    if(str(sys.argv[2]) != "Genome" and str(sys.argv[2]) != "Variants" and
            str(sys.argv[2]) != "Alignment" and str(sys.argv[2]) != "Annotation"):
        print("***WARNING: The provided dataset was not recognised: " +
              str(sys.argv[2]) + "\nNow aborting.")
    # Otherwise it is a valid dataset
    else:
        extract_dataset(dataset=str(sys.argv[2]), seqID=str(
            sys.argv[3]), region=str(sys.argv[4]), convertToBam=str(sys.argv[5]))

# Otherwise the user only wants to reconstruct a file of intererst
else:
    # Obtain file name
    file_name = str(sys.argv[1])
    # Obtain dataset and line size
    dataset, directory, line_size, size = obtain_file(
        target=file_name, mode="filename")
    # If the file is not found in the repository, exit the program
    if(dataset == "error"):
        print("***WARNING: No file was found in the repository with the provided name: " +
              file_name + "\nNow aborting.")
        sys.exit()
    # If the dataset to be reconstructed is a genome dataset. If the line size is not specified, default it to 60
    if(dataset == "Genome"):
        # Create index variable
        index = 0
        # Inform the user
        print("\nNow reconstructing Genome dataset contained in file: " +
              file_name + ".fa")
        # Use reconstruct_dataset function to recreate the file
        reconstruct_dataset(size=int(line_size), directory=directory,
                            output_file="../" + file_name, mode="Genome", seqID="0", region="0")
        # Inform the user
        print("\nReconstruction completed.")
    # Otherwise the dataset to be reconstructed is a dependent file
    else:
        # Reconstruct the variants dataset
        if(dataset == "Variants"):
            # Inform the user
            print("\nNow reconstructing Variants dataset contained in file: " + file_name)
            reconstruct_dataset(size=1, directory=directory, output_file="../" + file_name +
                                ".temporary", mode="Variants", seqID=str(sys.argv[3]), region=str(sys.argv[4]))
            # Sort the obtained file
            ShellCommand = Popen('(grep "^#" ../' + file_name + '.temporary; grep -v "^#" ../' + file_name +
                                 ".temporary | sort -V -k1,1 -k2,2n) | cat > ../" + file_name +
                                 "; rm ../" + file_name + ".temporary", shell=True).wait()
            # Inform the user
            print("\nReconstruction completed.")
        # Reconstruct the annotation dataset
        elif(dataset == "Annotation"):
            # Inform the user
            print(
                "\nNow reconstructing Annotation dataset contained in file: " + file_name)
            reconstruct_dataset(size=1, directory=directory, output_file="../" + file_name +
                                ".temporary", mode="Annotation", seqID=str(sys.argv[3]), region=str(sys.argv[4]))
            # Sort the obtained file
            ShellCommand = Popen('(grep "^#" ../' + file_name + '.temporary; grep -v "^#" ../' + file_name +
                                 ".temporary | sort -k1,1 -k4,4n) | cat > ../" + file_name +
                                 "; rm ../" + file_name + ".temporary", shell=True).wait()
            # Inform the user
            print("\nReconstruction completed.")
        # Reconstruct the alignment dataset
        elif(dataset == "Alignment"):
            # Inform the user
            print(
                "\nNow reconstructing Alignment dataset contained in file: " + file_name)
            reconstruct_dataset(size=1, directory=directory, output_file="../" + file_name,
                                mode="Alignment", seqID=str(sys.argv[3]), region=str(sys.argv[4]), convertToBam=str(sys.argv[5]))
            # Inform the user
            print("\nReconstruction completed.")

        else:
            print("***INTERNAL ERROR: Dataset not recognised.")
