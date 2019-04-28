#!/usr/bin/env python

# USAGE: The first argument must indicate the dataset to be parsed: Genome, Variants, Annotation.
# The second argument must indicate the path to the file be parsed.

####
# ***PART 1. IMPORT FUNCTIONS (SUBROUTINES)***
####

# Make the imports
import sys
from parse_functions import parse_dataset

#####
# ***PART 2. PARSE DATASET***
#####

parse_dataset(dataset=str(sys.argv[1]), input_path=str(
    sys.argv[2]), size=str(sys.argv[3]), update=int(sys.argv[4]))
