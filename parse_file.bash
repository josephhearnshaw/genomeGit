#!/bin/bash

##Arguments: dataset file threads file_size

#Obtain the script location
source="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#If the first argument indicates corresponds with genome dataset, check if a dataset already exists
if [ "$1" = "Genome" ]; then
   #If there is a Genome folder already in the repo, evaluate if there are dependent files associated to it.
   if [ -d "Genome" ]; then
      #Make update=1
      update=1
      #If there are dependent files that can be affected, this is not just parsing the genome file, this is an update.
      if [ -d "Variants" ] || [ -d "Annotation" ] || [ -d "Alignment" ]; then
         #Before starting the auto-update, it is necessary that the user has a valid version of mummer and tabix installed.
         command -v mummer >/dev/null 2>&1 || { echo "GenomeGit requires MUMmer for liftover, but no MUMmer installation was found in this machine.  Now aborting."; exit 0;}
         command -v tabix >/dev/null 2>&1 || { echo "GenomeGit requires Tabix for liftover, but no Tabix installation was found in this machine.  Now aborting."; exit 0;}
         #If for any reason a temporal temporary_directory directory already exists, delete it.
         if [ -d "temporary_directory" ]; then
            rm -r temporary_directory
         fi
         #Call the updating script, inform of which datasets to be updated. Arguments: file threads size tlength.
         python $source/update.py $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12}
         #Delete the temporary directory
         rm -r ./temporary_directory
      fi
      #When this is done, the genome parsing can go as normal. Inform the user.
      echo ""
      echo "A genome dataset has been detected in the repository. Now overwriting with new dataset."
      rm -r $1
      #Otherwise make update=0
   else
      update=0
   fi
   #Create a new folder and get inside
   mkdir $1
   cd $1
   ##Create seqIDs file inside the directory
   grep -o -E "^>.*" $2 > SeqIDs.txt
   #Determine the number of sequences
   sequence_number=$(grep -c ">" $2)
   if [ "$sequence_number" = "0" ]; then
      echo ""
      echo "***INPUT ERROR: The input file contains no sequences. Please make sure to input a valid FASTA file.***"
      echo ""
      rm -r ../$1
      exit 0
   fi
   #Inform the user
   echo "Now processing $2 containing a Genome dataset with $sequence_number sequences:"
   date
   #Call python script to parse all the sequences in the input file. Arguments: class filename filesize update
   python $source/parse.py $1 $2 $4 $update
   #Inform the user
   echo "Genome parsing finished:"
   date
   echo ""


   #If it does not correspond with a genome dataset, it is a dependent file.
else
   #If this is the first file of the dataset, create a new directory
   if ! [ -d $1 ]; then
      mkdir $1
   fi
   #Save the file name
   if [[ $2 == *".bam" ]]; then
      fileBasename=$(basename $2)
      FileName=${fileBasename/.bam/'.sam'}
   else
      FileName=$(basename $2)
   fi

   #If a file with the same name is stored, delete its contents and make update=1
   if [ -d ./$1/$FileName ]; then
      rm -r ./$1/$FileName/*
      update=1
      #Create a new directory inside the dataset directory, if it does not exist, and make update=0
   else
      mkdir ./$1/$FileName
      update=0
   fi
   #Enter the direcotry
   cd ./$1/$FileName
   #Inform the user
   echo "Now processing $1 dataset contained in $FileName:"
   date
   #Call the python script ot parse the dependent file. Arguments: class file size update
   python $source/parse.py $1 $2 $4 $update
   #Inform the user
   echo "Processing completed:"
   date
   echo ""
fi
