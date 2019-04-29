# genomeGit 3.0: A distributed version control system for fast and efficient updating of genome assembly data.
genomeGit 3.0 is a distributed version control system which utilizes Git for the quick storage and management of genomic data. It enables quick 'lifting-over' of genomic depedent files, while being storage efficient. GenomeGit 3.0 can currently deal with the following datasets:

* Genome assemblies (FASTA)
* Variant Calling Files (VCF)
* Annotation files (GFF/GFF3)
* Alignment files (SAM/BAM)

GenomeGit 3.0 is currently only compatiable with Unix Operating Systems.

# What's changed with genomeGit 3.0?
genomeGit 3.0 is now able to make use of BAM files. There is also a new hybrid alignment option avaliable using MashMap2 and Nucmer4. 
In addition, genomeGit 3.0 is faster, more accurate, can handle inversions, and handle both splits and merges within genomic data. 

### Prerequisites for genomeGit 3.0
The following dependencies are required: 
* [Python v 2.7+](https://www.python.org/)
* [Git 2.7.4](https://git-scm.com/downloads)
* [MUMmer 4.0](https://mummer4.github.io/)
* [MashMap 2.0](https://github.com/marbl/MashMap)
* [Tabix 1.9](http://www.htslib.org/doc/tabix.html)

genomeGit 3.0 also makes use of the following Python modules: 
 * [pytabix 0.1](https://pypi.org/project/pytabix/)
 * [pyfaidx 0.5.5.2](https://pypi.org/project/pyfaidx/)

### Installation
In order to be able to run the program anywhere on your system, run ```genomegit_install```. This wil lcreate a symlink for genomegit in the /usr/bin directory. You may need to make the main script executable using ```chmod u+x <path_to_genomeGit>/genomeigt```.
#### Temporary solution
You can temporarily append the *genomegit* directory to the *$PATH* variable by executing ```PATH=$PATH:directory```, where ```directory``` repersents the location of the *genomeGit* directory. The user may need to make the main script executable using ```chmod u+x <path_to_genomeGit>/genomegit```. 

#### Uninstalling genomeGit 3.0

To uninstall genomeGit, remove the genomeGit directory from your PC and then remove the symlink using the following: ```rm -rf ./genomegit``` and  then ```sudo rm /usr/bin/genomegit```.

## Running genomeGit

To display the genomeGit welcome message, execute ```genomegit```.

Please note that genomeGit can make use of Git commands, executed as follows: ```genomegit <git_command>```. [See git documentation for further information on how to use git.](https://git-scm.com/doc)

To get the list of all available genomeGit commands, execute ```genomegit help```.

### 1. Initializing the repository
The repository can be initialized by executing ```genomegit init```, creating a *.gnmgit* directory. 

This repository will store all of your genomic data and the *.git* repository. You can clone an existing repository by executing ```genomegit clone <url>```.


### 2. Adding new files to the repository

To add files into the repository, execute ```genomegit add <file>``` and ```genomegit commit -m <message>```. 
#### Additional ```add``` arguments for lift-over
Additional arguments can be passed to ```genomegit add <file>```, such at the number of threads (```--t=<x>``` or ```--thread=<x>```),  and the aligner you wish to use (```--a=<1 or 2>``` or ```--aligner=<1 or 2>```, where 1 will run the hybrid aligner and 2 runs Nucmer4 only). 
#### Specific aligner-related flags
Flags specific to NUCmer4 (```--c=<x>``` or ```--mincluster=<x>```) can be used ([see the NUCmer documentation for information regarding these flags](http://mummer.sourceforge.net/manual/#nucmer)). 

Likewise, the flags ```--k=<x>``` or ```--kmer=<x>``` and ```--s=<x>``` or ```--segLength=<x>``` and ```--pi=<x>``` can be used for MashMap2 ([see the MashMap GitHub page for information regarding these flags](https://github.com/marbl/MashMap)). 

#### Splits and merges related flags
For the hyrid alignment, the user can also use the flag ```--ms=<x>``` as either 1 or 2, where 2 will also detect merges (but not splits). By default, it'll detect splits (```--ms=1```). 

Using only NUCmer4 will result in the detection of both splits and merges. 

#### Obtaining a report of lifted over assemblies
genomeGit 3.0 will automatically classify the file inputted and parse it into it's respective Git-compatible sub-files, within the Git repository. A summary of the charactersticis of the data within the repository can be visualised using the command ```genomegit report```. 

When a user already has a genome assembly present within the repository and wishes to update it, genomeGit 3.0 will automatically migrate the coordinates of the stored dependent files. This is called lift-over. This process can be computationally demanding and it is thus recommended to use the optional ```--t=<x>``` or ```--threads=<x>``` parameter to choose the number of threads used during lift-over.

### 3. Creating a remote repository
genomeGit 3.0 enable users to perform updates within thier local repository and push this to a central repository for all users to use. This can be performed by the command ```genomegit init --bare <remote_name>```.

#### Remote repository within the same machine
To access a remote repository, the remote repository address needs to be added. This can be done as follows: ```genomegit remote add <remote_name> <remote_location>``` where ```<remote_location>``` is the absolute path to the repository of interested, if located within the same machine it's executed from.

#### Remote repository from a server 
To obtain a repository from a server, a username and server IP address is required, i.e. ```your_usernmae@xxx.xxx.xx.x```. 
The ```<remote_name>``` parameter is the 'nick name' of the repository that the user provides when *pushing* and *pulling* from a local repository. 

To update the remote repository, the user can fetch the remote repositories data and *push* it into their local repository by typing ```genomegit pull <remote_name> <branch_name>```. Any changes that were introduced into the local repository can be pushed into the remote by executing ```genomegit push <remote_name>```. 


### 4. Assembly version log, checking out a Assembly version of interest, and listing all the files of a particular dataset
#### Switching to a stored assembly version
You can switch to any of your stored assembly versions by using ```genomegit checkout <commit_hash>```, where ```<commit_hash>``` repersents the SHA-1 commit hash. This can be obtained via executing ```genomegit log```. 

#### Reconstructing Git-compatible files
To reconstruct any Git-compatible files, such as the extracted VCF data, execute ```genomegit get --dataset --sequence --region --commit-hash --message <filename>```, where ```<filename>``` and ```<get --dataset>``` are required arguments. 
The ```<filename>``` argument will require the user to enter the file of interest that they wish to reconstruct. 

If a user wishes to reconstruct a Variants datatype, they can execute ```genomegit get --dataset=Variants```. The optional parameters ```--sequence``` and ```--region``` can be used to extract regions of a sequence which are contained within the file of interested. The region must be specified as a range in the form of two integers, seperated by a dash ("-"), e.g. ```--region=1-5000```. If the file of interest belongs to a previous repository version, then ```--commit-hash``` and ```--message``` can be used to specify the version's commit hash, or its message. 

#### List all the files present in the repository
The command ```genomegit list``` will obtain a list with all the file names present within the repository. 
To revert back to the main branch, you can type ```genomegit checkout <branch name>```. 

### 5. Reporting changes that occurred between versions
To view the difference between two versions of the genomic data present within the repository, execute ```genomegit diff --message=<message> <hash1> <hash2>```, where ```<hash1> <hash2>``` repersent the hashes of the commits given to the user following lift-over. 

The user can alternatively use ```--message=<message>``` if they used a commit message instead. The number of threads can be provided using ```genomegit diff --threads``` or ```--t=<x>```, as previously described. This might prove useful when comparing non-consecutive versions, as comparisons of assembly versions may need to be performed. ```genomegit log``` will provide a list of commits with their hashes.

The user can also execute ```genomegit report``` to obtain a report output *GenomeGit_Report.txt*, which contains the all the information regarding the differences between assembly versions.


## Example commands
#### Installing genomeGit
```bash /path/to/directory/genomegit_install```
#### Initiating an empty repository, or cloning it
```genomegit init```

```genomegit clone username@138.250.31.98:/home/user/repository```
#### Creating a remote repository, connecting to it, and push and pull 
```genomegit init --bare <remote_name>```

```genomegit remote add origin username@138.250.31.98```

```genomegit pull origin master```

```genomegit push origin```
#### Adding and committing files into the repository
```genomegit add /path/to/directory/*``` 

```genomegit add --t=32 --a=2 --c=3000 NewAssembly.fa```

```genomegit commit -m "My_First_Commit"```
#### Extracting a file out of the repository
```genomegit get MyAssembly.fa```

```genomegit get MyAnnotation.gff --sequence=Ch02 --region=1-2000000```

```genomegit get --dataset=Variants --sequence=Ch03 --region=1-4000000 --message=My_First_Commit```
#### Produce reports and summaries of changes between versions
```genomegit report	--commit-hash 352d69d2327dc95b58d6ec10130366ddc760bd3d```

```genomegit diff --message My_First_Commit My_Second_Commit```
