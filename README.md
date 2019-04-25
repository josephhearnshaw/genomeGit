# GenomeGit: A distributed version control system for genome assembly data.
GenomeGit is a distributed version control system based on the creation of git repositories for the storage and management of genomic data. Currently, the program is able to deal with the following datasets:

* Genome assemblies (FASTA)
* Variant Calling Files (VCF)
* Annotation files (GFF/GFF3)
* Alignment files (SAM/BAM)

Liftover of dependent genomic VCF, SAM/BAM, and GFF files following the update of a genome assembly is supported. Currently, it is only compatible with Unix systems.

# What's new in GenomeGit 3.0?
GenomeGit 3.0 is now able to make use of BAM files. There is also a new hybrid alignment option avaliable using MashMap2 and Nucmer4. 

In addition, GenomeGit 3.0 is now faster, more accurate, can handle inversions, and handle both splits and merges within genomic data. 

## First steps
### Prerequisites
GenomeGit 3.0 requires installation of the following dependencies: 
* [Python v 2.7+](https://www.python.org/)
* [pip](https://github.com/pypa/pip)
* [Git](https://git-scm.com/downloads)
* [MUMmer 4.0](https://mummer4.github.io/)
* [MashMap 2.0](https://github.com/marbl/MashMap)
* [Tabix](http://www.htslib.org/doc/tabix.html)

It makes use of the following Python modules: 
 * [pytabix 0.1](https://pypi.org/project/pytabix/)
 * [pyfaidx 0.5.5.2](https://pypi.org/project/pyfaidx/)

The genomeGit installation will ask if you want to install thsi via pip. This assumes you have Python 2.7+ only. Please install these modules to your Python distro if you use a different Python distro.

### Installation
In order to be able to run the program, the location of the *GenomeGit* directory is required to be in the *$PATH* variable. You can install it by either one of the following two methods:
#### Temporally append the location of the directory into the *$PATH* variable (not recommended)
Temporally append the *GenomeGit* directory to the *$PATH* variable by executing ```PATH=$PATH:directory```, where ```directory``` repersents the location of the *GenomeGit* directory. The user may need to make the main script executable using ```chmod u+x <path_to_GenomeGit>/genomegit```. 

However, please note that using this method is only a short-term solution; you will need to perform this operation eachtime a terminal is closed and reopened.
#### Permanently append the location of the directory into the *$PATH* variable (recommended)
In order to easily append  *GenomeGit* directory to the *$PATH* variable permanenlty, a ```genomegit_install``` executable file has been provided. Simply execute this file to create the symlink to genomegit in the ```/usr/bin/ directory```. You may need to make genomegit executable via chmod (i.e. ```chmod u+x /usr/bin/genomegit```).

To uninstall genomeGit, remove the genomeGit directory from your PC and then remove the symlink using the following: ```sudo rm /usr/bin/genomegit```.

## Running GenomeGit

To display the GenomeGit welcome message execute ```genomegit```.

GenomeGit adapts and uses Git commands, which can be ran by execute ```genomegit <git_command>```. [See git documentation for more information on how to use git.](https://git-scm.com/doc)

To get the list of genomeGit commands, execute ```genomegit help```.


### 1. Initializing the repository
The repository can be initialized by executing ```genomegit init```, which will create a *.gnmgit* directory. 

This repository will store all of your genomic data and the *.git* repository. You can clone an existing repository by executing ```genomegit clone <url>```.


### 2. Adding new files to the repository
To add files into the repository, execute ```genomegit add <file>``` and ```genomegit commit -m <message>```. Additional arguments can be passed to ```genomegit add <file>```, such at the number of threads (```--t=<x>``` or ```--thread=<x>```), the aligner you wish to use (```--a=<1 or 2>``` or ```--aligner=<1 or 2>```) where 1 will run the hybrid aligner and 2 runs Nucmer4 only. 

Flags specific to Nucmer4 (```--b=<x>``` or ```--breaklen=<x>``` and ```--c=<x>``` or ```--mincluster=<x>```) can be used ([see the NUCmer documentation for information regarding these flags](http://mummer.sourceforge.net/manual/#nucmer)). Likewise, the flags ```--k=<x>``` or ```--kmer=<x>``` and ```--s=<x>``` or ```--segLength=<x>``` and ```--pi=<x>``` can be used in MashMap2 ([see the MashMap GitHub page for information regarding these flags](https://github.com/marbl/MashMap)). 

GenomeGit 3.0 will automatically classify the file inputted and parse it into it's respective Git-compatible sub-files, within the Git repository. A summary of the charactersticis of the data within the repository can be visualised using the command ```genomegit report```. 

When a user already has a genome assembly present within the repository and wishes to update it, GenomeGit will automatically migrate the coordinates of the stored dependent files. This is called liftover. This process can be computationally demanding and it is thus recommended to use the optional ```--t=<x>``` or ```--threads=<x>``` parameter to choose the number of threads used during liftover.

### 3. Creating a remote repository
Given that GenomeGit 3.0 is a Distributed Version Control System, it is able to enable users to perform updates in thier local repository and then push this to a central repository for all users to use. This can be performed by the command ```genomegit init --bare <remote_name>```.

To access a remote repository, the remote repository address needs to be added. This can be done as follows: ```genomegit remote add <remote_name> <remote_location>``` where ```<remote_location>``` is the absolute path to the repository of interested, if located within the same machine it's executed from, or the username and server IP address if stored on a server i.e. ```your_usernmae@xxx.xxx.xx.x```. The ```<remote_name>``` parameter is the 'nick name' of the repository that the user provides when *pushing* and *pulling* from a local repository. To update the remote repository, the user can fetch the remote repositories data and *push* it into their local repository by typing ```genomegit pull <remote_name> <branch_name>```. Any changes that were introduced into the local repository can be pushed into the remote by executing ```genomegit push <remote_name>```. 


### 4. Assembly version log, checking out a Assembly version of interest, and listing all the files of a particular dataset
You can switch to any of your stored assembly versions by using ```genomegit checkout <commit_hash>```, where ```<commit_hash>``` stands for the SHA-1 commit hash. This can be obtained by using ```genomegit log```. 

To reconstruct any Git-compatible files, such as the extracted VCF data, execute ```genomegit get --dataset --sequence --region --commit-hash --message <filename>```, where ```<filename>``` and ```<get --dataset>``` are mandatory arguments. The ```<filename>``` argument will require the user to enter the file of interest that they wish to reconstruct. 

If a user wishes to reconstruct a Variants datatype, they can execute ```genomegit get --dataset=Variants```. The optional parameters ```--sequence``` and ```--region``` can be used to extract regions of a sequence which are contained within the file of interested. The region must be specified as a range in the form of two integers, seperated by a dash ("-"), e.g. ```--region=1-5000```. If the file of interest belongs to a previous repository version, then ```--commit-hash``` and ```--message``` can be used to specify the version's commit hash, or its message. 

The command ```genomegit list``` will obtain a list with all the file names present within the repository. To revert back to the main branch, you can type ```genomegit checkout <branch name>```. 

### 5. Reporting changes that occurred between versions
To view the difference between two versions of the genomic data present within the repository, execute ```genomegit diff --message=<message> <hash1> <hash2>```, where ```<hash1> <hash2>``` repersent the hashes of the commits given to the user following lift-over. The user can alternatively use ```--message=<message>``` if they used a commit message instead. The number of threads can be provided using ```genomegit diff --threads``` or ```--t=<x>```, as previously described. This might prove useful when comparing non-consecutive versions, as comparisons of assembly versions may need to be performed. ```genomegit log``` will provide a list of commits with their hashes.

## Some common commands
#### Install GenomeGit
```bash /path/to/directory/genomegit_install```
#### Initiate empty repository or clone it
```genomegit init```

```genomegit clone username@123.456.78.9:/home/user/repository```
#### Create a remote repository, connect to it, push and pull 
```genomegit init --bare Remote```

```genomegit remote add origin username@123.456.78.9```

```genomegit pull origin master```

```genomegit push origin```
#### Add and commit files into the repository
```genomegit add /path/to/directory/*``` 

```genomegit add --t=8 --a=2 --b=1 --c=1500 NewAssembly.fasta```

```genomegit commit -m "First_Commit"```
#### Extract a file out of the repository
```genomegit get MyAssembly.fa```

```genomegit get MyAnnotation.gff --sequence=Ch02 --region=1-1000000```

```genomegit get --dataset=Variants --sequence=Ch02 --region=1-1000000 --message=my_first_commit```
#### Produce reports and summaries of changes between versions
```genomegit report	--commit-hash 430CE34D020724ED75A196DFC2AD67C77772D169```

```genomegit diff --message my_first_commit my_second_commit```




 
