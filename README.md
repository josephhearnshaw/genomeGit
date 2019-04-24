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
### Installation
In order to be able to run the program, the location of the *GenomeGit* directory is required to be in the *$PATH* variable. You can install it by either one of the following two methods:
#### Temporally append the location of the directory into the *$PATH* variable (not recommended)
Temporally append the *GenomeGit* directory to the *$PATH* variable by executing ```PATH=$PATH:directory```, where ```directory``` repersents the location of the *GenomeGit* directory. The user may need to make the main script executable using ```chmod u+x <path_to_GenomeGit>/genomegit```. 

However, please note that using this method is only a short-term solution; you will need to perform this operation eachtime a terminal is closed and reopened.
#### Permanently append the location of the directory into the *$PATH* variable (recommended)
In order to easily append  *GenomeGit* directory to the *$PATH* variable permanenlty, a ```genomegit_install``` executable file has been provided. Simply execute this file to create the symlink to genomegit in the ```/usr/bin/ directory```. You may need to make genomegit executable via chmod (i.e. ```chmod u+x /usr/bin/genomegit```).

To uninstall genomeGit, remove the genomeGit directory from your PC and then remove the symlink using the following: ```sudo rm /usr/bin/genomegit```.

## Running GenomeGit

To display the GenomeGit welcome message type ```genomegit```.

GenomeGit adapts and uses Git commands, which can be ran by typing ```genomegit <git_command>```. [See git documentation for more information on how to use git.](https://git-scm.com/doc)

To get the list of genomeGit commands, type ```genomegit help```.


### 1. Initializing the repository
The repository can be initialized by typing ```genomegit init```, which will create a *.gnmgit* directory. 

This repository will store all of your genomic data and the *.git* repository. You can clone an existing repository by typing ```genomegit clone <url>```.


### 2. Adding new files to the repository
To add files into the repository, type ```genomegit add <file>``` and ```genomegit commit -m <message>```. Additional arguments can be passed to ```genomegit add <file>```, such at the number of threads (```--t=<x>``` or ```--thread=<x>```), the aligner you wish to use (```--a=<1 or 2>``` or ```--aligner=<1 or 2>```) where 1 will run the hybrid aligner and 2 runs Nucmer4 only. 

Flags specific to Nucmer4 (```--b=<x>``` or ```--breaklen=<x>``` and ```--c=<x>``` or ```--mincluster=<x>```) can be used ([see the NUCmer documentation for information regarding these flags](http://mummer.sourceforge.net/manual/#nucmer)). Likewise, the flags ```--k=<x>``` or ```--kmer=<x>``` and ```--s=<x>``` or ```--segLength=<x>``` can be used in MashMap2. 

GenomeGit 3.0 will automatically classify the file inputted and parse it into it's respective Git-compatible sub-files, within the Git repository. A summary of the charactersticis of the data within the repository can be visualised using the command ```genomegit report```. 

When a user already has a genome assembly present within the repository and wishes to update it, GenomeGit will automatically migrate the coordinates of the stored dependent files. This is called liftover. This process can be computationally demanding and it is thus recommended to use the optional ```--t=<x>``` or ```--threads=<x>``` parameter to choose the number of threads used during liftover.

### 3. Creating a remote repository
Given that GenomeGit 3.0 is a Distributed Version Control System, it is able to enable users to perform updates in thier local repository and then push this to a central repository for all users to use. This can be performed by the command ```genomegit init --bare <remote_name>```.

To access a remote repository, the remote repository address needs to be added. This can be done as follows: ```genomegit remote add <remote_name> <remote_location>``` where ```<remote_location>``` is the absolute path to the repository of interested, if located within the same machine it's executed from, or the username and server IP address if stored on a server i.e. ```your_usernmae@xxx.xxx.xx.x```. The ```<remote_name>``` parameter is the 'nick name' of the repository that the user provides when *pushing* and *pulling* from a local repository. To update the remote repository, the user can fetch the remote repositories data and *push* it into their local repository by typing ```genomegit pull <remote_name> <branch_name>```. Any changes that were introduced into the local repository can be pushed into the remote by typing ```genomegit push <remote_name>```. 


### 5. Version log, checking out a version of interest and listing all the files of a particular dataset
User can switch to any of the stored assembly versions at any moment by typing ```genomegit checkout <commit_hash>```, where ```<commit_hash>``` stands for the SHA-1 commit hash of the version of interest. To obtain this hash, a review of versions can be viewed by typing ```genomegit log```. Additionally, in order to reconstruct the git-compatible objects containing the information for any of the datasets contained in the repository, type ```genomegit get --dataset --sequence --region --commit-hash --message <filename>```. The mandatory argument ```<filename>``` stands for the name of the file that the user wants to obtain back. Additionally, optional parameters ```--sequence``` and ```--region``` can be used in order to extract only a particular region of a sequence contained in the file of interest (region must be specified in form of two integers separated by a "-", for example ```--region 1-1000```). If the file to be retrieved belongs to a previous version of the repository, optional parameters ```--commit-hash``` and ```--message``` may be used alternatively in order to specify the version's commit hash or message respectively (for example ```--message=Version_A```). Using the command ```genomegit list``` it is possible to obtain a list with the names of the files stored in the repository. In order to revert back to the main branch, type ```genomegit checkout <branch_name>```.

### 6. Report of changes occurred between versions.
The user is able to see the differences between two chosen versions of the data contained in the repository at any time. This can be done by typing ```genomegit diff --threads --message <hash1> <hash2>```, where  ```<hash1> <hash2>``` stand for the hashes of the commits that correspond to the chosen genome versions; or alternatively the commit message if ```--message``` optional parameter is provided. Optional parameter ```--threads``` can be used to provide the number of threads to be used for the obtention of the sequence aligment (this alignment migth be necessary when comparing non-consecutive versions). The list of commits and their hashes can be seen by typing ```genomegit log```. As a result of this command, a summary of the differences between versions will be displayed.

## Some common commands
#### Install GenomeGit
```bash /path/to/directory/genomegit_install```
#### Initiate empty repository or clone it
```genomegit init```

```genomegit clone username@123.456.78.9:/home/user/repository```
#### Create a remote repository, connect to it, psuh and pull changes
```genomegit init --bare RemoteExample```

```genomegit remote add origin username@123.456.78.9```

```genomegit pull origin master```

```genomegit push origin```
#### Add and commit files into the repository
```genomegit add /path/to/directory/*```

```genomegit add --t=8 --a=2 --b=1 --c=1500 NewAssembly.fasta```

```genomegit commit -m "First_Commit"```
#### Extract a file out of the repository
```genomegit get MyAssembly.fasta```

```genomegit get MyAnnotation.gff --sequence=Ch01 --region=1-1000000```

```genomegit get --dataset=Variants --sequence=Ch01 --region=1-1000000 --message=First_Commit```
#### Produce reports and summaries of changes between versions
```genomegit report	--commit-hash e7201f04231d039a48ea41f1f786d4f447361175```

```genomegit diff --message First_Commit Second_Commit```




 
