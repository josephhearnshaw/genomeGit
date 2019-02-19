# GenomeGit 2: A distributed version control system for genome assembly data.
GenomeGit is a distributed version control system based on the creation of git repositories for the storage and management of genomic data. Currently, the program is able to deal with the following datasets:

* Genome assemblies (FASTA)
* Variant Calling Files (VCF)
* Annotation files (GFF/GFF3)
* Alignment files (SAM)

Liftover of VCF, SAM and GFF files upon update of genomic assembly is supported. Currently, it is only compatible with Unix systems.

## First steps
### Prerequisites
The application requires installation of the following dependencies: 
* [Python v 2.7+](https://www.python.org/)
* [Git](https://git-scm.com/downloads)
* [MUMmer 4.0](https://mummer4.github.io/)
* [Tabix](http://www.htslib.org/doc/tabix.html)
### Installation
In order to be able to run the program, the location of the *GenomeGit* directory is required to be in the *$PATH* variable. This can be done in two different ways:
#### Temporally append the location of the directory into the *$PATH* variable (not recommended)
To temporally append *GenomeGit* directory to the *$PATH* variable execute ```PATH=$PATH:directory``` in the shell, where ```directory``` stands for the location of the *GenomeGit* directory. The user may need to make the main script executable using ```chmod u+x <path_to_GenomeGit>/genomegit```. Please note that using this method is only a short-term solution, as it will require the user to perform this operation everytime a terminal is closed and opened.
#### Permanently append the location of the directory into the *$PATH* variable (recommended)
In order to easily append  *GenomeGit* directory to the *$PATH* variable permanenlty, a ```genomegit_install``` executable file has been provided. When this file is executed from the shell, it will add a symlink in the ```/usr/bin/ directory``` pointing at the location of the GenomeGit directory. Additionally, the user may need to make the main script executable using ```chmod u+x /usr/bin/genomegit```. This method constitutes a long-term solution, as the user won't need to repeat this operation everytime a terminal is closed and opened or even if the computer is restarted. To unistall GenomeGit, simply remove the GenomeGit directory from your computer and remove the symlink executing ```sudo rm /usr/bin/genomegit```.

## Running GenomeGit

To display the GenomeGit welcome message type ```genomegit```.
GenomeGit adapts regular git commands that can be run by typing ```genomegit <git_command>```. [See git documentation for more information on how to use git.](https://git-scm.com/doc)
Additional commands have been created to deal with the requirements of files containing genomic data.
To get the list of the most common commands type ```genomegit help```.

### 1. Initializing the repository
The repository can be initialized by typing ```genomegit init```, which will create a *.gnmgit* directory. This directory will contain all the genomic data stored in the repository, including the *.git* repository itself. Additionally, it is possible to clone an existing repository using ```genomegit clone <url>```, which will create a *.gnmgit* with the cloned contents inside.

### 2. Recording changes made into the repository
Currently GenomeGit 2 is able to store and manage four types of datasets:

* Genome assembly (FASTA).
* Variant Calling Files (VCF).
* Annotation files (GFF/GFF3).
* Sequence Alignment Map files (SAM).

To add files into the repository, simply run ```genomegit add <file>``` and ```genomegit commit -m <message>```. GenomeGit 2 will automatically clasify the input file (FASTA, VCF, SAM or GFF) and parse it into sub-files to be stored in the Git repository. Please note that the GenomeGit 1 command ```genomegit parse <dataset> <file>``` is no longer on use in GenomeGit 2. Additionally, if all the files of interest are located in the same directory, the user may provide the location of this directory: ```genomegit add /path/to/directory/*```. By doing this, GenomeGit will automatically add all the files in the specified directory. It is possible to visualize a summary of the characteristics of the data contained in the repository at any given point by using the command ```genomegit report```. In those cases when a genome assembly is already stored in the repository (with its corresponding VCF, GFF or SAM files) and the user desires to add a newer version of the assembly, GenomeGit will automatically migrate the coordinates of any of the stored SAM, GFF or VCF files to match those of the corresponding new assembly. This procedure, known as liftover, allows to keep the coherence of the data stored in the repository, and can be computationally demanding. Beacuse of this, the optional parameter ```--threads``` can be used to increment the number of threads (deafult *--threads=1*). This parameter is only taken in consideration during liftover, as it will be ignored in any other situation.

### 3. Creating a remote repository
Distributed Version Control Systems (DVCS) have the characteristic of allowing a group of the users to work simultaneously in their local repositories, and eventually sharing this work with the rest of users by adding the changes into a remote repository. Being a DVCS, GenomeGit provides with the same feature to the user, who can create a remote repository at any moment by typing ```genomegit init --bare <RepositoryName>```, where ```<RepositoryName>``` stands for the name given to the bare repository created.

### 4. Remote repository access
In order to acces a remote repository, it is first needed to add a remote repository address. This can be done by typing ```genomegit remote add <remote_name> <remote_location>```. The parameter ```<remote_location>``` stands for the absolute path to the remote repository of interest if it is located within the same machine where it is being executed, or the username and server IP adress if it is located in a server (```username@xxx.xxx.xx.x```). Please note that ```<remote_name>``` does not need to have the same value as the ```<RepositoryName>``` parameter used in previous section *3. Creating a remote repository*:  ```<remote_name>``` stands simply for the nickname the user wants to refer to this repository when *pushing* and *pulling* from the local repository. To be up to date with the remote repository, the user needs to fetch the remote's data and integrate it to the local repository: ```genomegit pull <remote_name> <branch_name>```. Afterwards, user can introduce changes in the local repository and push them into the remote: ```genomegit push <remote_name>```.

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

```genomegit add --threads=8 NewAssembly.fasta```

```genomegit commit -m "First_Commit"```
#### Extract a file out of the repository
```genomegit get MyAssembly.fasta```

```genomegit get MyAnnotation.gff --sequence=Ch01 --region=1-1000000```

```genomegit get --dataset=Variants --sequence=Ch01 --region=1-1000000 --message=First_Commit```
#### Produce reports and summaries of changes between versions
```genomegit report	--commit-hash e7201f04231d039a48ea41f1f786d4f447361175```

```genomegit diff --message First_Commit Second_Commit```




 
