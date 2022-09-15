# SECTION 1: Prepare the environment
## 1.1. Data:
All samples must be in fastq.gz format. All of them must be placed in the 
same directory (see section 1.2). 

## 1.2. Data summary:

A tabular file containing information about the samples is necessary to 
run the pipeline (only the samples described in this tabular file will be analysed). 
This file has been written manually, it must be named "samples_definition.tbl" and has the following format:

```
#SAMPLE_ID	ILLUMINA_ID	TECHNOLOGY	PAIRED	SAMPLE_FACTOR	METHOD_FACTOR	DESCRIPTION
R01_C09_G01	G1_C_S9		Illumina	PE	Bat_guano	Capture	Sample 9 desc
R01_C10_G02	G2_C_S10	R01_C10_G02	Illumina	PE	Bat_guano	Capture	Sample 10 desc
```

Only the first two columns are used in the pipeline, but it is recomended to add some metadata information in the tabular file for further data visualisations.
Illumina ID is the identifiers given in the samples files and corresponds 
exactly to the first part of the fasta file name (i.e.: followed by R{1,2}_001.fastq.gz). Sample ID corresponds to the sample identifier that will be used in all the pipeline steps. If you are not interested in modifying the ids, use the same code in both columns.

## 1.3. Projectvars
## 2. projectvars.sh file
Includes all environment variables necessary to run the pipeline. When 
this script is executed the whole directories filesystem is created.
Every run must have its projectvars.sh inside its base folder.

* Root directory ($RDIR) is the top level directory of the project.
* Base dirctory ($BDIR) refers to the root folder where all files related 
to the current run will be stored.
* Nextflow directory ($NXFDIR) is the directory where virwaste pipeline 
is stored.
This variables are set in the _projectvars.sh_ file by the user.
Recomended filesystem structure would be:
    
    RDIR
      | - NXFDIR
      | - BDIR (Run1)
            | - projectvars.sh
      |
      | - BDIR (Run2)
            | - projectvars.sh

Afer running _projectvars.sh_:
* fastq.gz files must be placed (or linked) to the directory 
"$BDIR/rawseqs_fastq"
* samples_definition.tbl must be placed directly on the $BDIR
* reference sequences (in fa.gz) format must be placed (or linked) in the refseqs directory.


# SECTION 2: Detailed description of the pipeline steps
#2.1. QUality:

## 2.2. Cleanning (+ quality)

## 2.2. Filter and discard non viral (+ quality)

In the first step (reads cleaning) of the pipeline, the names are 
changed (thus, raw reads keep their original ID while clean reads (and 
everything that comes after) is renamed into de sample ID.
