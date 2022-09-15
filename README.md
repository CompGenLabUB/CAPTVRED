# virwaste
Analysis pipeline for the VIRWASTE viral metagenomics project.

## Run from cmdline:
```{.sh}
nextflow $NXFDIR/main.nf -with-report $RPTDR/Nextflow_execution_report.html
```

### Optional parameters of the virwaste pipeline:
```{.sh}
--assembler upper case string. Available options are MEGAHIT(default), TRINITY (time and resource consuming) and SPADES(?).
--NCPUS integer. Default is 32.
```
### Optional nextflow parameters of interest:
```{.sh}
-resume Execute the script using the cached results, useful to continue executions that was stopped by an error.
-entry  Entry workflow name to be executed.
```
All allowed commands can be found in:  _Nextflow documentation (https://www.nextflow.io/docs/latest/index.html) > Command line inteface(CLI) > Commands > run_

## Files:
In the Nextflow repository for the VIRWASTE pipeline one can find the following files:
* projectvars.sh -> this file is not part of the nextflow environment, is used to set the local directories before launchong nextflow.
* nextflow.config -> Environment variables. This file is the same for all the runs. This file is read by default by nextflow (it must be placed in the project directory or in the base directory).
* main.nf -> Controls workflow


Nextflow Modules:
* rawfq_clean.nf -> BBDuk implementation for reads cleanning
* seq_stats.nf -> fastQC and  MultiQC implementation. This quality control programms are run at different points in the wotkflow to keep track of the quality of the data. 
* reads_align.nf -> Bowtie processes: create index and run the alignment ().
* reads_assembly.nf -> megahit, trinity and spades
* contigs_align.nf -> blast processes: create database and run the alignment.
* contigs_taxonomy.nf -> Taxonomic analysis of assembled contigs using kaiju.

## Getting started:
Before running the pipeline, file system must be prepared as follows:

__FIRST TIME:__<br /> 

__a.__ Create a project directory (Root directory). The pipeline and all run will be saved inside this folder.<br />
__b.1.__ Inside the root directory download/install the pipeline (i.e. folder can be named: virwaste).  ... description of how to do that.<br />
__b.2.__  Inside the root directory create a directory with the run name ID (Run directory or base directory).<br />

__FOR EACH RUN or EXPERIMENT:__<br /> 

__c.1.__ Place samples_definition.tbl file in the run direcoty. It must be a tabular file containing sequencing aID and processing ID for each sample (see documentation for more info). Additional metadata fields can be added if you want.<br />
__c.2.__ Place projectvars template in the run folder and add the variables $RDIR (root directory), $NXFDIR (nextflow directory), $BDIR (Run dierctory or base directory) and $RUNID (run identifier). It is recomended to rename this file as "projectvars.sh".<br />
__d.__ Run projectvars.sh. (see documentation for more info)<br />
```{.sh}
source projectvars.sh
```
__e.__ Move or link) the raw fastq files to the "$RAWFQ" directory.<br />
__f.__ Check the configuration file ($NXFDIR/nextflow.config) and, in case you are interested, change the parameters.<br />
__g.__ Filesystem is ready to run the pipeline.<br />
```{.sh}
nextflow $NXFDIR/main.nf -with-report $RPTDR/Nextflow_execution_report.html
```


### To take into consideration:
If assembly is performed using trinity and the prograam cannot assembly any contig the pipleni will ignore the programme fail, raise a warning and continue wwith the next file.
