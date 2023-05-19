<img src="./docs/captvred_logo.png" title="CAPTVRED" alt="CAPTVRED" width="350px" align="center" />

CAPTVRED PIPELINE is a pipeline designed to analyse viral NGS data obtained from targeted sequencing. More specifically, targeted sequencing with a set of probes designed based on a set of reference genomic sequences of interest. This pipeline provides a set of analysis for the viral identification and discovery through alignment, assembly and taxoniomic classification of the data. It is assumed that the proves were designed basedon the reference set  with the aim to find this sequences and other alike sequences.

## Run from cmdline:
```{.sh}
nextflow $NXFDIR/main.nf -with-report $RPTDR/Nextflow_execution_report.html
```

### Optional parameters of the CAPTVRED pipeline:
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
In the Nextflow repository for the CAPTVRED pipeline one can find the following files:
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
__0.__ If you don't have nextflow fostware install it following the given instructions in its documentation (https://www.nextflow.io/docs/latest/getstarted.html). <br />
__a.__ Create a project directory (Root directory). The pipeline and all run will be saved inside this folder.<br />
__b.1.__ Inside the root directory download the repository of the pipeline (folder can be named: CAPTVRED). <br />
__b.2.__  Inside the root directory create a directory with the run name ID (Run directory or base directory). <br />
__c.__ Prepare reference sequences (fasta file with all reference sequences and gff files for all genomes) and kaiju databases (see documentation for more info). <br />

__FOR EACH RUN or EXPERIMENT:__<br /> 

__d.1.__ Place samples_definition.tbl file in the run direcoty. It must be a tabular file containing sequencing aID and processing ID for each sample (see documentation for more info). Additional metadata fields can be added if you want.<br />
__d.2.__ Place projectvars template in the run folder and add the following variables: <br />
-  $RDIR (path to the root directory)<br />
-  $NXFDIR (path to the nextflow directory) <br /> 
-  $BDIR (path to the run dierctory or base directory) <br />
-  $RUNID (run identifier)<br />
-  $AMPSQD (reference sequences fasta directory)<br />
-  $AMPSQFA (reference sequences fasta filename)<br />
-  $KAIDB (kaiju databases directory)<br />
-  $AMPGFFD (reference sequences gff)<br />
It is recomended to rename this file as "projectvars.sh".<br />

__e.__ Run projectvars.sh. (see documentation for more info)<br />
```{.sh}
source projectvars.sh
```
__f.__ Move or link) the raw fastq files to the "$RAWFQ" directory.<br />
__g.__ Check the configuration file ($NXFDIR/nextflow.config) and, in case you are interested, change the parameters.<br />
__h.__ Filesystem is ready to run the pipeline.<br />
```{.sh}
nextflow $NXFDIR/main.nf -with-report $RPTDR/Nextflow_execution_report.html
```


### To take into consideration:
If assembly is performed using trinity and the prograam cannot assembly any contig the pipleni will ignore the programme fail, raise a warning and continue with the next file.
