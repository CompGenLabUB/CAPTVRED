# virwaste
Analysis pipeline for the VIRWASTE viral metagenomics project.

## Run from cmdline:
```{.sh}
source projectvars.sh
nextflow $NXFDIR/main.nf -with-report $RPTDR/Nextflow_execution_report.html
```

## Files:
In the Nextflow repository for the VIRWASTE pipeline one can find the following files:
* projectvars.sh -> this file is not part of the nextflow environment, is used to set the local directories before launchong nextflow.
* nextflow.config -> Environment variables. This file is the same for all the runs. This file is read by default by nextflow (it must be placed in the project directory or in the base directory).
* main.nf -> Controls workflow
* rawfq_clean.nf -> BBDuk implementation for reads cleanning
* seq_stats.nf -> fastQC and  MultiQC implementation. This quality control programms are run at different points in the wotkflow to keep track of the quality of the data. 

## Getting started:
Before running the pipeline, file system must be prepared as follows:

Virwaste:

 - nextflow.config 
 - main.nf\newline
 - rawfq_clean.nf
 - check_inputs.nf
   
RUN_XX:

   - projectvars.sh  
      modified with specific run parameters
   - samples_definition.tbl
      Tabular file with the following fields:  
      "SAMPLE_ID	ILLUMINA_ID	TECHNOLOGY	PAIRED	SAMPLE_FACTOR	METHOD_FACTOR	DESCRIPTION"
   - rawseqs_fastq 
      Directory containing ".fastq.gz" files for all paired end samples
   - refseqs
      Directory containing a ".fasta.gz" of reference sequences used to design the library 
      You may want to use a link to reference to the data directory in other location
      Bowtie index for the aligning step will be created in the directory where the fasta file is located
  

Once $projectvars.sh$ is run, "reports", "cleanseqs" and  "amplicons_alignment" directories are automatically created.
