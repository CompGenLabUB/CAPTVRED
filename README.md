<img src="./docs/captvred_logo.png" title="CAPTVRED" alt="CAPTVRED" width="350px" align="center" />

CAPTVRED PIPELINE is a pipeline designed to analyze viral NGS data obtained from targeted sequencing. More specifically, targeted sequencing with a set of probes designed based on a set of reference genomic sequences of interest. This pipeline provides an analysis for the viral identification and discovery through alignment, assembly, and taxonomic classification of the sequenced reads. It is assumed that the probes were designed based on the reference set  with the aim to enrich the samples with reads from these species and similar ones.
**In the following lines we will refer to the set of genomic reference sequences as REFSEQS**

## Run from the command line:
```{.sh}
nextflow $NXFDIR/main.nf -with-report $RPTDR/Nextflow_execution_report.html
```

### Optional parameters of the CAPTVRED pipeline:
```{.sh}
--assembler upper case string. Available options are MEGAHIT(default) and METASPADES.
--NCPUS integer. Default value is 32.
```
### Optional nextflow parameters of interest:
```{.sh}
-resume Execute the script using the cached results, useful to continue executions that were stopped by an error.
-entry  Entry workflow name to be executed.
```
All allowed commands can be found in:  _Nextflow documentation (https://www.nextflow.io/docs/latest/index.html) > Command line inteface(CLI) > Commands > run_

## Files:
In the Nextflow repository for the CAPTVRED pipeline, one can find the following files:
* projectvars.sh: &rarr; this file is not part of the nextflow environment, it is used to set the local directories before launching nextflow.
* nextflow.config : &rarr; Environment variables. This file is the same for all the runs. This file is read by default by nextflow (it must be placed in the project directory or the base directory).
* main.nf : &rarr; Controls workflow


Nextflow Modules:
* rawfq_clean.nf &rarr; BBDuk implementation for reads cleanning
* seq_stats.nf &rarr; fastQC and  MultiQC implementation. These quality control programs are run at different points in the workflow to keep track of the data quality. 
* reads_align.nf &rarr; Bowtie processes: create an index and run the alignment.
* reads_assembly.nf &rarr; megahit, trinity and spades
* contigs_align.nf &rarr; blast processes: prepare the database and run the alignment.
* contigs_taxonomy.nf &rarr; Taxonomic analysis of assembled contigs using kaiju.

## Getting started:
Before running the pipeline, the file system must be prepared as follows:

__FIRST TIME:__<br /> 
__0.__ If you don't have nextflow software: install it following the given instructions in its documentation (https://www.nextflow.io/docs/latest/getstarted.html). <br />
__a.__ Create a project directory (Root directory). The pipeline and all files related to the run will be saved inside this folder.<br />
```{.sh}
mkdir MYPROJECT
cd MYPROJECT
```
__b.1.__ Clone the repository of the pipeline (folder can be named: CAPTVRED). <br />
```{.sh}
git clone https://github.com/JosepFAbril/CAPTVRED.git
```
__b.2.__  Inside the root directory create a directory with the run name ID (Run directory or base directory). <br />
```{.sh}
cd MYPROJECT
mkdir DATASET_00
```
__c.__ Prepare reference sequences (fasta file with all reference sequences and _gff_ files for all genomes) and kaiju databases (see [documentation](https://github.com/JosepFAbril/CAPTVRED/blob/main/docs/readme_DOCUMENTATION_virwaste.md) for more info). <br />

__FOR EACH RUN or EXPERIMENT:__<br /> 

__d.1.__ Place samples_definition.tbl file in the run directory. It must be a tabular file containing the sequencing ID and the processing ID for each sample (see [documentation](https://github.com/JosepFAbril/CAPTVRED/blob/main/docs/readme_DOCUMENTATION_virwaste.md) for more info). Additional metadata fields can be added if you want.<br />
__d.2.__ Place projectvars template in the run/experiment folder and define variable names: <br />

```{.sh}
cp CAPTVRED/projectvars_template.sh  DATASET_01/projectvars.sh
```
Variables description
-  $RDIR (path to the root directory)<br />
-  $NXFDIR (path to the CAPTVRED pipeline directory) <br /> 
-  $BDIR (path to the run directory or base directory) <br />
-  $RUNID (run/experiment identifier)<br />
-  $AMPSQD (reference sequences fasta directory)<br />
-  $AMPSQFA (reference sequences fasta filename)<br />
-  $KAIDB (kaiju databases directory)<br />
-  $AMPGFFD (reference sequences _gff_)<br />
It is recommended to rename this file as "projectvars.sh".<br />

__e.__ Run projectvars.sh. (see [documentation](https://github.com/JosepFAbril/CAPTVRED/blob/main/docs/readme_DOCUMENTATION_virwaste.md) for more info)<br />
```{.sh}
cd DATASET_01
source projectvars.sh
```
__f.__ Place (or link) the raw fastq files to the "rawseqs_fastq" directory.<br />
__g.__ Filesystem is ready to run the pipeline.<br />See parameters.
```{.sh}
nextflow $NXFDIR/main.nf -with-report $RPTDR/Nextflow_execution_report.html
```

## Data:
* The data for the test set is provided in the following [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/CAPTVRED_testset.tar.gz). The folder contains 15 test samples (3 real metagenomics samples and 12 synthetic samples), and the script used for data generation.
* The data used for the assessment of the PANDEVIR capture panel are available in the following  [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/PANDEVIR_assess_testset.tar.gz). The folder contains the raw reads for all the samples and a tabular file with the samples name relation.


