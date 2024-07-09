<img src="./docs/captvred_logo.png" title="CAPTVRED" alt="CAPTVRED" width="350px" align="center" />

CAPTVRED PIPELINE is designed to analyze viral NGS data. Despite, any metagenomics dataset sequenced with illumina run can be analyzed using this workflow, it has been specifically designed to asses Targeted Enrichment Sequencing (TES) datasets, with a set of probes designed based on a set of reference genomic sequences of interest. This pipeline provides an analysis for the viral identification and discovery through alignment, assembly, and taxonomic classification of the sequenced reads. It is assumed that the probes were designed based on the reference set  to enrich the samples with reads from these species and similar ones.
**In the following lines we will refer to the set of genomic reference sequences as Viral Candidates**

# Getting started:
Before running the pipeline, the file system must be prepared as follows:

### Prepare the files <br />

##### A) Viral Candidates fasta:</u> <br />
A fasta file containing the sequences of viral candidates. It  is assumed that capture probes were designed based on this set of genomic sequences, however, any set of sequences of interest will be appropriate to include. It must be a gzipped fasta. The sequence headers must contain **the identifier code followed by a space**, after the space any other information can be added if desired.<br />
<br />
 ##### B) Samples description tabular file:</u><br />
A template for this tabular file is provided (). It must be completed entering one sample per row. Additional metadata fields can be added if appropriate.<br />
<br />
 ##### C) Sequenced fastq files:</u><br />
 All sequenced fastq files must be placed (or linked) in the same directory, the IDs must correspond to the ones in the first column of the sample description files.
 
### Install Conda environment:<br />
If mamba (or Conda) software is not installed, follow the installation instructions provided in the official documentation [MISSING LINK]().
A ```environment.yaml``` file is provided to create the conda environment.

```
mamba env create -f environment.yaml
mamba activate captvred
```
<details>
 <summary>Conda environment details</summary>
 Main programs installed in conda enviroment are described here:
 
| Program                      | Version       | Channel           |
|------------------------------|---------------|-------------------|
| perl                         | latest        | defaults          |
| python                       | 3.9.2         | defaults          |
| biopython                    | latest        | conda-forge       |
| bbmap                        | latest        | bioconda          |
| fastqc                       | latest        | bioconda          |
| multiqc                      | latest        | bioconda          |
| bowtie2                      | latest        | bioconda          |
| samtools                     | latest        | bioconda          |
| seqkit                       | latest        | bioconda          |
| megahit                      | latest        | bioconda          |
| spades                       | latest        | bioconda          |
| blast                        | latest        | bioconda          |
| gawk                         | latest        | conda-forge       |
| kaiju                        | 1.9.0         | bioconda          |
| r-base                       | 4.0.5         | r                 |
| r-ggplot2                    | latest        | r                 |
| r-tidyverse                  | latest        | r                 |
| r-plyr                       | latest        | r                 |
| r-gridExtra                  | latest        | r                 |
| bioconductor-rtracklayer     | latest        | bioconda          |
| bioconductor-GenomicFeatures | latest        | bioconda          |
| bioconductor-Rsamtools       | latest        | bioconda          |
| bioconductor-GenomicAlignments | latest      | bioconda          |
| bioconductor-VariantAnnotation | latest      | bioconda          |
| bioconductor-ggbio           | latest        | bioconda          |

 </details>

### Create a project directory (*Root directory*). <br />
The pipeline and all  related files will be placed in this location.<br />
```{.sh}
mkdir MYPROJECT
cd MYPROJECT
```
### Get the pipeline.<br />
Pipeline can be downloaded via github clone repository:

   ```{.sh}
   git clone https://github.com/CompGenLabUB/CAPTVRED.git
   ```
or via nextflow pull command:
  ```{.sh}
  nextflow pull CompGenLabUB/CAPTVRED
  ```
### Set Up <br />

```{.sh}
cd MYPROJECT
nextflow conf/main.nf              \
         --set_seqs /path/to/Viral_candidates_fasta.gz  \
         --setname "MY_VIRAL_CANDIDATES"
```
> [!NOTE]
> **Some considerations:**<br />
> - Please ensure to lauch this command from the project directory.<br />
> - This step might take some time since it needs to download and process the reference database.<br />
<details>
  <summary><bl>More details about the references setup<bl></summary>
  <br />
 
  This module prepares the file setup and databases to run the pipeline afterward. In the case of databases, the main steps are:  <br />
1. __Database download__: Viral reference database [link]() most recent version is downloaded. If the database is already downloaded in the desired location it will not be downloaded again. This step can be forced by using the flag:  ```--rvdb_update```. <br />
2.  __Merge database with Viral Candidate sequences__: In this step, sequences from the viral candidates set that are not present in the reference database are included. If the merge has been done previously it will not be repeated. The merge can be forced again by using the flag ```--merge_update```.<br />
3. __Split datbase__: The full database is split into a subset containing only the sequences classified in the families of interest. Another subset containing the rest of the species is created at the same time. If the subset is already created it will not be created again. It can be forced by using the flag ```--dbsplit_update ```.<br />

Note that by redoing any of the described flags, all the downstream steps will be repeated as well. Thus, by activating the flag ```--rvdb_update```, ```--merge_update``` is automatically activated; and by activating ```--merge_update```, ```--dbsplit_update ``` is activated as well.
If the run is interrupted for any reason, remember that the ```-resume``` nextflow option will restart the pipeline from where it left off in the previous execution.

</details>

# Run CAPTVRED from the command line:

```{.sh}
nextflow $NXFDIR/main.nf
    --samp /path/to/samples_definition.tbl       \
    --fastq_dir /path/to/fastq/files/directory   \
    --runID "RUN_ID"
     -with-report $RPTDR/Nextflow_execution_report.html
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



## Data:
* The data for the test set is provided in the following [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/CAPTVRED_testset.tar.gz). The folder contains 15 test samples (3 real metagenomics samples and 12 synthetic samples), and the script used for data generation.
* The data used for the assessment of the PANDEVIR capture panel are available in the following  [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/PANDEVIR_assess_testset.tar.gz). The folder contains the raw reads for all the samples and a tabular file with the samples name relation.
