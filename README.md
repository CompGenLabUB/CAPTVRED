<img src="./docs/captvred_logo.png" title="CAPTVRED" alt="CAPTVRED" width="350px" align="center" />

CAPTVRED PIPELINE is designed to analyze viral NGS data. Any metagenomics dataset sequenced with illumina run can be analyzed using this workflow, however it has been specifically designed to asses targeted sequencing with a set of probes designed based on a set of reference genomic sequences of interest. This pipeline provides an analysis for the viral identification and discovery through alignment, assembly, and taxonomic classification of the sequenced reads. It is assumed that the probes were designed based on the reference set  with the aim to enrich the samples with reads from these species and similar ones.
**In the following lines we will refer to the set of genomic reference sequences as Viral Candidates**

# Getting started:
Before running the pipeline, the file system must be prepared as follows:


### __0.__ Install Nextflow software:<br />
If you don't have nextflow software: install it following the given instructions in its documentation (https://www.nextflow.io/docs/latest/getstarted.html). <br />

### __1.__ Create a project directory (*Root directory*). <br />
The pipeline and all files related to the run will be saved inside this folder.<br />
```{.sh}
mkdir MYPROJECT
```
### __2.__ Clone the repository of the pipeline.<br />
Folder can be named: CAPTVRED. <br />
```{.sh}
cd MYPROJECT
git clone https://github.com/CompGenLabUB/CAPTVRED.git
```

### __3.__ Prepare reference sequences. <br />


 #### __3.a.__ Preset reference sequences: <br /> 
Preset reference sequences are provided as _tar.gz_ file in the following [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/REFSEQS.tar.gz). It includes RVDB_NT database for the contigs taxonomy, nr_euk database for fast preliminar classification and illumina adapters fasta file. To use this set of preset sequences download and extract: <br />

```{.sh}
cd MYPROJECT
wget https://compgen.bio.ub.edu/datasets/CAPTVRED/CAPTVRED_refseqs.tar.gz
tar –xvzf REFSEQS.tar.gz
```
If you plan to use the capture-based approach: <br />
> __I.__ Move (or link) fasta file of viral candidates in the _"REFSEQS/viral_candidates"_ directory. <br />
>```{.sh}
>mv  /dir/to/viral_candidates.fa.gz  ./REFSEQS/viral_candidates/Viral_candidates.fa.gz
>```
>__II.__ Move (or link) information file of viral candidates in the _"REFSEQS/viral_candidates"_ directory.<br />    
>```{.sh}
>mv  /dir/to/viral_candidates_info.tsv  ./REFSEQS/viral_candidates_info.tsv 
>```
>Info file has the following format (tab separated): <br />
>
>
>|**##Famlily** |**Specie**      |**SpecTaxonId**|**TaxonId**|**Host**  |**SeqID**  |**Region**    |**Size**	   |**Name**                          |**decription**  |
>|--------------|----------------|--------------|---------|------|--------------|-----------------|---------|--------------------------------------|----------------|
>|Filoviridae	 |Cuevavirus	    |1513237	     |1513237	 |Bat	  |NC_016144.1	 |CompleteGenome	 |18927	   |Cuevavirus_Lloviu_cuevavirus_isolate  |Lloviu virus/M.schreibersii-wt/ESP/2003/Asturias-Bat86, complete genome |
>|Bunyavididae	|Crimean-congo_virus |1980519	   |1980519	 |Human	|NC_005301.3	 |Lsegment	       |12108	   |Crimean-congo_virus_Lseg	            |Crimean-Congo hemorrhagic fever >virus segment L, complete sequence |
> 
>__III.__ Move (or link) gff files of all viral candidates into _"REFSEQS/viral_candidates/gff_genomes"_ directory. <br />
>
>```{.sh}
>mv  /dir/to/gff_genomes_dir/*.gff   ./REFSEQS/gff_genomes 
>```
<br />

#### __3.b.__ Customized reference sequences: <br /> 
More information if provided in the [documentation file](https://github.com/JosepFAbril/CAPTVRED/blob/main/docs/readme_DOCUMENTATION_virwaste.md) 


# For each new run or experiment:


### __4__ Prepare the filesystem for each run or experiment:<br /> 


#### __4.1.__ Create directory for the new run:
```{.sh}
mkdir MYPROJECT/RUN_00
```



#### __4.2.__ Fill samples definition file:
Place samples_definition.tbl template in the run directory and fill the required information.
```{.sh}
cp MYPROJECT/CAPTVRED/samples_definition_template.tbl  MYPROJECT/RUN_00/samples_definition.tbl
```
Samples definition file must be completed entering one sample per row. Additional metadata fields can be added if appropiate.<br />
<br />


#### __4.3.__ Set project variables and parameters: 
Place projectvars template in the run/experiment folder and define variable names: <br />

```{.sh}
cp MYPROJECT/CAPTVRED/projectvars_template.sh  MYPROJECT/RUN_00/projectvars.sh
```
Projectvars file is divided in three main blocks: *Run information*, *Options and variables*, and *Paths and filesystems*. Most of the parameters are set to the default values and can be optionally modified, however, there are some variables that have to be set by the user (marked with a star ( :star: ))

* Run information
  
In this block, ***RUNID*** ( :star: ) must be filled with the name of the run/exepriment. i.e.RUN_00.
***R1*** and ***R2*** ( :star: ) samples suffixes must be set too, it refers to the suffix in the fastq filenames to identify R1 and R2 files. Other variables must stay fixed, please do not modify them.

* Options and variables
  
Despite all parameters can be modified directly from the commandline, some of them can be changed from the project vars as well. Default options are already set. 
  
* Paths and filesystem
The analyses filesystem is described in this block. 
The most important ones are:
  -  ***$RDIR*** ( :star: )- Path to the root directory (i.e./data/metagenomics/MYPROJECT). <br />
  -  ***$NXFDIR*** - path to the CAPTVRED pipeline directory, by default set to *$RDIR/CAPTVRED* <br /> 
  -  ***$BDIR*** - path to the run directory or base directory, by default set to *${RDIR}/${RUNID}* <br />
  - ***$REFSQD*** - path to the reference sequences directory. 
  -  ***$AMPSQD***  viral candidates seqs directory, ***$AMSQFA***  ( :star: ) (viral candidates reference sequences in fasta format), ***$AMPSQINFO***  ( :star: ) (information tabular file of viral candidates), and ***AMPGFFD***   (direcory of viral candidates genomes gff files) can be modified with the name of the user files if necessary. <br />
  
For running the pipeline with the provided reference sequences no other variables need to be modified. Blast and kaiju databases or contamination filtering references can be modified in this section as well. for nore information see [documentation](https://github.com/JosepFAbril/CAPTVRED/blob/main/docs/readme_DOCUMENTATION_virwaste.md)

> [!IMPORTANT]  
> The pipeline have some package and program dependencies, make sure the location of all this dependencies are set in the variable "$PATH", otherwise, CAPTVRED will not find the and the workflow wit crash. To make sure it does not happen, add the required paths to the "PATH" variable in this section (see an example in the template).
 
__e.__ Run projectvars.sh. (see [documentation](https://github.com/JosepFAbril/CAPTVRED/blob/main/docs/readme_DOCUMENTATION_virwaste.md) for more info)<br />
```{.sh}
cd MYPROJECT/RUN_00
source projectvars.sh
```

__f.__ Place (or link) the raw fastq files to the "rawseqs_fastq" directory.<br />
```{.sh}
mv dir/to/my/fastq/files/*.fq.gz  MYPROJECT/RUN_00/rawseqs_fastq
```

__g.__ Filesystem is ready to run the pipeline.<br />See parameters.
```{.sh}
nextflow $NXFDIR/main.nf -with-report $RPTDR/Nextflow_execution_report.html
```




# Run CAPTVRED from the command line:
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



## Data:
* The data for the test set is provided in the following [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/CAPTVRED_testset.tar.gz). The folder contains 15 test samples (3 real metagenomics samples and 12 synthetic samples), and the script used for data generation.
* The data used for the assessment of the PANDEVIR capture panel are available in the following  [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/PANDEVIR_assess_testset.tar.gz). The folder contains the raw reads for all the samples and a tabular file with the samples name relation.
