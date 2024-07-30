<img src="./docs/captvred_logo.png" title="CAPTVRED" alt="CAPTVRED" width="350px" align="center" />

CAPTVRED PIPELINE is designed to analyze viral metagenomics datasets from Target Enrichment Sequencing (or Capture-based Metagenomics). This pipeline provides an analysis for viral identification through alignment, assembly, and taxonomic classification of the sequenced reads. The analyses focus on the set of species of interest, for which the dataset has been enriched, and other related sequences from the same taxonomic family. **In the following lines we will refer to the set of genomic sequences of interest as Viral Candidates.**

# Getting started:
Before running the pipeline, the file system must be prepared as follows:

### Prepare the files <br />

##### A) Viral Candidates fasta:</u> <br />
A fasta file containing the sequences of viral candidates. It  is assumed that capture probes were designed based on this set of genomic sequences, however, any set of sequences of interest will be appropriate to include. It must be a gzipped fasta. The sequence headers must contain **the identifier code followed by a space**, after the space any other information can be added if desired.<br />
<br />

##### B) Samples description tabular file:</u><br />
A template for this tabular file is provided (samples_definition_template.sh). It must be completed entering one sample per row. Additional metadata fields can be added if appropriate.<br />

<details>
 <summary>Samples description details</summary>
 This tabular file has two required fields:
 
 * Sample ID: This is the identifier that the pipeline will use to name all files and in the final report.
   
 * IlluminaID: This is the file name prefix containing raw data. Suffix of the samples indicating R1 and R2 files can be modified in the nextflow.config file or in the commandline using ```--R1``` [default: R1_001] and ```--R2``` [default: R2_001].
 The rest of the fields proposed in the template are recommendations and will probably be used in future versions of the workflow.
 </details>
 
 ##### C) Sequenced fastq files:</u><br />
 All sequenced fastq files must be placed (or linked) in the same directory, the IDs must correspond to the ones in the first column of the sample description files.
 
### Get the pipeline<br />
Pipeline can be downloaded via github clone repository:

   ```{.sh}
   git clone https://github.com/CompGenLabUB/CAPTVRED.git
   ```
or via nextflow pull command:
  ```{.sh}
  nextflow pull CompGenLabUB/CAPTVRED
  ```

### Install Conda environment:<br />
If mamba (or Conda) software is not installed, follow the installation instructions provided in the [documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
A ```environment.yml``` file is provided in this repository to create the conda environment.

```
cd CAPTVRED
mamba env create -f environment.yaml
mamba activate captvred
```
<details>
 <summary>Conda environment details</summary>
 The main programs installed in conda environment are described here:
 
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

### Create a project directory <br />
The pipeline and all  related files will be placed in this location.<br />
```{.sh}
mkdir -vp MYPROJECT
cd MYPROJECT
```

### Set Up <br />

```{.sh}
cd CAPTVRED/conf
nextflow main.nf              \
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
1. __Database download__: Viral reference database ([RVDB](https://rvdb.dbi.udel.edu/)) most recent version is downloaded. If the database is already downloaded in the desired location it will not be downloaded again. This step can be forced by using the flag:  ```--db_update```. <br />
2.  __Merge database with Viral Candidate sequences__: In this step, sequences from the viral candidates set that are not present in the reference database are included. If the merge has been done previously it will not be repeated. The merge can be forced again by using the flag ```--merge_update```.<br />
3. __Split datbase__: The full database is split into a subset containing only the sequences classified in the families of interest. Another subset containing the rest of the species is created at the same time. If the subset is already created it will not be created again. It can be forced by using the flag ```--dbsplit_update ```.<br />

Note that by redoing any of the described flags, all the downstream steps will be repeated as well. Thus, by activating the flag ```--db_update```, ```--merge_update``` is automatically activated; and by activating ```--merge_update```, ```--dbsplit_update ``` is activated as well.
If the run is interrupted for any reason, remember that the ```-resume``` nextflow option will restart the pipeline from where it left off in the previous execution.

</details>

# Run CAPTVRED:

```{.sh}
cd CAPTVRED
nextflow main.nf
          --samp /path/to/samples_definition.tbl       \
          --fastq_dir /path/to/fastq/files/directory   \
          --runID "RUN_ID"
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

# Output:
CAPTVRED produces an HTML report summarizing key findings to facilitate the visualization and interpretation of the results. From this page, the user can access the quality and the computational performance reports. It also includes summary tables for metadata and sequence recovery. 
For each sample three tables are generated for viral assignations:

* Read level:
  
| READ_ID   | TAG | READ_LENGTH | BESTHIT_LEN | BESTHIT_COVERAGE | REFSEQ_ID   | KAIJU_SCORE | NREADS_MAPED | TAXONID   | SPECIES                                   | FAMILY            |
|-----------|-----|-------------|-------------|------------------|-------------|-------------|--------------|-----------|-------------------------------------------|-------------------|
| k127_80   | B   | 2189        | 2189        | 100.000          | NC_002023.1 | NA          | 3178         | 11320     | Influenza_A_virus                         | Orthomyxoviridae  |
| k127_67   | B   | 5731        | 5731        | 100.000          | AF208067.1  | NA          | 10819        | 694005    | Murine_coronavirus                        | Coronaviridae     |
| k127_95   | B   | 10571       | 10571       | 100.000          | NC_001474.2 | NA          | 19942        | 12637     | Dengue_virus                              | Flaviviridae      |
| k127_19   | B   | 1608        | 1608        | 100.000          | FJ390061.2  | NA          | 2016         | 11320     | Influenza_A_virus                         | Orthomyxoviridae  |
| k127_23   | B   | 1983        | 1983        | 100.000          | KX377335.1  | NA          | 3323         | 64320     | ZIKV                                      | Flaviviridae      |
| k127_111  | B   | 2081        | 2081        | 100.000          | KJ633807.1  | NA          | 2962         | 11320     | Influenza_A_virus                         | Orthomyxoviridae  |
| k127_74   | B   | 425         | 425         | 100.000          | KJ633811.1  | NA          | 275          | 11320     

<details>
 <summary>Fields description:</summary>
 * **READ_ID**: Uniq identifier for each read/contig.
 * **TAG**: B for blastn, T for tblastx and K for kaiju.
 * **READ_LENGTH**: Read or contig length in bp.
 * **BESTHIT_LEN**: Length of the best hit.
 * **BESTHIT_COVERAGE**: Coverage of the best hit.
 * **REFSEQ_ID**: Assignation specie sequence id.
 * **KAIJU_SCORE**: Score reported by kaiju NA if blastn (default) option is running.
 * **NREADS_MAPED**: Number of raw reads mapped to this seqid. 
 * **TAXONID**: Sequence taxon id.
 * **SPECIES**: Species name.
 * **FAMILY**: Species family taxonomic classification.
</details>
  
* Sequence level:

| SEQUENCE_ID | TAXON_ID | SPECIES              | FAMILY            | SEQ_LENGTH | NUCS_ALN | COVERAGE_PCT | BHIT_IDENTITY | BESTHSP_COUNT | NREADS_MAPED |
|-------------|----------|----------------------|-------------------|------------|----------|--------------|---------------|---------------|--------------|
| NC_001474.2 | 12637    | Dengue_virus         | Flaviviridae      | 10723      | 10571    | 98.582       | 100.000       | 1             | 19942        |
| KJ633811.1  | 11320    | Influenza_A_virus    | Orthomyxoviridae  | 1027       | 850      | 82.765       | 100.000       | 2             | 550          |
| NC_007357.1 | 11320    | Influenza_A_virus    | Orthomyxoviridae  | 2341       | 2189     | 93.507       | 100.000       | 1             | 3178         |
| DQ415901.1  | 290028   | Human_CoV/HKU1       | Coronaviridae     | 30097      | 29945    | 99.495       | 99.588        | 2             | 58330        |
| OK017853.1  | 2833184  | Sarbecovirus_sp.     | Coronaviridae     | 29369      | 29369    | 100.000      | 95.796        | 1             | 57998        |
| AF255742.1  | 11320    | Influenza_A_virus    | Orthomyxoviridae  | 1565       | 1411     | 90.160       | 99.929        | 1             | 1610         |

* Species level:

| TAXON_ID | SPECIES                                              | FAMILY         | MIN_COV | MAX_COV | MEAN_COV | MIN_PID | MAX_PID | MEAN_PID | N_SEQS | N_CONTIGS | NREADS_MAPED | INFO                                                                                           |
|----------|------------------------------------------------------|----------------|---------|---------|----------|---------|---------|----------|--------|-----------|--------------|------------------------------------------------------------------------------------------------|
| 1335626  | Middle_East_respiratory_syndrome-related_coronavirus | Coronaviridae  | 99.393  | 99.393  | 99.393   | 99.349  | 99.349  | 99.349   | 1      | 1         | 58734        | SEQ:MK129253.1,LEN:30150,NUCALN:29967,COV:99.393,BHIDENT:99.349,N:1                            |
| 694009   | Severe_acute_respiratory_syndrome-related_coronavirus| Coronaviridae  | 99.492  | 99.492  | 99.492   | 99.946  | 99.946  | 99.946   | 1      | 1         | 116308       | SEQ:NC_045512.2,LEN:29903,NUCALN:29751,COV:99.492,BHIDENT:99.946,N:1                           |
| 694000   | Miniopterus_bat_coronavirus_1                        | Coronaviridae  | 99.463  | 99.463  | 99.463   | 100.000 | 100.000 | 100      | 1      | 1         | 55148        | SEQ:NC_010437.1,LEN:28326,NUCALN:28174,COV:99.463,BHIDENT:100.000,N:1                          |
| 37124    | Chikungunya_virus                                     | Togaviridae    | 97.966  | 97.966  | 97.966   | 99.734  | 99.734  | 99.734   | 1      | 1         | 22148        | SEQ:MG280943.1,LEN:11896,NUCALN:11654,COV:97.966,BHIDENT:99.734,N:1                            |
| 694014   | Avian_coronavirus                                    | Coronaviridae  | 99.352  | 99.352  | 99.352   | 99.993  | 99.993  | 99.993   | 1      | 1         | 53712        | SEQ:AJ311317.1,LEN:27635,NUCALN:27456,COV:99.352,BHIDENT:99.993,N:1                            |
| 2496529  | Mengla_virus                                         | Filoviridae    | 99.169  | 99.169  | 99.169   | 100.000 | 100.000 | 100      | 1      | 2         | 35000        | SEQ:NC_055510.1,LEN:18300,NUCALN:18148,COV:99.169,BHIDENT:100.000,N:2                          |


# Data:
* The data for the test set is provided in the following [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/CAPTVRED_testset.tar.gz). The folder contains 15 test samples (3 real metagenomics samples and 12 synthetic samples), and the script used for data generation.
* The data used for the assessment of the PANDEVIR capture panel are available in the following  [link](https://compgen.bio.ub.edu/datasets/CAPTVRED/PANDEVIR_assess_testset.tar.gz). The folder contains the raw reads for all the samples and a tabular file with the samples name relation.
