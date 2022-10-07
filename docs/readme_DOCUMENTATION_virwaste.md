<style>
body { counter-reset: h1counter h2counter h3counter h4counter h5counter h6counter; }

h1 { counter-reset: h2counter; }
h2 { counter-reset: h3counter; }
h3 { counter-reset: h4counter; }
h4 {}

pre {
  background-color: #CCCCCC;
}

h1:before {
    counter-increment: h1counter;
    content: counter(h1counter) ".\0000a0\0000a0";
}

h2:before {
    counter-increment: h2counter;
    content: counter(h1counter) "." counter(h2counter) ".\0000a0\0000a0";
}

h3:before {
    counter-increment: h3counter;
    content: counter(h1counter) "." counter(h2counter) "." counter(h3counter) ".\0000a0\0000a0";
}

h4:before {
    counter-increment: h4counter;
    content: counter(h1counter) "." counter(h2counter) "." counter(h3counter) "." counter(h4counter) ".\0000a0\0000a0";
}


</style>

# Introduction

# Run the pipeline

## Use the github link

In this case you have to download the *projectvars.sh* file and, if necessary, 
the configuration file ( *nextflow.config* ). 
En omplir el projectvars pots deixar la variable $NXFDIR empty. __!! REVISAR !!__

```{.sh}

```

## Download and install manually the pipeline

Descarregar tot el repository de github.

```{.sh}

```

## Rerun

# Prepare the environment


## Data summary:

A tabular file containing information about the samples is necessary to 
run the pipeline (only the samples described in this tabular file will 
be analysed). This file has been written manually, it must be named 
"samples_definition.tbl" and has the following format:

```
#SAMPLE_ID	ILLUMINA_ID	TECHNOLOGY	PAIRED	SAMPLE_FACTOR	METHOD_FACTOR	DESCRIPTION
R01_C09_G01	G1_C_S9	Illumina	PE	Bat_guano	Capture	Sample 9 desc
R01_C10_G02	G2_C_S10	R01_C10_G02	Illumina	PE	Bat_guano	Capture	Sample 10 desc
```

Only the first two columns are used in the pipeline, but it is recomended 
to add some metadata information in the tabular file for further data visualisations.
Illumina ID is the identifiers given in the samples files and corresponds 
exactly to the first part of the fasta file name (i.e.: followed by R{1,2}_001.fastq.gz). 
Sample ID corresponds to the sample identifier that will be used in all the pipeline steps.
(In the first step of the pipeline (reads cleaning), the names are 
changed in such a way that raw reads keep their original ID while clean reads (and 
everything that comes after) is renamed into de sample ID.) 
If you are not interested in modifying the ids, use the same code in both columns.

## Input data:

All samples must be in fastq.gz format. All of them must be placed in the 
same directory (see section 1.2). 

It is assumed that input sequence files are paired files named following
a pattern. The root name for each sample must be described in the [data summary](#data-summary),
the pattern for the suffix refering to each pair of reads must be described
in the [ *projectvars.sh* file](#projectvarssh-file).

i.e. if you set of samples is:

```
sample_AA_april_2021_r01.fastq.gz
sample_AA_april_2021_r02.fastq.gz
sample_bb_april_2021_r01.fastq.gz
sample_bb_april_2021_r02.fastq.gz
sample_C_april2021_r01.fastq.gz
sample_C_april2021_r02.fastq.gz
```

Your data summary could be similar to:

```
sample_A    sample_AA_april_2021
sample_B    sample_bb_april_2021
sample_C    sample_C_april2021
```
In this manner you simplify and standarize the name of the files generated 
all along the pipeline. 

And, finally, you sould define pattern in the projectvars as `R1="_r01"`, since this
is the pattern for the fist file of each pair, and `R1="_r02"`, since this
is the pattern for the second file of each pair.


## Reference sequences

__Tabular info file__  __!!! REVISAR !!!__<br />
Descriptive tabular file including the NCBI id, full name, short name and 
taxon ID for each of the reference sequences in the reference set.  

camps= ##Famlily	Specie	Host	SeqID	Region	Size	Name	decription
```
```

__Fasta file__  <br />
Fasta file containing the reference sequences. Identifier must be the NCBI ID 
and optionally other information (separated with a space!). 
This sequences will be used as reference database for alignments (performed 
with bowtie for the reads and blast for the contigs) and for the taxonomic 
classification (see "Kaiju databases"  subsection).


__GFF files__  __!! REVISAR !!__<br />
A separated GFF file 
This files can be obtained by getting the genebank information of the genmome, 
extracting the gff3 and selecting the coding regions.

This is the command we used to fetch our gff sequences:
```{.sh}
gawk '{ print $1 }' Viral_candidates_zoonosis_refseqs_idrelation.tsv | \
  while read n;
    do {
      echo "# $n" 1>&2;
      efetch -db nuccore -format gb -id "$n" > "$n".gb;
      /usr/bin/bp_genbank2gff3 refseqs/"$n".gb -o stdout | \
        gawk -v pfx="$n" -v isgff=1 '
             $0 ~ /^##FASTA/ { isgff=0; }
             { ofl= isgff ? ".gff" : ".gff.fasta";
               print $0 > pfx""ofl;
               }' -;
        gawk '$3 ~ /region|CDS/' $n.gff \      ## potser no cal
                               > $n.feats.gff;
    }; done
```
Where _Viral_candidates_zoonosis_refseqs_idrelation.tsv_ is the tabular
file mentioned in previous subsection.

__Kaiju Databases__ <br />

Taxonomic classification with kaiju is performed with multiple kaiju 
databases. By default the virwaste pipeline runs the taxonomic analysis 
on 5 diferent databases: nr_euk, refseq, viruses, rvdb and a database 
created with the set of amplicon reference sequences (named capref).

In the Kaiju documentation, the instructions to download the availabe
source databases and to create a new database from a fasta file 
(https://github.com/bioinformatics-centre/kaiju). 

## Projectvars.sh file

Includes all environment variables necessary to run the pipeline. When 
this script is executed the whole directories filesystem is created.
Every run must have its projectvars.sh inside its base folder.

* Run identifier ($RUNID) is the especific name of the sequencing experiment.  <br />
* $R1 and $R2, which corresponds to the pattern used to describe fist and decond pair of
files respectively. (see [Input Data section](#input-data))

* Root directory ($RDIR) is the top level directory of the project.
* Base dirctory ($BDIR) refers to the root folder where all files created 
in each of the steps in the pipeline will be stored.
* Nextflow directory ($NXFDIR) is the directory where virwaste pipeline 
is stored.
* Directory that contains the fasta file with the reference sequences and 
the tabular info file ($AMPSQD). The pipeline will create  bowtie and
blast databases and indexes of the fasta in this folder. <br />
* Name of the file contianing the amplicon reference sequences ($AMPSQFA).
Bowtie and blast databases will be created with the same prefix name. <br />
* Directory where the kaiju databases are stored ($KAIDB) each database 
must be stored in its own directory. __!! REVISAR !!__<br />

This variables are set in the _projectvars.sh_ file by the user.
Recomended filesystem structure would be:

```{=latex}
    RDIR<br />
      | - NXFDIR<br />
      | - BDIR (Run1)<br />
            | - projectvars.sh<br />
      |<br />
      | - BDIR (Run2)<br />
            | - projectvars.sh<br />
```
 
After running _projectvars.sh_:
* fastq.gz files must be placed (or linked) to the directory 
"$BDIR/rawseqs_fastq"
* samples_definition.tbl must be placed directly on the $BDIR
* reference sequences (in fa.gz) format must be placed (or linked) in the refseqs directory.

## Configuration file:

The configuration file is named *nextflow.config*, you can find it in the workflow main directory. 
It contains the default variables used in each of the steps of the pipeline. 
If you are interested in changing any of the parameters befor running the pipeline 
you can:

> **a) Modify the document**
All parameters can be modified directly on the configuration file. 
Preferably, do not modify the filename keep it on the same directory. 
If you change the filename or directory of the configuration file, or if
you are running the pipeline using the github link, you can provide the name 
of the configuration file as a commandline argument using the `-params-file` option.


> **b) Add the parameter as a commandline option**

If you are interested in adjusting only one or few parameters it might be more
convinient to provide this options directly as command line options. You can do so 
by simply writting: `--ParameterOfInterest NewValue`.
i.e. If you are interested in changing the default number of CPUS to 16 you
can run the pipeline as follows:

```{.sh}
nextflow run https://github.com/JosepFAbril/virwaste --NCPUS 16
```

You can find more information on how to custom and provide the configuration file in the
 [nextflow config documentation](https://www.nextflow.io/docs/latest/config.html) .

# Detailed description of the pipeline steps

## Quality:

Reads quality is performed at four different stages during the pipeline, 
(1) in the raw reads (before any maniputalion of the data), (2) in the 
clean reads (after bbduk) and (3) in the filtered reads (after discarding 
non-viral reads). [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
is run with default parameters. In all three 
cases the quality assesment is performed using fastQC individually for 
each file. After that all fastqc reports in the run are summarized in 
a single report using [multiQC](https://multiqc.info/)(v1.9). 
MultiQC is also used to summarize the bam files quality after alignment of
the reads into the reference sequences using bowtie (see section 2.3).
In both cases MultiQC is run with default parameters. The final _html_ 
reports are copied to the reports directory for the final data visualization.

No parameters can be modified in the quality stages of the pipeline.

More information:
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
* [MultiQC](https://multiqc.info/)

## Cleaning (+ quality)

Reads cleaning is performed using 
[BBDuk(BBMap version 38.96)](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), 
it takes the paired files and process them (always coulpled), it outputs 
a pair of mated files (with the suffix *_pe1* and *_pe2* respectively) 
and also a singletons file (with the suffix *_sgl*). As mentioned above, 
the quality of these files is subsequently assessed with fastQC and mulitQC programs. 

All output files (paired and singleton fastq compressed files and quality 
reports) are saved in the `$BDIR/cleanseqs` directory.

At this point, the reads identifiers is translated from the ilumina ID to
the sample ID (see [data summary section](#data-summary)).

**READS CLEANING PARAMETERS**

* `params.bbdukREF`

By default, the reference file used in the cleaning step contains the reference
illumina adapters and is given as a BBDUK resource. If you use a different 
technology and/or adapters, provide them in a comma separated 
file.  

* `params.bbdukMINLEN`

Reads shorter than this after trimming will be discarded, bbduk uses 10 by
default, in this pipeline the minlength is set to 32. 

* `params.bbduqMAQ`

Minimum average quality, reads with average quality (after trimming) below 
this will be discarded. BBDuk, bu default has no theshols value (maq=0), 
in this pipeline it is set to 10.

(Remember that you can [modify any of the parameters](#configuration-file) 
in the config file or directly from the commandline.)


More information: [BBDUK guide](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)


## Filter and discard non viral reads

Before mapping or assembling steps, we try to discard reads identified as
non-viral species (such as euchariotic, archeea or bacterial DNA). To do so, 
a taxonomical classification is performed using [Kaiju](https://kaiju.binf.ku.dk/)(v.1.9.0).
Paired and single reads are assigned to a taxonomic grup separately and,
later, merged in a same file to generate a tab-separated summary and a 
krona plot in html format (you can find this files in `$BDIR/taxonomy/kaiju/"sample_ID"/reads_taxon`). 
All reads assigned to non-viral are discarded (using seqkit v.0.15.0) and a new set of fastq files 
is saved in the `$BDIR/cleaneqs` directory with the suffix `*.filtered.fastq.gz`). As mentioned above, 
the quality of these files is subsequently assessed with fastQC and mulitQC programs. 

**KAIJU (ON READS) PARAMETERS**

* `.params.kaijuDBRAW`:

The reference database used to assign a taxonomic range to the reads. By 
default it is used *nr_euk* database, it can be changed by *refseq* or even by
a cutomized one. You can find more information about different databases in 
the [kaiju github repository](https://github.com/bioinformatics-centre/kaiju).

* `params.kaijuMAX_FORKS`, `params.kaijuMAX_RETRIES` and `params.kaijuNCPUS`:

The kaiju is the higher resource-consuming program in the pipeline, when 
running too many kaiju proccesses in parallel it can lead to a pipeline crash, 
to avoid so, specific parallelisation parameters are set for the kaju in both reads taxonomy and 
[contigs taxonomy](#taxonomic-classification). `kaijuMAX_FORKS` 
refers to the  maximum number of process instances that can be executed 
in parallel, it is set by default to 4. `kaijuNCPUS` refers to the maximum
number of CPUs used in this process, it is set to 26 by default. The error 
strategy is set to "retry" up to 2 times ( `kaijuMAX_RETRIES`). Tou can find more information about 
directives in [this page of the Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#directives).

(Remember that you can [modify any of the parameters](#configuration-file) 
in the config file or directly from the commandline.)

More information: 
* Kaiju [documentation](https://kaiju.binf.ku.dk/) and 
[github repository](https://github.com/bioinformatics-centre/kaiju).
* Seqkit [documentation](https://bioinf.shenwei.me/seqkit/)


## Reads alignment on reference sequences:

The final set of reads is aligned against the set of reference sequences (the
ones used to designs the capture probes) to see which species were found 
and what is the coverage. The alignment is performed using the 
[bowtie2 software (v.2.4.2)]((http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)) 
and after [samtools (v.1.1)](http://www.htslib.org/doc/samtools.html) is used to 
discard unmapped and low quality reads. This step is performed
separately in paired ends reads and single end reads with the same exact parameters.
The bowtie is run always with the *end-to-end* alignment mode and *sensitive* 
presets, this parameters cannot be modified in this pipeline since many 
false positives and/or low quality hits have been detected with any other approach. 


**ALIGNMENT PARAMETERS**

* `params.alignMINQ`
This parameters refers to the minimum acceptable quality (in phred score format) 
of the aligned reads when filtering with samtools. By default it is set 
to 13 (Q=13), which is translated in an error prbability of 0.05 
(E=10<sup>-0.1Q<sup>=10<sup>-1.3<sup>=0.05) or a 95% of accuracy.

More information: 
* [Bowtie documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
* [Samtools manual](http://www.htslib.org/doc/samtools.html)

## Reads assembly:

The final set of reads is assembled into contigs for the further taxonomic
analysis. 

**ASSEMBLY PARAMETERS**

* `params.assembler`:

This pipeline is developed to support two assemblers: MEGAHIT and TRINITY.
By default the assembly is performed using megahit, which is notably faster 
and less resource-consuming.

*CONSIDERATIONS ON TRINIY ASSEMBLIER:  Trinity is a highly resource-consuming
program and requires higher number of reads than other programs to assembly a contig.
For this reason the maximum memory allowed for this  is this option is 250G (you 
can modify this using the parameter* `params.trinityMAXM`. *To avoid the 
crash of the whole pipeline in the cases when the programme is not able to 
assembly any contig an exception handler have been added to the trinity module; 
if this happens nextflow will raise a warning, for any further information 
you will have to look over log files. Please take into account that this pipeline
has been developed to work with environmental samples (from cities sewage), 
the DNA may be degradated and in some cases trinity may not be able to 
assembly any contig.*

* `params.megahit_errorHandler`:
This paraemeters allows the user to add an error handler in case that 
megahit crashes. This might happen when the program is not able to 
assembly any reads, specially in negative controls or degradated samples.
To avoid this, the error handling strategy can be changed to 'ignore' 
in the configuration file or directly as a CL option.


* `params.assemblyMINCONLEN`:

No contigs shorter than this number will be assemblied. The same parameter
is valid either for megahit or trinity assembliers. Considering the reads
lengths and the origin of the samples, by default, the minimum contig 
length is set to 100.
After assembling the contigs, the reads are aligned using the contigs
set as reference (In this case parameters are --end-to-end --fast) using 
bowtie and samtools. All sequences longer than 100nt that do not align 
to any contig, are considered singletons and added to the final assembly.


More information: 
* [Trinity documentation](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
* [Megahit documentation](https://github.com/voutcn/megahit)


## Blast of contigs:

Contigs are aligned 

**_Program version:_**

More information: 

## Taxonomic classification: 

**KAIJU (ON CONTIGS) PARAMETERS**

* `params.kaijuMAX_FORKS`, `params.kaijuMAX_RETRIES` and `params.kaijuNCPUS`:

The kaiju is the higher resource-consuming program in the pipeline, when 
running too many kaiju proccesses in parallel it can lead to a pipeline crash, 
to avoid so, specific parallelisation parameters are set for the kaju in both [reads taxonomy](#filter-and-discard-non-viral-reads) and 
contigs taxonomy. `kaijuMAX_FORKS` 
refers to the  maximum number of process instances that can be executed 
in parallel, it is set by default to 4. `kaijuNCPUS` refers to the maximum
number of CPUs used in this process, it is set to 26 by default. The error 
strategy is set to "retry" up to 2 times ( `kaijuMAX_RETRIES`). Tou can find more information about 
directives in [this page of the Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#directives).


(Remember that you can [modify any of the parameters](#configuration-file) 
in the config file or directly from the commandline.)

More information: 

* 
* 

## Figures

**Mapped reads summary**

**Coverage Figures**

**KronaPlots**
<!--
Besides the aligned and filtered bam files, this step also produces some summary
plots showing the profile of each sample: A barplot containing the counts of reads 
aligned to each reference sequence is saved in the reports directory and a coverage figure

els outputs d'aquest pas son: coverage figures (see assembly) i 
uns barplots que ens diran quants reads han mapat a cada genoma per mostra (diferenciarem pe de sgl).
-->

<!--
els grafics daquest pas son: els coverage plots: explicar + posar un exemple.
-->

# Dependencies:


* FastQC v0.11.9
* MultiQC v1.9
* BBMap v.38.96
* Kaiju v.1.9.0
* Bowtie2 v.2.4.2
* samtools v.1.11
* htslib v.1.11-4


<!--
Compiling commands:

* MARKDOWN TO HTML

pandoc --from=gfm --to=html --output=readme_DOCUMENTATION_virwaste.html readme_DOCUMENTATION_virwaste.md

* MARKDOWN TO PDF

pandoc --from=gfm --to=latex                          \
       --template=template_de_la_maria.tex            \
       --variable=papersize:a4 --variable=toc:true    \
       --number-section                               \
       --output=readme_DOCUMENTATION_virwaste.tex     \
       readme_DOCUMENTATION_virwaste.md

pdflatex readme_DOCUMENTATION_virwaste.tex

-->
