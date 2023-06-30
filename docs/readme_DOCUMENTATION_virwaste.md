# Prepare the environment
## Projectvars.sh file

Includes all environment variables necessary to run the pipeline. When 
this script is executed the whole directories filesystem is created.
Every run must have its projectvars.sh inside its base folder.

* Run identifier ($RUNID) is the specific name of the sequencing experiment.  <br />
* $R1 and $R2, which correspond to the pattern used to describe the first and second pair of
files respectively. (see [Input Data section](#input-data))

* Root directory ($RDIR) is the top-level directory of the project.
* Base directory ($BDIR) refers to the root folder where all files created 
in each of the steps in the pipeline will be stored.
* Nextflow directory ($NXFDIR) is the directory where CAPTVRED pipeline 
is stored.
* Directory that contains the fasta file with the reference sequences and 
the tabular info file ($AMPSQD). The pipeline will create  bowtie and
blast databases and indexes of the fasta in this folder. <br />
* Name of the file containing the amplicon reference sequences ($AMPSQFA).
Bowtie and blast databases will be created with the same prefix name. <br />
* Directories where the databases are stored: $BDBD (Blast reference database directory), $AMPSQD (Directory that contains the fasta file of REFSEQS), $KAIDBD(Directory where the kaiju databases are installed).

These variables are set in the _projectvars.sh_ file by the user.
Recommended filesystem structure would be:

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
* reference sequences (in fa.gz) format must be placed (or linked) in the REFSEQS directory.

## Configuration file:

The configuration file is named *nextflow.config*, you can find it in the workflow main directory. 
It contains the default variables used in each of the steps of the pipeline. 
If you are interested in changing any of the parameters before running the pipeline 
you can:

> **a) Modify the document**
All parameters can be modified directly on the configuration file. 
Preferably, do not modify the filename keep it in the same directory. 
If you change the filename or directory of the configuration file, or if
you are running the pipeline using the GitHub link, you can provide the name 
of the configuration file as a command line argument using the `-params-file` option.


> **b) Add the parameter as a command line option**

If you are interested in adjusting only one or a few parameters it might be more
convenient to provide these options directly as command line options. You can do so 
by simply writing: `--ParameterOfInterest NewValue`.
i.e. If you are interested in changing the default number of CPUS to 16 you
can run the pipeline as follows:

```{.sh}
nextflow run https://github.com/JosepFAbril/virwaste --NCPUS 16
```

You can find more information on how to custom and provide the configuration file in the
 [nextflow config documentation](https://www.nextflow.io/docs/latest/config.html) .

## Data summary:

A tabular file containing information about the samples is necessary to 
run the pipeline (only the samples described in this tabular file will 
be analyzed). This file has been written manually, it must be named 
"samples_definition.tbl" and has the following format:

```
#SAMPLE_ID	   ILLUMINA_ID	   TECHNOLOGY	  PAIRED	 SAMPLE_FACTOR	  METHOD_FACTOR	   DESCRIPTION
R01_C09_G01	   G1_C_S9	       Illumina	    PE	     Bat_guano	      Capture	         Sample 9 desc
R01_C10_G02	   G2_C_S10	       Illumina	    PE	     Bat_guano	      Capture	         Sample 10 desc
```

Only the first two columns are used in the pipeline, but it is recommended 
to add some metadata information in the tabular file for further data visualizations.
**Illumina ID** is the identifier given in the sample files and corresponds 
exactly to the root of the raw fastq file name. 
**Sample ID** corresponds to the sample identifier that will be used in all the pipeline steps.
(In the first step of the pipeline (reads cleaning), the names are 
changed in such a way that raw reads keep their original ID while clean reads (and 
everything that comes after) are renamed into de sample ID.) 
If you are not interested in modifying the ids, use the same code in both columns.

## Input data:

All samples must be in fastq.gz format. All of them must be placed in the 
same directory. 

It is assumed that input sequence files are paired files named following
a pattern. The root name for each sample must be described in the [data summary](#data-summary),
and the pattern for the suffix referring to each pair of reads must be described
in the [ *projectvars.sh* file](#projectvarssh-file) ($R1 and $R2 variables).

e.g. if your set of samples is:

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
In this manner, you simplify and standardize the name of the files generated 
all along the pipeline. 

And, finally, you should define the pattern in the projectvars as `R1="_r01"`, since this
is the pattern for the first file of each pair, and `R1="_r02"`, since this
is the pattern for the second file of each pair.


## Reference sequences

__Tabular info file__ <br />
Descriptive tabular file including the NCBI id, full name, short name and 
taxon ID for each of the reference sequences in the reference set.  

camps= ##Famlily	Specie	Host	SeqID	Region	Size	Name	Description
```
```

__Fasta file__  <br />
Fasta file containing the reference sequences. Identifier must be the NCBI ID 
and optionally other information (separated with a space!). 
These sequences will be used as a reference database for alignments (performed 
with bowtie for the reads and blast for the contigs) and for the taxonomic 
classification (see "Kaiju databases"  subsection).


__GFF files__  __!! REVISAR !!__<br />
A separate GFF file 
These files can be obtained by getting the genebank information of the genome, 
extracting the gff3 and selecting the coding regions.

This is the command we used to fetch our _gff_ sequences:
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
Where _Viral_candidates_zoonosis_refseqs_idrelation.tsv_ is the tabular info file mentioned in the [previous subsection](#reference-sequences).

__Kaiju Databases__ <br />

Taxonomic classification with kaiju is performed by default using the RVDB database as default. 
However, it can be changed to any other set of sequences. "nr_euk", "refseq" and "viruses" databases are provided ready to use in the [kaiju materials website](https://bioinformatics-centre.github.io/kaiju/downloads.html)

# Run the pipeline

## Use the GitHub link

```{.sh}
nextflow run https://github.com/JosepFAbril/virwaste 
```

## Download and install manually the pipeline

Clone the repository from GitHub

```{.sh}
nextflow run CAPTVRED [optional parameters] -with-report $RPTDR/Nextflow_execution_report.html
```

## Rerun
If the pipeline crashes at some point it can continue the execution from the last cached results by using the argument *-resume*.

```{.sh}
nextflow run CAPTVRED [optional parameters] -with-report $RPTDR/Nextflow_execution_report.html -resume
```


# Detailed description of the pipeline steps

## Quality:

Reads quality is performed at four different stages during the pipeline, 
(1) in the raw reads (before any manipulation of the data), (2) in the 
clean reads (after the cleaning step) and (3) in the filtered reads (after discarding 
non-viral reads). [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
is run with default parameters. In all three 
cases, the quality assessment is performed using FastQC individually for 
each file. All FastQC reports in the run are summarized in 
a single report using [multiQC](https://multiqc.info/)(v1.9). 
MultiQC is also used to summarize the bam files quality after the alignment of
the reads into the reference sequences using bowtie (see section 2.3).
In both cases, MultiQC is run with default parameters. The final _html_ 
reports are copied to the reports directory for the final data visualization.

No parameters can be modified in the quality stage of the pipeline.

More information:
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
* [MultiQC](https://multiqc.info/)

## Cleaning

Reads cleaning is performed using 
[BBDuk(BBMap version 38.96)](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/), 
it takes the paired files and processes them (always coupled), and it outputs 
a pair of mated files (with the suffix *_pe1* and *_pe2* respectively) 
and also a singletons file (with the suffix *_sgl*). As mentioned above, 
the quality of these files is subsequently assessed with FastQC and MulitQC programs. 

All output files (paired and singleton fastq compressed files and quality 
reports) are saved in the `$BDIR/cleanseqs` directory.

At this point, the reads identifiers are translated from the Illumina ID to
the sample ID (see [data summary section](#data-summary)).

**READS CLEANING PARAMETERS**

* `params.bbdukREF`

By default, the reference file used in the cleaning step contains the reference
Illumina adapters and is given as a BBDUK resource. If you use a different 
technology and/or adapters, provide them in a comma-separated 
file.  

* `params.bbdukMINLEN`

Reads shorter than this after trimming will be discarded, BBDuk uses 10 by
default, in this pipeline the minlength is set to 32. 

* `params.bbduqMAQ`

Minimum average quality, reads with average quality (after trimming) below 
this will be discarded. BBDuk, by default has no threshold value (maq=0), 
in this pipeline, it is set to 10.

(Remember that you can [modify any of the parameters](#configuration-file) 
in the config file or directly from the command line.)

More information: [BBDUK guide](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)


## Filter and discard non-viral reads

Before mapping or assembling steps, we try to discard reads identified as
non-viral species (such as eukaryotic, archaea, or bacterial DNA). To do so, 
a taxonomical classification is performed using [Kaiju](https://kaiju.binf.ku.dk/)(v.1.9.0).
Paired and single reads are assigned to a taxonomic group separately and,
later, merged in the same file to generate a tab-separated summary and a 
krona plot in HTML format (you can find these files in `$BDIR/taxonomy/kaiju/"sample_ID"/reads_taxon` directory). 
All reads assigned to non-viral are discarded (using seqkit v.0.15.0) and a new set of FASTQ files 
is saved in the `$BDIR/cleaneqs` directory with the suffix `*.filtered.fastq.gz`). As mentioned above, 
the quality of these files is subsequently assessed with FastQC and MulitQC programs. 

**KAIJU (ON READS) PARAMETERS**

* `.params.kaijuDBRAW`:

The reference database is used to assign a taxonomic range to the reads. By 
default it is used *nr_euk* database, it can be changed by *refseq* or even by
a customized one. You can find more information about different databases in 
the [kaiju GitHub repository](https://github.com/bioinformatics-centre/kaiju).

* `params.kaijuMAX_FORKS`, `params.kaijuMAX_RETRIES` and `params.kaijuNCPUS`:

Kaiju is the higher resource-consuming program in the pipeline, when 
running too many kaiju processes in parallel the pipeline may crash, 
to avoid so, specific parallelization parameters are set for the kaju in both reads taxonomy and 
[contigs taxonomy](#taxonomic-classification). `kaijuMAX_FORKS` 
refers to the  maximum number of process instances that can be executed 
in parallel, it is set by default to 4. `kaijuNCPUS` refers to the maximum
number of CPUs used in this process, it is set to 26 by default. The error 
strategy is set to "retry" up to 2 times ( `kaijuMAX_RETRIES`). You can find more information about 
directives on [this page of the Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#directives).

(Remember that you can [modify any of the parameters](#configuration-file) 
in the config file or directly from the command line.)

More information: 
* Kaiju [documentation](https://kaiju.binf.ku.dk/) and 
[github repository](https://github.com/bioinformatics-centre/kaiju).
* Seqkit [documentation](https://bioinf.shenwei.me/seqkit/)


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



## Taxonomic classification: 

### Nucleotide-level classification:
### Protein-level classification:

## Reads alignment on reference sequences:

The final set of contigs is aligned against the set of reference sequences (the
ones used to design the capture probes) to see which species were found 
and what is the coverage. The alignment is performed using the 
[bowtie2 software (v.2.4.2)]((http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)) 
and after [samtools (v.1.1)](http://www.htslib.org/doc/samtools.html) is used to 
discard unmapped and low quality reads. This step is performed
separately in paired ends reads and single-end reads with the same exact parameters.
The bowtie is run always with the *end-to-end* alignment mode and *sensitive* 
presets, these parameters cannot be modified in this pipeline since many 
false positives and/or low-quality hits have been detected with any other approach. 


**ALIGNMENT PARAMETERS**

* `params.alignMINQ`
This parameter refers to the minimum acceptable quality (in phred score format) 
of the aligned reads when filtering with *samtools*. By default, it is set 
to 13 (Q=13), which is translated in an error probability of 0.05 
(E=10<sup>-0.1Q<sup>=10<sup>-1.3<sup>=0.05) or a 95% of accuracy.

More information: 
* [Bowtie documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
* [Samtools manual](http://www.htslib.org/doc/samtools.html)

<!--
## Figures

**Mapped reads summary**
**Coverage Figures**
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
