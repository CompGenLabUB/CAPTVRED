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
exactly to the first part of the fasta file name (i.e.: followed by R{1,2}_001.fastq.gz). 
Sample ID corresponds to the sample identifier that will be used in all the pipeline steps. 
If you are not interested in modifying the ids, use the same code in both columns.

## 1.3. Reference sequences
__Tabular info file__  __!!! REVISAR !!!__<br />
Descriptive tabular file including the NCBI id, full name, short name and 
taxon ID for each of the reference sequences in the reference set.  


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

__Kaiju Databases__ <br />

Taxonomic classification with kaiju is performed with multiple kaiju 
databases. By default the virwaste pipeline runs the taxonomic analysis 
on 5 diferent databases: nr_euk, refseq, viruses, rvdb and a database 
created with the set of amplicon reference sequences (named capref).

In the Kaiju documentation, the instructions to download the availabe
source databases and to create a new database from a fasta file 
(https://github.com/bioinformatics-centre/kaiju). 



## 1.4. Projectvars.sh file
Includes all environment variables necessary to run the pipeline. When 
this script is executed the whole directories filesystem is created.
Every run must have its projectvars.sh inside its base folder.

* Root directory ($RDIR) is the top level directory of the project.
* Base dirctory ($BDIR) refers to the root folder where all files created 
in each of the steps in the pipeline will be stored.
* Nextflow directory ($NXFDIR) is the directory where virwaste pipeline 
is stored.
* Run identifier ($RUNID) is the especific name of the sequencing experiment. A <br />
* Directory that contains the fasta file with the reference sequences and 
the tabular info file ($AMPSQD). The pipeline will create  bowtie and
blast databases and indexes of the fasta in this folder. <br />
* Name of the file contianing the amplicon reference sequences ($AMPSQFA).
Bowtie and blast databases will be created with the same prefix name. <br />
* Directory where the kaiju databases are stored ($KAIDB) each database 
must be stored in its own directory. __!! REVISAR !!__<br />

This variables are set in the _projectvars.sh_ file by the user.
Recomended filesystem structure would be:

    RDIR<br />
   .   | - NXFDIR<br />
   .   | - BDIR (Run1)<br />
   .   .      | - projectvars.sh<br />
   .   |<br />
   .   | - BDIR (Run2)<br />
   .   .      | - projectvars.sh<br />

After running _projectvars.sh_:
* fastq.gz files must be placed (or linked) to the directory 
"$BDIR/rawseqs_fastq"
* samples_definition.tbl must be placed directly on the $BDIR
* reference sequences (in fa.gz) format must be placed (or linked) in the refseqs directory.


# SECTION 2: Detailed description of the pipeline steps

## 2.1. Quality:

## 2.2. Cleanning (+ quality)

## 2.2. Filter and discard non viral (+ quality)

In the first step (reads cleaning) of the pipeline, the names are 
changed (thus, raw reads keep their original ID while clean reads (and 
everything that comes after) is renamed into de sample ID.
