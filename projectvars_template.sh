## 

  #####  RUN INFORMATION  #####

# Run identifier:
export RUNID="";

# Samples suffix:
export R1="_R1_001"
export R2="_R2_001"

  ##### PROJECT PATHS AND FILESYSTEM #####

# Base "root" dir
export RDIR="";

# directory where to find nextflow scripts and config files.
export NXFDIR=""; 
export BIND="$NXFDIR/bin"

# Directory with the data: samples_definition.tbl, rawseqs_fastq dir, ...
export BDIR=""; #NEXTFLOW MUST BE RUN FROM THIS DIRECTORY
    #samples_definition.tbl must be placed in this folder.
    #projectvars.sh must be placed in this folder

#Reference sequences and databases:
export REFSQSD="${RDIR}/refseqs"

export BL_DB_D="${REFSQSD}/rvdb_nt" #Blast ref database directory
export BL_DB_NM="C-RVDBvCurrent.fasta.gz"
export BL_DB_IF="C-RVDB_allentries_info.txt"

export KAI_DB_D="${REFSQSD}/kaijudb";  #Directory where the kaiju databases are installed
export KAI_DB_NM="refseq";

export BBDUK_RES="${REFSQSD}/bbmap/resources" # resources required for bbduk.sh to work

export AMPSQD="${REFSQSD}/ampliconseqs"; #directory that contains the fasta file of refseqs.
export AMPSQFA="Viral_candidates_zoonosis.fa.gz"; #Name of the compressed fasta file that contians the sequences of all amplicons used to design the probes.
export AMPSQINFO="Viral_candidates_zoonosis_full_info.tsv";
export AMPGFFD=""; #directory that contains the gff file of refseqs.



   #####  OPTIONS AND VARIABLES  #####
# Computation:
export NCPUS=32

# Assembly:
export minconlen=100;   #Minimum contigs length

#blast
export DBsize=5000000000;
# Handle potential contamination:
contamin_ids=""; #NCBI identifiers of the refseqs of the potential contaminant sequence.
contamin_faa="";  #Protein fasta file of the potentially contaminant sequence.

# Filesystem required by workflow:
#Do not modify!!
export TMPD="$BDIR/tmp"
mkdir -vp $TMPD

export RAWFQ="$BDIR/rawseqs_fastq"
mkdir -vp $RAWFQ
     #fastq data mult be placed here before running nextflow.

export CLNDIR="$BDIR/cleanseqs"
mkdir -vp $CLNDIR

export AMPALD="$BDIR/amplicons_alignment";  #amplicons alignment directory
mkdir -vp $AMPALD

export ASSBLD="$BDIR/assembly";  #assembly directory
mkdir -vp $ASSBLD/megahit
mkdir -vp $ASSBLD/trinity

export CBLASTD="$BDIR/contigs_blast";
mkdir -vp $CBLASTD/blastn
mkdir -vp $CBLASTD/tblastx

export TAXDIR="$BDIR/taxonomy";
mkdir -vp $TAXDIR/reads_taxon
mkdir -vp $TAXDIR/contigs_taxon/blastn
mkdir -vp $TAXDIR/contigs_taxon/tblastx
mkdir -vp $TAXDIR/contigs_taxon/kaiju
mkdir -vp $TAXDIR/taxon_viral_candidates

export RPTD="$BDIR/reports";
mkdir -vp $RPTD;
mkdir -vp $RPTD/coverage_figures

export LOGD="$BDIR/logs";
mkdir -vp $LOGD
#Do not modify!!


##
## NextFlow



# directory where working files are stored
export NXF_HOME="$BDIR/.nextflow";
export NXF_WORK="$BDIR/work";
export NXF_HTML="$RPTD/$RUNID.html";
export NXF_TIMELINE="$RPTD/${RUNID}_timeline.html";
export NXF_DAG="$RPTD/${RUNID}_dag.dot";
# mkdir -vp $NXF_HOME $NXF_WORK; # those folders are created automatically by NextFlow

# directories where used programes are installed:
#export BBDUK_PATH="/usr/local/install/bbmap";

##
## Global variables
#export FASTAQ="rawseqs_fastq/{ILLUMINA_ID}_R[12]_001.fastq.gz";
