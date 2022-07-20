#! /usr/bin/env nextflow

process bbduk_clean {

	input:
  
	  val dep
          tuple val(sampid), val(illuid)

	output:
	  val "${illuid}_pe1.fastq.gz", emit: outPE1
	  val "${illuid}_pe2.fastq.gz", emit: outPE2
	  val "${illuid}_sgl.fastq.gz", emit: outSGL

	script:
	
	
	"""
	mkdir -vp $params.clnfq_dir;
	# Adapters Trimming:
	
	\$BBDUK_PATH/bbduk.sh  \
	         in=${sampid}_R1_001.fastq.gz \
		    in2=${sampid}_R2_001.fastq.gz \
		    out=${illuid}_pe1.fastq.gz \
		   out2=${illuid}_pe2.fastq.gz \
		   outs=${illuid}_sgl.fastq.gz \
            ref=\$BBDUK_PATH/resources/adapters.fa \
              k=13 ktrim=r useshortkmers=t mink=5 \
          qtrim=t trimq=20 minlength=32           \
        threads=$params.NCPUS overwrite=true maq=10 \
        stats=${illuid}_stats.out \
             2> ${illuid}.bbduk_clean.log 1>&2 ;
    touch ${illuid}.bbduk_clean.ok; 
	"""

}


