#! /usr/bin/env nextflow

process bbduk_clean {

	input:
  
	  val dep
          tuple val(illuid), val(sampid)

	output:
	  val "${sampid}_pe1.fastq.gz", emit: outPE1
	  val "${sampid}_pe2.fastq.gz", emit: outPE2
	  val "${sampid}_sgl.fastq.gz", emit: outSGL

	script:
  logfl=sampid.toString().split("/")[-1]
	logfl="${params.logs_dir}/${logfl}.bbduk_clean.log"
  
	"""
	# Adapters Trimming:
	
	\$BBDUK_PATH/bbduk.sh  \
	         in=${illuid}${R1}.fastq.gz \
		    in2=${illuid}${R2}.fastq.gz \
		    out=${sampid}_pe1.fastq.gz \
		   out2=${sampid}_pe2.fastq.gz \
		   outs=${sampid}_sgl.fastq.gz \
            ref=${params.bbdukREF} \
              k=13 ktrim=r useshortkmers=t mink=5 \
          qtrim=t trimq=20 minlength=${params.bbdukMINLEN}           \
        threads=$params.NCPUS overwrite=true maq=${params.bbdukMAQ} \
        stats=${sampid}_stats.out \
             2> $logfl 1>&2 ;
    touch ${sampid}.bbduk_clean.ok; 
	"""

}


