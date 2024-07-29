#! /usr/bin/env nextflow

process bbduk_clean {

	input:
  
	  val dep
      tuple val(illuid), val(sampid)
	  val logdir
	  val REF

	output:
	  val "${sampid}_pe1.fastq.gz", emit: outPE1
	  val "${sampid}_pe2.fastq.gz", emit: outPE2
	  val "${sampid}_sgl.fastq.gz", emit: outSGL



	script:

	logfl=sampid.toString().split("/")[-1]
	logfl="${logdir}/${logfl}.bbduk_clean.log"
	R1=params.R1
	R2=params.R2

	println "### bbduk_clean is running!"
  
	"""
	extension=""
	if [ -f ${illuid}${R1}.fastq.gz ]; then extension=".fastq.gz"; fi;
	if [ -f ${illuid}${R1}.fq.gz ]; then extension=".fq.gz"; fi;

	# Adapters Trimming:
	bbduk.sh  in=${illuid}${R1}\${extension} \
          in2=${illuid}${R2}\${extension} \
          out=${sampid}_pe1.fastq.gz \
          out2=${sampid}_pe2.fastq.gz \
          outs=${sampid}_sgl.fastq.gz \
          ref=${REF} \
          k=13 ktrim=r useshortkmers=t mink=5 \
          qtrim=t trimq=20 minlength=${params.bbdukMINLEN}           \
          threads=${params.NCPUS} overwrite=true maq=${params.bbdukMAQ} \
          stats=${sampid}_stats.out \
         2> $logfl 1>&2 ;
    touch ${sampid}.bbduk_clean.ok; 
	"""

}

process samps_idtranslate {
     input:
       val dep
       tuple val(illuid), val(sampid)

     output:
	
	val "${sampid}_pe1.fastq.gz", emit: outPE1
	val "${sampid}_pe2.fastq.gz", emit: outPE2
	val "${sampid}_sgl.fastq.gz", emit: outSGL

      """
	extension=""
	if [ -f ${illuid}${R1}.fastq.gz ]; then extension=".fastq.gz"; fi;
	if [ -f ${illuid}${R1}.fq.gz ]; then extension=".fq.gz"; fi;
	
        cp  ${illuid}${R1}\${extension}  ${sampid}_pe1.fastq.gz;
        cp  ${illuid}${R2}\${extension}  ${sampid}_pe2.fastq.gz;
        touch ${sampid}_sgl.fastq.gz;
       
      """
}
