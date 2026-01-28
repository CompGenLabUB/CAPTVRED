#! /usr/bin/env nextflow

process bbduk_clean {
    memory '8.GB' 
	tag "$meta.id"

	input:
	  tuple val(meta), path(reads)

	output:
	  tuple val(meta), path("${meta.id}_trim_{pe1,pe2,sgl}.fastq.gz"), emit: cleanReads
	  // path("${meta.id}_trim_sgl.fastq.gz")                     , emit: cleanSgl
	  path "${meta.id}_stats.out"                              , emit: stats
      path "${meta.id}.bbduk.log"                              , emit: log
	  


	script:
	
	def mem = task.memory ? "-Xmx${task.memory.toGiga()-1}g" : ""
  
	"""

	# Adapters Trimming:
	bbduk.sh  in=${reads[0]} \
          in2=${reads[1]} \
          out=${meta.id}_trim_pe1.fastq.gz \
          out2=${meta.id}_trim_pe2.fastq.gz \
          outs=${meta.id}_trim_sgl.fastq.gz \
          ref=${params.bbdukREF } \
          k=13 ktrim=r useshortkmers=t mink=5 \
          qtrim=t trimq=20 minlength=${params.bbdukMINLEN}           \
          threads=${task.cpus} overwrite=true maq=${params.bbdukMAQ} \
          stats=${meta.id}_stats.out \
         2> ${meta.id}.bbduk.log 1>&2 ;
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
