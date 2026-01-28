#! /usr/bin/env nextflow

process fastQC {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    // val samps
    // val rawfdir
    // val logsdir
     
    output:
      // val "$odir", emit: QCout
        tuple val(meta), path("*.html"), emit: html
        tuple val(meta), path("*.zip") , emit: zip
        path "*.log"                   , emit: log

    script:
        // println "## FASTQC on sample: $samps"
        // ofl=samps.toString().split("/")[-1].replaceAll(".gz|.sorted.mapped.bam", ".QC")
        // odir="$rawfdir/$ofl"
        // logfl=samps.toString().split("/")[-1].replaceAll(".gz|.sorted.mapped.bam", ".log")
        // logfl="$logsdir/fastQC_$logfl"


    """

     fastqc --quiet --threads ${task.cpus} $reads 2> ${meta.id}.fastqc.log
    """

}



process multiQC {
    tag "all_samples"

    input:
     path qc_files
     val type
    
    output:
    path "*.html", emit: html
    path "*.log" , emit: log
     
    script:
      // sampdirs=sps.join("  ")
      //input is a list, transform it into a string

    """
         multiqc .  --title "$params.runID MultiQC in ${type} reads"  \
                    --filename ${params.runID}_multiqc_${type}.html   \
                    2>  ${params.runID}_multiqc_${type}.log  
    """
}


process fq2fasta {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    // Mantenemos el meta para que FCSGX sepa qué muestra es
    tuple val(meta), path("${meta.id}.fasta"), emit: fasta

    script:
    """
    seqkit fq2fa ${reads.join(' ')} | sed 's/ /_/'  > ${meta.id}.fasta
    """
}


process discard_contaminants {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(contaminants)

    output:
    tuple val(meta), path("*_cleaned.fastq.gz"), emit: cleaned_reads

    script:
    def r1_src = reads[0]
    def r2_src = reads[1]
    def sg_src = reads[2]
    """
    awk '\$5=="EXCLUDE"{split(\$1, id, "_"); print id[1]}' $contaminants \
              > ${meta.id}_nonviral.ids
    if [[ -s $r1_src ]]; then
        seqkit grep -v -f ${meta.id}_nonviral.ids ${r1_src} | gzip -9 > ${meta.id}_pe1_cleaned.fastq.gz
    else 
        cp $r1_src ${meta.id}_pe1_cleaned.fastq.gz;
    fi;

    if [[ -s $r2_src ]]; then 
        seqkit grep -v -f ${meta.id}_nonviral.ids ${r2_src} | gzip -9 > ${meta.id}_pe2_cleaned.fastq.gz
    else
        cp $r2_src ${meta.id}_pe2_cleaned.fastq.gz
    fi;

    if [[ -s $sg_src ]]; then
        seqkit grep -v -f ${meta.id}_nonviral.ids ${sg_src} | gzip -9 > ${meta.id}_sgl_cleaned.fastq.gz
    else
        cp $sg_src ${meta.id}_sgl_cleaned.fastq.gz
    fi;
    """

}

 //  process multiQC_raw {
 //      
 //      input:
 //       val sps
 //       val reportsdir
 //       val logdir
 //      
 //      script:
 //        sampdirs=sps.join("  ")
 //        //input is a list, transdorm it into a string so miliqc can use it as argument
 //     
 //       
 //      """
 //      multiqc $sampdirs \
 //          --title "$params.runID MultiQC in raw reads"    \
 //          --fullnames                                     \
 //          --force                                         \
 //          --filename ${params.runID}_multiqc_raw.html     \
 //          --outdir ${reportsdir}                          \
 //           2> ${logdir}/${params.runID}_multiqc_raw.log;
 //           
 //      """
 //  }
 //  
 //  process multiQC_clean {
 //      input:
 //       val sps
 //       val repdir
 //       val logdir
 //  
 //      script:
 //        sampdirs=sps.join("  ")
 //        //input is a list, transdorm it into a string so miliqc can use it as argument
 //     
 //      """
 //      multiqc $sampdirs \
 //          --title "$params.runID MultiQC in clean reads"  \
 //          --fullnames                                     \
 //          --force                                         \
 //          --filename ${params.runID}_multiqc_clean.html   \
 //          --outdir ${repdir}                              \
 //           2> ${logdir}/${params.runID}_multiqc_clean.log;
 //      
 //      """
 //  }
 //  
 //  
 //  process multiQC_filt {
 //      input:
 //       val sps
 //      
 //      script:
 //        sampdirs=sps.join("  ")
 //        //input is a list, transdorm it into a string so miliqc can use it as argument
 //      
 //       
 //      """
 //      multiqc $sampdirs \
 //          --title "$params.runID MultiQC in filtered reads (kaiju non-viral removed)"    \
 //          --fullnames                                     \
 //          --force                                         \
 //          --filename ${params.runID}_multiqc_filt.html    \
 //          --outdir ${params.reports_dir}                  \
 //           2> ${params.logs_dir}/${params.runID}_multiqc_filt.log;
 //  
 //      """
 //  }
 //  
 //  
 //  process multiQC_bowtie { 
 //  
 //      input:
 //       val sps
 //       val reportsdir
 //       val logdir
 //       val type
 //      
 //      script:
 //       // samplogs=sps.join("  ").replaceAll(".sorted.mapped.bam", ".log")
 //      //input is a list, transform it into a string so multiqc can use it as argument
 //      sampdirs=sps.join("  ")
 //      
 //      """
 //      echo "AQUI ESTEMM" > kkfinsh
 //      multiqc ${sampdirs}  \
 //          --title "${params.runID} MultiQC in bam files after amplicon align with bowtie"    \
 //          --fullnames                                     \
 //          --force                                         \
 //          --filename ${params.runID}_multiqc_${type}.html   \
 //          --outdir ${reportsdir}  \
 //          2> ${logdir}/${params.runID}_multiqc_${type}.log;
 //      """
 //  }
 //  