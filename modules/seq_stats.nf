#! /usr/bin/env nextflow

process fastQC {
    input:
     val samps
     val rawfdir
     val logsdir
     
    output:
      val "$odir", emit: QCout
      
    script:
        // println "## FASTQC on sample: $samps"
        ofl=samps.toString().split("/")[-1].replaceAll(".gz|.sorted.mapped.bam", ".QC")
        odir="$rawfdir/$ofl"
        logfl=samps.toString().split("/")[-1].replaceAll(".gz|.sorted.mapped.bam", ".log")
        logfl="$logsdir/fastQC_$logfl"

    """
    mkdir -vp $odir
    fastqc $samps -o $odir 2> $logfl
    """

}



process multiQC {
    
    input:
     val sps
     val reportsdir
     val logdir
     val type
     
    script:
      sampdirs=sps.join("  ")
      //input is a list, transform it into a string

    """
    multiqc $sampdirs \
        --title "$params.runID MultiQC in ${type} reads"    \
        --fullnames                                         \
        --force                                             \
        --filename ${params.runID}_multiqc_${type}.html     \
        --outdir ${reportsdir}                              \
         2> ${logdir}/${params.runID}_multiqc_${type}.log;
         
    """
}


process multiQC_raw {
    
    input:
     val sps
     val reportsdir
     val logdir
    
    script:
      sampdirs=sps.join("  ")
      //input is a list, transdorm it into a string so miliqc can use it as argument
   
     
    """
    multiqc $sampdirs \
        --title "$params.runID MultiQC in raw reads"    \
        --fullnames                                     \
        --force                                         \
        --filename ${params.runID}_multiqc_raw.html     \
        --outdir ${reportsdir}                          \
         2> ${logdir}/${params.runID}_multiqc_raw.log;
         
    """
}

process multiQC_clean {
    input:
     val sps
     val repdir
     val logdir

    script:
      sampdirs=sps.join("  ")
      //input is a list, transdorm it into a string so miliqc can use it as argument
   
    """
    multiqc $sampdirs \
        --title "$params.runID MultiQC in clean reads"  \
        --fullnames                                     \
        --force                                         \
        --filename ${params.runID}_multiqc_clean.html   \
        --outdir ${repdir}                              \
         2> ${logdir}/${params.runID}_multiqc_clean.log;
    
    """
}


process multiQC_filt {
    input:
     val sps
    
    script:
      sampdirs=sps.join("  ")
      //input is a list, transdorm it into a string so miliqc can use it as argument
    
     
    """
    multiqc $sampdirs \
        --title "$params.runID MultiQC in filtered reads (kaiju non-viral removed)"    \
        --fullnames                                     \
        --force                                         \
        --filename ${params.runID}_multiqc_filt.html    \
        --outdir ${params.reports_dir}                  \
         2> ${params.logs_dir}/${params.runID}_multiqc_filt.log;

    """
}


process multiQC_bowtie { 

    input:
     val sps
     val reportsdir
     val logdir
     val type
    
    script:
     // samplogs=sps.join("  ").replaceAll(".sorted.mapped.bam", ".log")
    //input is a list, transform it into a string so multiqc can use it as argument
    sampdirs=sps.join("  ")
    
    """
    echo "AQUI ESTEMM" > kkfinsh
    multiqc ${sampdirs}  \
        --title "${params.runID} MultiQC in bam files after amplicon align with bowtie"    \
        --fullnames                                     \
        --force                                         \
        --filename ${params.runID}_multiqc_${type}.html   \
        --outdir ${reportsdir}  \
        2> ${logdir}/${params.runID}_multiqc_${type}.log;
    """
}
