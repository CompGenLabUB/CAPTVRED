#! /usr/bin/env nextflow

process fastQC {
    input:
     val x
     
    output:
      val "$odir", emit: QCout
      
    script:
        //println "## FASTQC on sample: $x"
        odir=x.toString().replaceAll(".gz", ".QC")
    """
    #aa=$x
    #odir=\${aa%%.gz}".QC"
    mkdir -vp $odir
    fastqc $x -o $odir
    """

}

process multiQC_raw {
    input:
     val x
    script:
   
      //input is a list, transdorm it into a string so miliqc can use it as argument
    sampdirs=x.join("  ")
    // println "## MultiQC on dirs: $sampdirs"
     
    """
    multiqc $sampdirs \
        --title "$params.runID MultiQC in raw reads"    \
        --fullnames                                     \
        --force                                         \
        --filename ${params.runID}_multiqc_raw.html   \
        --outdir ${params.reports_dir}   \
         2> ${params.reports_dir}/${params.runID}_multiqc_raw.html;
        
    """
}

process multiQC_clean {
    input:
     val x
    script:
   
      //input is a list, transdorm it into a string so miliqc can use it as argument
    sampdirs=x.join("  ")
    //println "## MultiQC on dirs: $sampdirs"
     
    """
    multiqc $sampdirs \
        --title "$params.runID MultiQC in clean reads"    \
        --fullnames                                     \
        --force                                         \
        --filename ${params.runID}_multiqc_clean.html   \
        --outdir ${params.reports_dir}  \
         2> ${params.reports_dir}/${params.runID}_multiqc_clean.log;
        
    """
}

process multiQC_bowtie_amp { 

    input:
     val x

    script:
   
    //input is a list, transform it into a string so multiqc can use it as argument
    //sampdirs=x.join("  ")
    samplogs=x.join("  ").replaceAll(".sorted.bam", ".log")
    println "## MultiQC on dirs: $samplogs"
    """
    multiqc $samplogs \
        --title "$params.runID MultiQC in bam files after amplicon align with bowtie"    \
        --fullnames                                     \
        --force                                         \
        --filename ${params.runID}_multiqc_bowtie_amplicons_alignment.html   \
        --outdir ${params.reports_dir}  \
        2> ${params.reports_dir}/${params.runID}_multiqc_bowtie_amplicons_alignment.log;
    """
}
