#! /usr/bin/env nextflow

process fastQC {
    input:
     val x
     
    //output:
    
    script:
        //println "## FASTQC on sample: $x"
    """
    aa=$x
    odir=\${aa%%.gz}".QC"
    mkdir -vp \$odir
    fastqc $x -o \$odir
    """

}

//process multiQC{

//}
