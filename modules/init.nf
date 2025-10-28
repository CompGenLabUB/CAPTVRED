#! /usr/bin/env nextflow

process create_logd {

input:
    val newdir

output:
    val logfl

script:

    println "AA"
    logfl="$newdir/create_filesystem.log"
    println "-- $logfl --"

    """
    
    #check if logfl exists and delete it
    [[ -f $logfl ]] && rm $logfl
    mkdir -vp $newdir
    touch $logfl
    """
}

process create_filesys {

input:
    val newdir
    val logfl

output:
    val logfl

script:

    """
    echo "## NAME IS $newdir!!!" >> $logfl;
    mkdir -vp $newdir >> $logfl 2>&1 
    """

}
