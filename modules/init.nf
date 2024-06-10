#! /usr/bin/env nextflow

process create_logd {
    input:
    val newdir

output:
    val logfl

script:

    logfl="$newdir/create_filesystem.log"

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
