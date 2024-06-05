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

/*
process db_for_kaiju (){

    input:
        val dbname

    output:

    script:

        logfl="$params.wdir/logs/prepare_db_for_kaijufilter_{$dbname}.log"
        """
        kaiju-makedb -s $dbname >  $logfl 2>&1

        """

}


process init_rvdb () {

    input:
        val link
        val dbname

    output:
        val fullpath

    script:
        logfl="${params.wdir}/logs/init_dbs.log" 
        full_link="${link}/${dbname}"
        fullpath="${params.refsqs_dir}/${dbname}"

        """
         echo "## ## DOWNLOADING DATABASE" >> $logfl
         wget $full_link -P $params.refsqs_dir 2>> $logfl
        """

}


process merge_rvdb_setref () {

    input:
        val rvdb_fa //Full path
        val set_fa  //Full path

    output:
        val merged_fa  //Full path
    
    script:
        //  merged_fa=rvdb_fa.toString().replaceAll(".fasta.gz|.fasta|.fa.gz|.fa", "+set.fasta.gz")
        merged_fa="${params.refsqs_dir}/rvdb+${params.setname}.fasta.gz"
        """
         cat $rvdb_fa $set_fa | gzip  -  > $merged_fa
        """

}
*/