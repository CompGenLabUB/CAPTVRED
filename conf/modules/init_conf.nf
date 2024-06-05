process create_logf {
input:
    val refseqsdir
output:
    val logfl

script:

    logfl="${refseqsdir}/init_captvred_dbs.log"
    def now = new Date().format('dd-MM-yyyy')

    """
    #check if logfl exists, if not, create new folder
    [[ -f $logfl ]] ||  ( touch $logfl )
    echo "### Preparing DATABASES FOR CAPTVRED : $now ###"
    """
}

process db_for_kaiju (){

    input:
        val dbname
        val logfl

    output:


    script:

        """
        kaiju-makedb -s $dbname >>  $logfl 2>&1
        """
}



process init_rvdb () {

    input:
        val odir
        val link
        val dbname
        val logfl

    output:
        val fullpath

    script:
        full_link="${link}/${dbname}"
        fullpath="${odir}/${dbname}"

        """
        mkdir -vp "${odir}";
        echo "   > Downloading $dbname database" >> $logfl;
        if $params.rvdb_update; then 
            wget $full_link -P ${odir}  2>> $logfl;
        else
            if [ -f $fullpath ]; then
                echo "    ... \"$dbname\" is already downloaded." >> $logfl;
            else 
                wget $full_link -P ${odir}  2>> $logfl;
            fi;
        fi;
        """
}

//-a $params.rvdb_update -eq "false"
process merge_rvdb_setref () {

    input:
        val odir
        val rvdb_fa //Full path
        val dbname
        val set_fa  //Full path
        val logfl

    output:
        val merged_fa  //Full path
    
    script:
        //  merged_fa=rvdb_fa.toString().replaceAll(".fasta.gz|.fasta|.fa.gz|.fa", "+set.fasta.gz")
        merged_fa="${odir}/rvdb+${params.setname}.fasta.gz"
        """
         echo "   > Mergeing $dbname + $params.setname" >> $logfl;
         if  $params.merge_update; then
            zcat $rvdb_fa $set_fa | gzip  -  > $merged_fa;
            echo "    ... done!";
         else
            if [  -f $merged_fa  ]; then
                echo "    ... \"$dbname\" + \"$params.setname\" is already merged." >> $logfl;
            else
                zcat $rvdb_fa $set_fa | gzip  -  > $merged_fa;
                echo "    ... done!";
            fi;
         fi;
        """
}