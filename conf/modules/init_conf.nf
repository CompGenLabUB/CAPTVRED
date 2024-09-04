process create_logf () {
input:
    val refseqsdir
output:
    val logfl

script:

    logfl="${refseqsdir}/init_captvred_dbs.log"
    def now = new Date().format('dd-MM-yyyy')

    """
    #check if logfl exists, if not, create new folder
    [[ -f $logfl ]] ||  ( touch $logfl );
    echo "### Preparing DATABASES FOR CAPTVRED : $now ###" > $logfl
    """
}

process db_for_kaiju_cus () {

    input:
        val dbname
        val logfl

    output:


    script:
        def now = new Date().format('dd-MM-yyyy')

        """ 
        echo "### Downloading $dbname for kaiju : $now ###" >> $logfl
        kaiju-makedb -s $dbname >>  $logfl 2>&1
        """
}

process db_for_kaiju_pred () {
     input:
        val odir
        val datab
       // val flname
        val logfl
    
    output:
        val fullpath

    script:
        // flname=link.toString().split("/")[-1]
        fullpath="${odir}/${datab}"
        println "DB is: $datab"
        println "Link is: $params.K_nr_euk_link"
        def link = ":o"
        if (datab == "nr_euk"){
            link = params.K_nr_euk_link
        } else if (datab == "refseq") {
            link = params.K_refseq_link
        } else if (datab == "viruses") {
            link = params.K_viruses_link
        } else if (datab == "rvdb") {
            link = params.K_rvdb_link
        } else {
            println "WARNING!! Unknown kaiju database! \n Available options are: nr_euk, refseqs, viruses, rvdb"
        }
        println "Link is: $link" 
        """
        mkdir -vp $odir;
        wget $link -O  $fullpath  2>  $logfl;
        """

}

process stdrd_link () {   // Create a link with stable name so it can be called from main workflow

    input:
        val origin_file
        val dest_link
        val logfl

    output:
        val dest_link  //Full path
    
    script:
      """
        ln -fs $origin_file  $dest_link >> $logfl 2>&1
        echo "    ... link done!" >> $logfl;
      """
}