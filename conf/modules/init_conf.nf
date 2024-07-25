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

process db_for_kaiju () {

    input:
        val dbname
        val logfl

    output:


    script:
        def now = new Date().format('dd-MM-yyyy')

        """ 
        echo "### Downloading $dbname for kaiju : $now ###" > $logfl
        kaiju-makedb -t $params.NCPUS -s $dbname >>  $logfl 2>&1
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
      println "   > Linking $origin_file to $dest_link :)"
      """
         if [ -f $origin_file ]; then
            ln -fs $origin_file  $dest_link >> $logfl 2>&1
            echo "    ... link done!" >> $logfl;
         else
            echo "## ERROR!! $origin_file does not exist!!";
         fi;
      """
}