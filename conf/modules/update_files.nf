
// DOWNLOAD RVDB
process get_rvdb () {

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
            wget $full_link -O ${fullpath}  2>> $logfl;
        else
            if [ -f $fullpath ]; then
                echo "    ... \"$dbname\" is already downloaded." >> $logfl;
            else 
                wget $full_link -O ${fullpath}  2>> $logfl;
            fi;
        fi;
        """
}


process get_names_and_nodes () {
    
    input:
        val odir
        val link
        val logfl
    
    output:
        val odir

    script:
        flname=link.toString().split("/")[-1]
        fullpath="${odir}/${flname}"

        """
        mkdir -vp "${odir}";
        echo "   > Downloading $flname database" >> $logfl;
        if $params.taxon_update; then 
            wget $link -O  $fullpath  2>> $logfl;
            tar -xzvf $fullpath  2>> $logfl;
        else 
            if [ -f $fullpath ]; then
                echo "    ... \"$flname\" is already downloaded." >> $logfl;
            else 
                wget $link -O $fullpath  2>> $logfl;
                tar -xzvf $fullpath  2>> $logfl;
            fi;
        fi;
        """

}