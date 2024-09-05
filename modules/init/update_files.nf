
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
        wget $full_link -O ${fullpath}  2>> $logfl;
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
        wget $link -O  $fullpath  2>> $logfl;
        tar -xzvf $fullpath  -C ${odir} 2>> $logfl;
        """

}


process get_accession2taxid () {
    
    input:
        val odir
        val link
        val logfl
    
    output:
        val fullpath

    script:
        flname=link.toString().split("/")[-1]
        fullpath="${odir}/${flname}"

        """
        wget $link -O  $fullpath  2>> $logfl;
        """

}