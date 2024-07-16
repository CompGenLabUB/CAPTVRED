//-a $params.rvdb_update -eq "false"
process merge_rvdb_setref () {

    input:
        val merged_fa
        val rvdb_fa //Full path
        val dbname
        val set_fa  //Full path
        val logfl

    output:
        val merged_fa  //Full path
    
    script:
        //  merged_fa=rvdb_fa.toString().replaceAll(".fasta.gz|.fasta|.fa.gz|.fa", "+set.fasta.gz")
       
      //merged_fa="${odir}/rvdb+${params.setname}.fasta.gz"
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

process chek_setref_ids () {
    
    input:
        val odir
        val rvdb_fa //Full path
        val set_fa  //Full path
        val setname
        
    output:
        val missing_fa  //Full path

    script:
     
     setids="${odir}/${params.setname}.ids"
     missing_ids="${odir}/${params.setname}.missing.ids"
     missing_fa="${odir}/${params.setname}.missing.fa.gz"

     """
     zcat $set_fa  | grep '^>' - | awk '{print \$1}' - | sed 's/>//' > $setids

     zcat $rvdb_fa | awk -vFS="|" ' NR==FNR{
                                      set[\$1];next
                                    }(\$3 in set){
                                        delete set[\$3]
                                    }END{
                                        for (i in set) print i
                                    }' $setids  -  \
             > $missing_ids;

     seqkit grep -f $missing_ids  $set_fa   |\
            sed 's/>/>acc|$setname|/'       |\
            sed 's/ /|/'                    |\
            gzip    -  >  $missing_fa 
     """
}


process filter_FOI () {  //families of interest
    input:
        val odir
        val set_taxfl    //Full path
        val db_tax       //Full path
        val db_fa        //Full path
        val bin          // Bin directory

    output:
        val "$foi_subset_fa",   emit: FOI
        val "$other_subset_fa", emit: OTHER

    script:
        foi_lst=set_taxfl.replaceAll(".tax.gz", "_families.foi") 
        otherID=set_taxfl.replaceAll(".tax.gz", "_families.TEST")
        foi_subset_fa=db_fa.replaceAll(".fasta.gz", "_foi_subset.fasta.gz")
        other_subset_fa=db_fa.replaceAll(".fasta.gz", "_other_subset.fasta.gz")
        logfile=db_fa.replaceAll(".fasta.gz", "_filter.log")

        """
        # Get families of interest list
        zegrep -o  "\\sfamily:[A-Za-z]*:[0-9]*"  $set_taxfl |\
                sort                                        |\
                uniq > $foi_lst;
        ## echo otherID
        #Split fasta 
        $bin/filter_by_family_taxon.py             \
                --fasta_file     $db_fa            \
                --fam_list_file  $foi_lst          \
                --fulltaxon_file $db_tax           \
                --FOI_seqs       $foi_subset_fa    \
                --OTHER_seqs     $other_subset_fa # \
            # 2> $logfile;

        """

}