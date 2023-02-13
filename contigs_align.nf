#! /usr/bin/env nextflow

process make_db_for_blast {
    input:
        val(refseqs_fasta)
        val(redodb)

    output:

        val "${out_db}", emit: DB
        val "${is_db}", emit: CTRL

    script:
    
        out_db=refseqs_fasta.toString().replaceAll(".fa.gz|.fa", "_blastdb")
        dblist=out_db.split('/')
        db_name=dblist[-1]
        cdfile=[dblist[0..dblist.size()-2].join('/'), "blastdb_created_genome.cdate"].join('/')
        spid=refseqs_fasta.toString().split('/')[-1].split('[.]')[0]
       // fa_fl= new File("${refseqs_fasta}")
       """
        echo "MAKEDB  -- $refseqs_fasta   \n $out_db " >> ${params.wdir}/kk.$spid;
        if [ -s ${refseqs_fasta} ]; then
              if [ -e ${cdfile} -a ${redodb} == "FALSE" ]; then
                 cat ${cdfile} 1>&2;
                 echo "AAA" >> ${params.wdir}/kk;
                    
              else
                     if [[ \$(echo ${refseqs_fasta}) =~ .gz\$ ]]; then
                         gunzip -c ${refseqs_fasta}               | \
                           makeblastdb -in -  -dbtype nucl         \
                                       -title \"${db_name}\"       \
                                       -out ${out_db}              \
                                       2> ${out_db}.log  1>&2;
                                       
                           date +"DB created on %Y/%m/%d %T %Z %s"    \
                              > ${cdfile};
                           echo "BBB" >> ${params.wdir}/kk;
                     else
                         makeblastdb -in ${refseqs_fasta}  -dbtype nucl    \
                                      -title \"${db_name}\"                 \
                                      -out ${out_db}                        \
                                      2> ${out_db}.log  1>&2;
                                           
                         date +"DB created on %Y/%m/%d %T %Z %s"       \
                            > ${cdfile};
                         echo "CCC" >> ${params.wdir}/kk;
                     fi;
              fi;
         else
             echo "${refseqs_fasta} is empty... Cannot create a blastdb\n" > ${out_db}.log
             echo "DDD" >> ${params.wdir}/kk;
        fi;
        echo "EEE" >> ${params.wdir}/kk;
        """
}


process do_blastn {

    input:
        val query // with whole path
        val rdb // with whole path
        val out_dir // output directory
        
    output:
        val "${blast_dir}/${blast_algn}"

    script:

       spid=query.toString().split('/')[-1].split('[.]')[0]
       blast_q=query.toString().split('/')[-1].replaceAll(".gz", "").replaceAll(".fa", "")
       blast_r=rdb.toString().split('/')[-1]
       blast_algn="${blast_q}_ON_${blast_r}.${params.blast_approach}.tbl"
       blast_dir="${out_dir}/${spid}"
       dbindex="${rdb}.nin";
       """
       echo "BLASTN $query  :: $blast_q \n $rdb :: $blast_r \n ${blast_dir}/${blast_algn}" >> ${params.wdir}/kk.$spid;
        if [ -d $blast_dir ]; then 
                echo "$blast_dir"; 
            else 
                mkdir $blast_dir; 
        fi;

        if [[ -e ${dbindex} ]]; then
            if [[ ${query} =~ .gz\$ ]]; then
              gzip -dc ${query} |\
                 blastn -query  -  -db ${rdb}                     \
                 -num_alignments 1 -perc_identity 50              \
                  -task blastn -out ${blast_dir}/${blast_algn}    \
                  -evalue 10e-10    -dbsize ${params.taxondbsize}      \
                  -outfmt \"${params.bl_outfmt}\"    \
                  -num_threads ${params.NCPUS}       \
                   2> ${blast_dir}/${blast_algn}.log 1>&2;
            else
              blastn -query ${query} -db ${rdb}          \
                  -num_alignments 1 -perc_identity 50             \
                  -task blastn -out ${blast_dir}/${blast_algn}    \
                  -evalue 10e-10    -dbsize ${params.taxondbsize}      \
                  -outfmt \"${params.bl_outfmt}\"    \
                  -num_threads ${params.NCPUS}       \
                   2> ${blast_dir}/${blast_algn}.log 1>&2;
            fi;
        else
            echo "No database found for blast." \
                 > ${blast_dir}/${blast_algn}.log;
        fi;
       """

}

process do_blast_kaiju {

    input:
        val query // with whole path
        val out_dir // output directory
        val dummy
        
    output:
        val "${blast_dir}/${blast_algn}"

    script:

       spid=query.toString().split('/')[-1].split('[.]')[0]
       blast_q=query.toString().split('/')[-1].replaceAll(".gz", "").replaceAll(".fa", "")
       blast_r="${blast_q}.reference_blastdb" // rdb.toString().split('/')[-1]
       db=query.toString().replaceAll(".gz", "").replaceAll(".fa", "")
       rdb="${db}.reference_blastdb"
       blast_algn="${blast_q}_ON_${blast_r}.${params.blast_approach}.tbl"
       blast_dir="${out_dir}/${spid}"
       dbindex="${rdb}.nin";
       """
       echo "BLASTN $query  :: $blast_q \n $rdb :: $blast_r \n ${blast_dir}/${blast_algn}" >> ${params.wdir}/kk.$spid;
       echo "$dummy";
        if [ -d $blast_dir ]; then 
                echo "$blast_dir"; 
            else 
                mkdir $blast_dir; 
        fi;

        if [[ -e ${dbindex} ]]; then
            if [[ ${query} =~ .gz\$ ]]; then
              gzip -dc ${query} |\
                 tblastx -query  -  -db ${rdb}                      \
                    -num_alignments 1               \
                    -out ${blast_dir}/${blast_algn}    \
                    -dbsize ${params.taxondbsize}                     \
                    -outfmt \"${params.bl_outfmt}\"    \
                    -num_threads ${params.NCPUS}       \
                    2> ${blast_dir}/${blast_algn}.log 1>&2;
            else
              tblastx -query ${query} -db ${rdb}          \
                  -num_alignments 1               \
                  -out ${blast_dir}/${blast_algn}    \
                  -dbsize ${params.taxondbsize}                     \
                  -outfmt \"${params.bl_outfmt}\"    \
                  -num_threads ${params.NCPUS}       \
                   2> ${blast_dir}/${blast_algn}.log 1>&2;
            fi;
        else
            echo "No database found for blast." \
                 > ${blast_dir}/${blast_algn}.log;
        fi;
       """

}


process do_tblastx {

    input:
        val query // with whole path
        val db // with whole path
        
    output:
        val "${blast_dir}/${blast_algn}"
        
    script:

       blast_q=query.toString().split('/')[-1].replaceAll(/.gz$/, "").replaceAll(/.fa$/, "")
       blast_r=db.toString().split('/')[-1]
       blast_algn="${blast_q}_ON_${blast_r}.${params.blast_approach}.tbl"
       blast_dir=params.contigs_blast_dir
       
       """
        if [[ \$(echo ${query}) =~ .gz\$ ]]; then
          gzip -dc ${query} |\
             tblastx -query  -  -db ${db}          \
              -out ${blast_dir}/${blast_algn}    \
              -num_alignments 1 -perc_identity 50              \
              -evalue 10e-10   -subject_besthit           \
              -outfmt \"${params.bl_outfmt}\"    \
              -num_threads ${params.NCPUS}       \
               2> ${blast_dir}/${blast_algn}.log;
        else
          tblastx -query ${query} -db ${db}          \
              -out ${blast_dir}/${blast_algn}    \
              -num_alignments 1 -perc_identity 50              \
              -evalue 10e-10      -subject_besthit       \
              -outfmt \"${params.bl_outfmt}\"    \
              -num_threads ${params.NCPUS}       \
               2> ${blast_dir}/${blast_algn}.log;
        fi;
               
       """

}


process best_reciprocal_hit {
    input:
        
        val(q_fa)  //query fasta ( contigs )
        val(r_fa)  //reference fasta 
        val(qor_blast)  // query onto reference blast (whole path)
        val(roq_blast)  // reference onto query blast (whole path)
    
    output:
        val "${out_brh}.rbbh.tbl"
        
    exec:
    
    out_brh=qor_blast.replaceAll(/.tbl$/, "")
    qor_fl= new File("${qor_blast}")
    roq_fl= new File("${roq_blast}")
    if( qor_fl.size() > 0 && roq_fl.size() > 0 ) {

        """
            python2 ${BIND}/rbbh.py       \
                ${qor_blast}   ${q_fa}        \
                ${roq_blast}   ${r_fa}        \
                1e-10    0                    \
              > ${out_brh}.rbbh.tbl    \
             2> ${out_brh}.rbbh.log;
        """
    }else{
        brh_fl= new File("${out_brh}.rbbh.tbl")
        brh_fl.write ""
    
    }
}

process merge_blast_outs {
    input: 
        val(bl1)
        val(bl2)
    
    output:
        val("${blast_all}"), emit: BALL
    
    script:
        splist=bl1.toString().split('[.]')
        sproot=splist[0..3].join(".")
        blast_all="${sproot}.clas+unc_blast.${params.blast_approach}.tbl"
        spid=splist[-1]
    """
    if [ -f $bl1 -a -f $bl2 ] ; then
        cat ${bl1} ${bl2} >  ${blast_all} 2> ${blast_all}.log;
        echo "FFF ${blast_all}" >> ${params.wdir}/kk.$spid;
    fi;
    echo "XXX $spid" >> ${params.wdir}/kk.$spid;
    """

}

process blast_sum_coverage {
    input:
        val(blastout)
        val(clids)
        val(unids)
        
    output:
        val("$blast_sumcov"), emit: TBL
      //  val("$unaling_ids"),emit: UNIDS
      
    script:
    blast_sumcov=blastout.replaceAll(".tbl",".coverage.tbl")
    coverage_log=blastout.replaceAll(".tbl",".coverage.tbl.log")
    
    blast_merge=blast_sumcov.replaceAll(".tbl",".merge")
    merge_log=blast_sumcov.replaceAll(".tbl",".merge.tbl.log")
    
    sampid=blast_sumcov.split('/')[-1].split('[.]')[0]
    // odir="${params.taxfastdir}/${sampid}"
   //  unaling_ids=allqids.replaceAll(".ids", ".unaligned.ids")
    """
    if [ -s $blastout ]; then 
    ## Get summary of blast out coverage:
        ${params.bindir}/coverage_blastshorttbl.pl \
              $blastout > $blast_sumcov 2> $coverage_log;
        
        ## merge:
        if ($clids==F); then 
                ${params.bindir}/virwaste_taxon_output_merge_blastonly.pl   \
                  -i ${params.blast_refseqs_dir}/C-RVDB_allentries_info.txt       \
                  -B $blast_sumcov    \
                  --mincov ${params.mincovpct} -o ${blast_merge}  2> ${merge_log} 1>&2;
        else
                ${params.bindir}/virwaste_taxon_output_merge_blastonly.pl   \
                  -i ${params.blast_refseqs_dir}/C-RVDB_allentries_info.txt       \
                  -B $blast_sumcov   -c ${clids}  -u ${unids}     \
                  --mincov ${params.mincovpct}  -o ${blast_merge}  2> ${merge_log} 1>&2;
        fi;
    fi;
    """
    // ## Get list of unclassified ids:
    // ## gawk 'BEGIN{ while(getline<ARGV[1]>0) G[\$1]=\$1; ARGV[1]=""; } 
    // ##       \$2 in G {} else { F[\$1]++; } 
    // ##       END{ for (f in F) print f; }
    // ##       ' $blastout $allqids > $unaling_ids;
    // ## grep -v  -f 
    

}

