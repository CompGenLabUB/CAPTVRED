#! /usr/bin/env nextflow

process make_db_for_blast {
    input:
        val(refseqs_fasta)

    output:

        val "${out_db}", emit: DB
        val "${is_db}", emit: CTRL

    exec:
        out_db=refseqs_fasta.toString().replaceAll(".fa.gz|.fa", "_blastdb")
        dblist=out_db.split('/')
        db_name=dblist[-1]
        fa_fl= new File("${refseqs_fasta}")
        if( fa_fl.size() > 0 ) {
            is_db=1
            cdfile=[dblist[0..dblist.size()-2].join('/'), "blastdb_created_genome.cdate"].join('/')
            cdfile2 = new File("$cdfile");
            if( cdfile2.exists() ) {
                """
                cat ${cdfile} 1>&2;
                """
                
            } else {

                if  (refseqs_fasta =~ /\.fa\.gz$/ ) {
                    """
                     gunzip -c ${refseqs_fasta}                |\
                        makeblastdb -in -  -dbtype nucl         \
                                    -title ${db_name}           \
                                    -out ${out_db}              \
                                  2> ${out_db}.log;
                                  
                     date +"DB created on %Y/%m/%d %T %Z %s"    \
                          > ${cdfile};
                    """
                
                }
                
                if  (refseqs_fasta =~ /\.fa$/ ) {
                    """
                     makeblastdb -in ${refseqs_fasta}  -dbtype nucl    \
                                 -title ${db_name}                     \
                                 -out ${out_db}                        \
                                 2> ${out_db}.log;
                                      
                     date +"DB created on %Y/%m/%d %T %Z %s"       \
                          > ${cdfile};
                    """
                }
           }
    }else{
     is_db=0
     log_fl= new File("${out_db}.log")
     log_fl.write "${refseqs_fasta} is empty... Cannot create a blastdb\n"  
    
    }

}


process do_blastn {

    input:
        val query // with whole path
        val db // with whole path
        
    output:
        val "${blast_dir}/${blast_algn}"
        
    script:

       blast_q=query.toString().split('/')[-1].replaceAll(".gz", "").replaceAll(".fa", "")
       blast_r=db.toString().split('/')[-1]
       blast_algn="${blast_q}_ON_${blast_r}.${params.blast_approach}.tbl"
       blast_dir=params.contigs_blast_dir
       
       """
        if [[ \$(echo ${query}) =~ .gz\$ ]]; then
          gzip -dc ${query} |\
             blastn -query  -  -db ${db}          \
              -task blastn -out ${blast_dir}/${blast_algn}    \
              -evalue 10e-10     -subject_besthit   \
              -outfmt \"${params.bl_outfmt}\"    \
              -num_threads ${params.NCPUS}       \
               2> ${blast_dir}/${blast_algn}.log;
        else
          blastn -query ${query} -db ${db}          \
              -task blastn -out ${blast_dir}/${blast_algn}    \
              -evalue 10e-10    -subject_besthit  \
              -outfmt \"${params.bl_outfmt}\"    \
              -num_threads ${params.NCPUS}       \
               2> ${blast_dir}/${blast_algn}.log;
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
              -evalue 10e-10   -subject_besthit           \
              -outfmt \"${params.bl_outfmt}\"    \
              -num_threads ${params.NCPUS}       \
               2> ${blast_dir}/${blast_algn}.log;
        else
          tblastx -query ${query} -db ${db}          \
              -out ${blast_dir}/${blast_algn}    \
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


process blast_sum_coverage {
    input:
        val(blastout)
        val(allqids)
        
    output:
        val("$blast_sumcov"), emit: TBL
        val("$unaling_ids"),emit: UNIDS
    exec:
    blast_sumcov=blastout.replaceAll(".tbl",".coverage.tbl")
    unaling_ids=allqids.replaceAll(".ids", ".unaligned.ids")
    
    """
    ## Get summary of blast out coverage:
    ${params.bindir}/coverage_blastshorttbl.pl \
                                $blastout > $blast_sumcov
                                
    ## Get list of unclassified ids:
    gawk 'BEGIN{ while(getline<ARGV[1]>0) G[\$1]=\$1; ARGV[1]=""; } 
          \$2 in G {} else { F[\$1]++; } 
          END{ for (f in F) print f; }
          ' $blastout $allqids > $unaling_ids;
    #grep -v  -f 
    """

}

