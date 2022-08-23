#! /usr/bin/env nextflow

process make_db_for_blast {
    input:
        val(refseqs_fasta)

    output:

        val "${out_db}"

    script:
        
        out_db=refseqs_fasta.toString().replaceAll(".fa.gz|.fa", "_blastdb")
        dblist=out_db.split('/')
        db_name=dblist[-1]
        cdfile=[dblist[0..dblist.size()-2].join('/'), "blastdb_created_genome.cdate"].join('/')
        cdfile2 = new File("$cdfile");
        if( cdfile2.exists() ) {
            """
            cat ${cdfile} 1>&2;
            
            """
            
        } else {
            if  (refseqs_fasta =~ /\.fa\.gz$/ ) 
                """
                 gunzip -c ${refseqs_fasta}                |\
                    makeblastdb -in -  -dbtype nucl         \
                                -title ${db_name}           \
                                -out ${out_db}              \
                              2> ${out_db}.log;
                              
                 date +"DB created on %Y/%m/%d %T %Z %s"    \
                      > ${cdfile};
                """
                
            if  (refseqs_fasta =~ /\.fa$/ ) 
                """
                 makeblastdb -in ${refseqs_fasta}  -dbtype nucl            \
                             -title ${db_name}              \
                             -out ${out_db} \
                             2> ${out_db}.log;
                                  
                 date +"DB created on %Y/%m/%d %T %Z %s"       \
                      > ${cdfile};
                    """
            
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
       blast_algn="${blast_q}_ON_${blast_r}.tbl"
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
       blast_algn="${blast_q}_ON_${blast_r}.tblastx.tbl"
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
        val "${out_brh}.blastn.rbbh.tbl"
        
    script:
    
    out_brh=qor_blast.replaceAll(/.tbl$/, "")
    
    """
        python2 ${BIND}/rbbh.py       \
            ${qor_blast}   ${q_fa}        \
            ${roq_blast}   ${r_fa}        \
            1e-25    80                    \
          > ${out_brh}.rbbh.tbl    \
         2> ${out_brh}.rbbh.log;
        
         
    """
}




