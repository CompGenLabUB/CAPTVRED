#! /usr/bin/env nextflow


process make_db_for_blast {
    input:
        val(refseqs_fasta)
        val(redodb)

    output:

        val "${out_db}", emit: DB
        val "${is_db}", emit: CTRL

    script:
    
        out_db=refseqs_fasta.toString().replaceAll(".fasta.gz|.fasta|.fa.gz|.fa", "_blastdb")
        dblist=out_db.split('/')
        db_name=dblist[-1]
        cdfile=[dblist[0..dblist.size()-2].join('/'), "blastdb_created_genome.cdate"].join('/')
        spid=refseqs_fasta.toString().split('/')[-1].split('[.]')[0]
       """
        if [ -s ${refseqs_fasta} ]; then
              if [ -e ${cdfile} -a ${redodb} == "FALSE" ]; then
                 cat ${cdfile} 1>&2;
                    
              else
                     if [[ \$(echo ${refseqs_fasta}) =~ .gz\$ ]]; then
                         gunzip -c ${refseqs_fasta}               | \
                           makeblastdb -in -  -dbtype nucl         \
                                       -title \"${db_name}\"       \
                                       -out ${out_db}              \
                                       2> ${out_db}.log  1>&2;
                                       
                           date +"DB created on %Y/%m/%d %T %Z %s"    \
                              > ${cdfile};
                     else
                         makeblastdb -in ${refseqs_fasta}  -dbtype nucl    \
                                      -title \"${db_name}\"                 \
                                      -out ${out_db}                        \
                                      2> ${out_db}.log  1>&2;
                                           
                         date +"DB created on %Y/%m/%d %T %Z %s"       \
                            > ${cdfile};
                     fi;
              fi;
         else
             echo "${refseqs_fasta} is empty... Cannot create a blastdb\n" > ${out_db}.log
             
        fi;
        """
}


process do_blastn {

    input:
        val query // with whole path
        val rdb // with whole path
        val out_dir // output directory
        
    output:
        val "${blast_dir}/${blast_algn}", emit:OUT
        val "${query}", emit: CONTIGS
    script:

       spid=query.toString().split('/')[-1].split('[.]')[0]
       blast_q=query.toString().split('/')[-1].replaceAll(".gz", "").replaceAll(".fa", "")
       blast_r=rdb.toString().split('/')[-1]
       blast_algn="${blast_q}_ON_${blast_r}.BLASTN.tbl"
       blast_dir="${out_dir}/${spid}"
       dbindex="${rdb}.nin";
       """

        if [ -d $blast_dir ]; then 
                echo "$blast_dir"; 
            else 
                mkdir $blast_dir; 
        fi;

        if [[ -e ${dbindex} ]]; then
            if [[ ${query} =~ .gz\$ ]]; then
              gzip -dc ${query} |\
                 blastn -query  -  -db ${rdb}         -dust no                   \
                 -num_alignments 1 -perc_identity   ${params.blast_pident}       \
                  -task blastn -out ${blast_dir}/${blast_algn}                   \
                  -evalue ${params.blast_eval}   -dbsize ${params.taxondbsize}   \
                  -outfmt \"${params.bl_outfmt}\"    \
                  -num_threads ${params.NCPUS}       \
                   2> ${blast_dir}/${blast_algn}.log 1>&2;
            else
              blastn -query ${query} -db ${rdb}        -dust no                       \
                  -num_alignments 1 -perc_identity ${params.blast_pident}             \
                  -task blastn -out ${blast_dir}/${blast_algn}                        \
                  -evalue  ${params.blast_eval}     -dbsize ${params.taxondbsize}     \
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
        val rdb // with whole path
        val out_dir // output directory
        
    output:
        val "${blast_dir}/${blast_algn}", emit:OUT
        val "${query}", emit: CONTIGS
        
    script:

       spid=query.toString().split('/')[-1].split('[.]')[0]
       blast_q=query.toString().split('/')[-1].replaceAll(/.gz$/, "").replaceAll(/.fa$/, "")
       blast_r=rdb.toString().split('/')[-1]
       blast_algn="${blast_q}_ON_${blast_r}.TBLASTX.tbl"
       blast_dir="${out_dir}/${spid}" //params.contigs_blast_dir
       dbindex="${rdb}.nin";
       
       """

        if [ -d $blast_dir ]; then 
                echo "$blast_dir"; 
            else 
                mkdir $blast_dir; 
        fi;

        
        if [[ -e ${dbindex} ]]; then
            if [[ \$(echo ${query}) =~ .gz\$ ]]; then
                gzip -dc ${query} | \
                    tblastx -query  -  -db ${rdb}                    \
                    -out ${blast_dir}/${blast_algn}                  \
                    -num_alignments 1                                \
                    -evalue ${params.blast_eval}                     \
                    -subject_besthit                                 \
                    -outfmt \"${params.bl_outfmt}\"                  \
                    -num_threads ${params.NCPUS}                     \
                    2> ${blast_dir}/${blast_algn}.log;
            else
                tblastx -query ${query} -db ${rdb}           \
                    -out ${blast_dir}/${blast_algn}          \
                    -num_alignments 1                        \
                    -evalue ${params.blast_eval}             \
                    -subject_besthit                         \
                    -outfmt \"${params.bl_outfmt}\"          \
                    -num_threads ${params.NCPUS}             \
                    2> ${blast_dir}/${blast_algn}.log;
            fi;
        else
            echo "No database found for blast." \
                 > ${blast_dir}/${blast_algn}.log;
        fi
               
       """

}

process blast_sum_coverage {
    input:
        val(blastout)
        val(clids)
        val(unids)
        
    output:
        val("$blast_sumcov"), emit: TBL
        val(byrd), emit: BYR
        val(bysq), emit: BYSQ
        val(bysp), emit: BYSP
        val(stats), emit: SUM
        val(stats), emit: SUM2
      
    script:
    blast_sumcov=blastout.replaceAll(/\.tbl/,".coverage.tbl")
    coverage_log=blastout.replaceAll(/\.tbl/,".coverage.tbl.log")
    
    blast_merge=blastout.replaceAll(/\.tbl/,".merge")
    merge_log=blast_sumcov.replaceAll(/\.tbl/,".merge.tbl.log")

    
    if (params.taxalg ==~ /(?)BLASTN/){blalg="blastn"};
    if (params.taxalg ==~ /(?)TBLASTX/){ blalg="tblastx"};

    
    sampid=blast_sumcov.split('/')[-1].split('[.]')[0]
    
     //reports dir
     byrd="${params.reports_dir}/${sampid}.${blalg}_taxonomysum_byread.tbl"
     bysq="${params.reports_dir}/${sampid}.${blalg}_taxonomysum_bysequence.tbl"
     bysp="${params.reports_dir}/${sampid}.${blalg}_taxonomysum_byspecie.tbl"
     stats="${params.reports_dir}/${sampid}.${blalg}.stats.out"

    """
    if [ -s $blastout ]; then 
    ## Get summary of blast out coverage:
        ${params.bindir}/coverage_blastshorttbl.pl \
              $blastout > $blast_sumcov 2> $coverage_log;
        
        ## merge:
        if ($clids==F);   then 
                ${params.bindir}/virwaste_taxon_output_merge_blastonly.pl         \
                  -i ${params.blast_refseqs_dir}/${params.blast_ref_db_info}       \
                  -B $blast_sumcov   --mincov ${params.mincovpct}                 \
                   --pe ${params.asbl_dir}/${params.assembler}/${sampid}/${sampid}_pe.bowtie_onto_contigs.maped.sorted.stats         \
                   --sg ${params.asbl_dir}/${params.assembler}/${sampid}/${sampid}_se.bowtie_onto_contigs.maped.sorted.stats         \
                  --samp ${sampid}   -o ${blast_merge}  2> ${merge_log} 1>&2;
        else
                ${params.bindir}/virwaste_taxon_output_merge_blastonly.pl                                  \
                  -i ${params.blast_refseqs_dir}/C-RVDB_allentries_info.txt                                \
                  -B $blast_sumcov   --mincov ${params.mincovpct}                                          \
                  --pe ${params.asbl_dir}/${params.assembler}/${sampid}/${sampid}_pe.bowtie_onto_contigs.maped.sorted.stats      \
                  --sg ${params.asbl_dir}/${params.assembler}/${sampid}/${sampid}_se.bowtie_onto_contigs.maped.sorted.stats      \
                   -c ${clids}  -u ${unids}                                                                \
                   -o ${blast_merge}  2> ${merge_log} 1>&2;
        fi;
    fi;

     cp  ${blast_merge}_taxonomysum_byread.tbl      ${byrd}
     cp  ${blast_merge}_taxonomysum_bysequence.tbl  ${bysq}
     cp  ${blast_merge}_taxonomysum_byspecie.tbl    ${bysp}
     cp  ${blast_merge}.stats.out                   ${stats}
    
    """ 
}

process do_cov_onrefseqs () {  
       // Get FASTA FILES (refseqs and contigs) for coverage of refseqs ans perform blast analysis.
    
    input:
        val(RefSqsDef)   // Refseqs info file (tsv)
        val(refsfa)      // Refseqs sequences fasta file
        val(contigsfa)   //  contigs fasta file
        val(assign)      // assignations file taxonomy_byread.tbl
        val(out_suffix)  // suffix of the final blast out file 

    output:
       val("$finalblout"), emit: BLOUT
       val("$blast_sumcov"), emit: COV

    script:
       spid=contigsfa.toString().split('/')[-1].split('[.]')[0]
       finalblout="${params.taxvircan}/${spid}/${spid}_${out_suffix}.tbl"
       blast_sumcov="${params.taxvircan}/${spid}/${spid}_${out_suffix}.coverage.tbl"
       coverage_log="${params.taxvircan}/${spid}/${spid}_${out_suffix}.coverage.log"
       
      """
            mkdir -p ${params.taxvircan}/${spid};
            [ -e $finalblout ] && rm -v $finalblout;
            touch $finalblout;
            for txn in \$(awk -vIFS='\t' '\$1!~/^#/ {print \$3}' $RefSqsDef | sort | uniq);
             do {
               while read ln;
                 do {
                   awk -vFS="\t" -vtaxid=\$txn '\$3==taxid{print \$6}' $RefSqsDef  > ${params.tmp_dir}/${spid}.\${txn}.sqids;
                 }; done < $RefSqsDef;
               
                # 0. init files
                [ -e ${params.tmp_dir}/${spid}.\${txn}.reference.fa  ] && rm -v ${params.tmp_dir}/${spid}.\${txn}.reference.fa;
                [ -e ${params.tmp_dir}/${spid}.\${txn}.contigs.ids   ] && rm -v ${params.tmp_dir}/${spid}.\${txn}.contigs.ids;
                [ -e ${params.tmp_dir}/${spid}.\${txn}.contigs.fa    ] && rm -v ${params.tmp_dir}/${spid}.\${txn}.contigs.fa;
                [ -e ${params.tmp_dir}/${spid}.\${txn}.blast_out.tbl ] && rm -v ${params.tmp_dir}/${spid}.\${txn}.blast_out.tbl;
                
                 # 1. # Extract fastas of refseqs:
                 zcat ${refsfa} | seqkit grep -f ${params.tmp_dir}/${spid}.\${txn}.sqids - \
                              > ${params.tmp_dir}/${spid}.\${txn}.reference.fa;
                      
                 # 2. # Get contigs ids (filter by taxonid). And extract fasta of contigs of interest.
                 awk -vtaxid="\${txn}" '\$9==taxid {print \$1}                       \
                        ' ${assign} > ${params.tmp_dir}/${spid}.\${txn}.contigs.ids;
                 if [ -s ${params.tmp_dir}/${spid}.\${txn}.contigs.ids ]; 
                       then 
                         seqkit grep -f ${params.tmp_dir}/${spid}.\${txn}.contigs.ids ${contigsfa} \
                                      > ${params.tmp_dir}/${spid}.\${txn}.contigs.fa;
                         #echo "2 \$txn   :  DONE:)" >> /data/virpand/pandemies/TEST_SET/DEV_TESTSET/pv.${spid};

                         # 3. # Do the blast and add results to "final blast" file.
                         egrep -Hc '^>' ${params.tmp_dir}/${spid}.\${txn}.reference.fa     \
                                        ${params.tmp_dir}/${spid}.\${txn}.contigs.fa       \
                                     >> /data/virpand/pandemies/TEST_SET/SIMDATA/pv.${spid};
                          blastn -query   ${params.tmp_dir}/${spid}.\${txn}.contigs.fa           \
                                 -subject ${params.tmp_dir}/${spid}.\${txn}.reference.fa         \
                                 -num_alignments 1                                               \
                                 -evalue ${params.vcan_eval}                                      \
                                 -outfmt \"${params.bl_outfmt}\"                                 \
                                 -out ${params.tmp_dir}/${spid}.\${txn}.blast_out.tbl            \
                                 2> ${params.tmp_dir}/${spid}.\${txn}.blast_out.log;
                           if [ -s ${params.tmp_dir}/${spid}.\${txn}.blast_out.tbl  ]; 
                               then 
                                   cat ${params.tmp_dir}/${spid}.\${txn}.blast_out.tbl >> ${finalblout};
                                   printf "%d HITS found\\n" \$(wc -l ${params.tmp_dir}/${spid}.\${txn}.blast_out.tbl | awk '{print \$1}') >> /data/virpand/pandemies/TEST_SET/DEV_TESTSET/pv.${spid};
                               fi;
                      fi;
           }; done;

          if [ -s ${finalblout} ]; then 
           ## Get summary of blast out coverage:
           ${params.bindir}/coverage_blastshorttbl.pl \
              ${finalblout} > $blast_sumcov 2> $coverage_log;
          else
            echo "BLAST MERGED FILE IS EMPTY "  > $blast_sumcov;
          fi;

      """
}


process do_cov_on_viralcandidates () { 
     input:
        val(RefSqsDef)   // Refseqs info file (tsv)
        val(vircandb)      // Refseqs sequences fasta file
        val(contigsfa)   //  contigs fasta file
        val(assign)      // assignations file taxonomy_byread.tbl
        val(out_suffix)  // suffix of the final blast out file 

    output:
       val("$finalblout"), emit: BLOUT
       val("$blast_sumcov"), emit: COV

    script:
       spid=contigsfa.toString().split('/')[-1].split('[.]')[0]
       finalblout="${params.taxvircan}/${spid}/${spid}_${out_suffix}.tbl"
       blast_sumcov="${params.taxvircan}/${spid}/${spid}_${out_suffix}.coverage.tbl"
       coverage_log="${params.taxvircan}/${spid}/${spid}_${out_suffix}.coverage.log"
        """
            mkdir -p ${params.taxvircan}/${spid};
            [ -e $finalblout ] && rm -v $finalblout;
            touch $finalblout;
            for txn in \$(awk -vIFS='\t' '\$1!~/^#/ {print \$3}' $RefSqsDef | sort | uniq);
             do {
                  # 0. # Init files:
                        [ -e ${params.tmp_dir}/${spid}.\${txn}.contigs.ids   ] && rm -v ${params.tmp_dir}/${spid}.\${txn}.contigs.ids;
                        [ -e ${params.tmp_dir}/${spid}.\${txn}.contigs.fa    ] && rm -v ${params.tmp_dir}/${spid}.\${txn}.contigs.fa;
                        
                  # 2. # Get contigs ids (filter by taxonid). And extraxt fasta of contigs of interest.
                        awk -vtaxid="\${txn}" '\$9==taxid {print \$1} \
                          ' ${assign} > ${params.tmp_dir}/${spid}.\${txn}.contigs.ids;
                   if [ -s ${params.tmp_dir}/${spid}.\${txn}.contigs.ids ]; 
                       then 
                         seqkit grep -f ${params.tmp_dir}/${spid}.\${txn}.contigs.ids ${contigsfa} \
                                      > ${params.tmp_dir}/${spid}.\${txn}.contigs.fa;
                         cat ${params.tmp_dir}/${spid}.\${txn}.contigs.fa >> ${params.taxvircan}/${spid}/all_class_contigs.fa;
                    fi;
             }; done;
                     # 3. # Do the blast
                     blastn -query   ${params.taxvircan}/${spid}/all_class_contigs.fa            \
                            -db  ${vircandb}                                                     \
                                 -num_alignments 1       -perc_identity ${params.vcan_pident}    \
                                 -evalue ${params.vcan_eval}                                     \
                                 -outfmt \"${params.bl_outfmt}\"                                 \
                                 -out  ${finalblout}                                             \
                                   2> ${finalblout}.log;
                                   
             if [ -s ${finalblout} ]; then 
                 ## Get summary of blast out coverage:
                 ${params.bindir}/coverage_blastshorttbl.pl \
                   ${finalblout} > $blast_sumcov 2> $coverage_log;
            
             fi;
        """
 }

