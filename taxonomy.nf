#! /usr/bin/env nextflow

process kaiju_contigs {
    label 'limit_kaiju'
    errorStrategy 'finish'
    
    input:
        tuple   val(db_id), val(contigs_fa)

    output:

        val("$wnames")
    
    script:
    
    samp_id=contigs_fa.split('/')[-1].split('[.]')[0]
    odir="${params.taxfastdir}/${samp_id}"
    wnames="${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.names.out"
    
   
    """
         ## Find taxonomy classification:
         [ -d ${odir} ] || mkdir -vp ${odir}
         kaiju -v -z ${task.cpus}                                   \
                  -t  ${params.kaijuDBD}/${db_id}/nodes.dmp         \
                  -f ${params.kaijuDBD}/kaiju_db_${db_id}.fmi       \
                  -i  ${contigs_fa}                                 \
                  -E 0.001    -s 65 -e 2                          \
                  -o ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out \
                  2> ${params.logs_dir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out
          
          ## Add species names
          kaiju-addTaxonNames -t ${params.kaijuDBD}/${db_id}/nodes.dmp   \
                -n ${params.kaijuDBD}/${db_id}/names.dmp                 \
                -r superkingdom,family,species                            \
                -i ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out                 \
                -o ${wnames};
        
        ## transform to tabular file.
        kaiju2table -t ${params.kaijuDBD}/${db_id}/nodes.dmp                          \
                    -n ${params.kaijuDBD}/${db_id}/names.dmp                          \
                    -r species                                                        \
                    -l superkingdom,phylum,class,order,family,genus,species           \
                    -o ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.summary.tsv      \
                    ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out;
        
        ## Create kronaplots
        kaiju2krona -u -v                                              \
                    -t ${params.kaijuDBD}/${db_id}/nodes.dmp           \
                    -n ${params.kaijuDBD}/${db_id}/names.dmp           \
                    -i ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out           \
                    -o ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out.krona;
        
        ktImportText -o ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out.krona.html \
                        ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out.krona
         echo "KAIJU $wnames " >> ${params.wdir}/kk.$samp_id;
    """

}


process extract_ids () {
    input:
        val(wnames)

    output: 
        // tuple  val("$txn_cl"), val("$outf_un"), val("$outf_cl")
         val("$txn_cl"), emit: CLT
         val("$outf_un"), emit: UNIDS
         val("$outf_cl"), emit: CLIDS
        
    script:
        outf_cl=wnames.replaceAll(".names.out", ".classified.ids");  //IDs of classified contigs
        facl=outf_cl.replaceAll(".ids", ".fa");
        outf_un=wnames.replaceAll(".names.out", ".unclassified.ids"); //IDs of unclassified contigs
        txn_cl=wnames.replaceAll(".names.out", ".classified.taxon" );  //TaxonIDs of classified contigs
        spid=wnames.toString().split('/')[-1].split('[.]')[0];
        """
       
        ## Get id list of classified and unclassified reads:
        awk 'BEGIN{FS="\\t|;"} {if (\$1 ~ /^C\$/ && \$8 ~ "Viruses" && \$(NF-1) !~ "NA") print \$2}'  $wnames > $outf_cl;
        awk 'BEGIN{FS="\\t|;"} {if ((\$1 ~ /^C\$/ && (\$8 !~ "Viruses" || \$(NF-1) ~ "NA")) || \$1 ~ /^U\$/) print \$2}'  $wnames > $outf_un;
         
        #taxonids of classified
        awk 'BEGIN{ FS="\\t|;" }
             { if ((\$1 ~ /C/) && (\$8 ~ "Viruses") && !(\$(NF-1) ~ "NA" )) print \$3 }
            ' $wnames | sort -n | uniq > $txn_cl;
        
        echo "EXTRACTIDS $txn_cl   \n $outf_un  \n $outf_cl" >> ${params.wdir}/kk.$spid;
        """

}

process taxonid_to_fasta{
    input:
        val(cl_txn)
        
    output:
        val("$ref_fa"), emit: "REF"
        
    script:

        ref_fa=cl_txn.replaceAll(/.taxon$/, ".reference.fa")
        cl_refid=cl_txn.replaceAll(/.taxon$/, ".reference.fa.ids")
        sampid=cl_txn.split('/')[-1].split('[.]')[0]
        odir="${params.taxfastdir}/${sampid}"
        DB_FA="${params.blast_refseqs_dir}/C-RVDBvCurrent.fasta.gz"
        ids_tbl="${params.blast_refseqs_dir}/C-RVDB_genebank_to_taxon.ids"
        
        """
         echo "TAXONID TO FASTA $cl_txn   \n $ref_fa" > ${params.wdir}/kk.$sampid;
        [ -d $odir ] || mkdir $odir;
        touch ${ref_fa};
        if [[ -s $DB_FA  &&  -s $cl_txn ]];
        then 
           
           gawk 'BEGIN{ while(getline<ARGV[1]>0) G[\$1]=\$1; ARGV[1]=""; } \$2 in G { F[\$1]++; } END{ for (f in F) print f; }' $cl_txn $ids_tbl > $cl_refid  \
                       2> ${params.logs_dir}/${sampid}.taxonid_to_fasta.gawk.log;
           seqkit grep -r -m 0 --pattern-file  ${cl_refid}        \
                    ${DB_FA}  > ${ref_fa} 2> ${params.logs_dir}/${sampid}.taxonid_to_fasta.seqkit.log;
        
        else
            echo " $DB_FA or $cl_txn is empty" >>  ${params.logs_dir}/${sampid}.taxonid_to_fasta.gawk.log;
        
        fi
        
        """
}

process readid_to_fasta{

    input:

        val(ids)
        
    output:
        val("$fasta"), emit: "FA"
        
    script:
        fasta=ids.replaceAll(".ids", ".fa");
        sampid=ids.split('/')[-1].split('[.]')[0]
        odir="${params.taxfastdir}/${sampid}"
        CONTIGS_FA="${params.asbl_dir}/megahit/${sampid}/${sampid}.contigs+singletons.fa"
        
        """
        echo "READID TO FASTA $ids   \n $fasta" >> ${params.wdir}/kk.$sampid;
        [ -d $odir ] || mkdir $odir;
        touch ${fasta};

       ## obtain fasta of the "unclassified" set of contigs
              #seqkit 
        if [ -s $CONTIGS_FA ] && [ -s  $ids ]; then 
             seqkit grep --pattern-file  ${ids}       \
                    ${CONTIGS_FA}  > ${fasta}            \
                    2> ${params.logs_dir}/seqkit_unclass.${sampid}.log;
        else
            touch ${fasta};
        fi;
        
        """
}

process kaiju_raw {
    label 'limit_kaiju'
    
    input:
        tuple  val(db_id), val(pe1),  val(pe2), val(sgle)
        // dbs: only name.  //reads: whole path

                
    output:
        val "${odir}/${samp_id}_all.kaiju.${db_id}.names.out"
    
    script:
    
    samp_id=sgle.split('/')[-1].replaceAll("_sgl.fastq.gz", "")
    odir="${params.taxdir}/reads_taxon/${samp_id}"
        
        """
            [ -e ${odir} ] ||  mkdir -vp ${odir};
            #PE
            kaiju -v -z ${task.cpus}                           \
                 -t  ${params.kaijuDBD}/${db_id}/nodes.dmp        \
                 -f  ${params.kaijuDBD}/kaiju_db_${db_id}.fmi     \
                 -i ${pe1}     -j ${pe2}                          \
                 -o  ${odir}/${samp_id}_pe.kaiju.${db_id}.out;

            kaiju-addTaxonNames -t ${params.kaijuDBD}/${db_id}/nodes.dmp   \
                -n ${params.kaijuDBD}/${db_id}/names.dmp                   \
                -r superkingdom,genus,species                              \
                -i ${odir}/${samp_id}_pe.kaiju.${db_id}.out                \
                -o ${odir}/${samp_id}_pe.kaiju.${db_id}.names.out;
            
            #SG
            kaiju -v -z ${task.cpus}                               \
                  -t  ${params.kaijuDBD}/${db_id}/nodes.dmp         \
                  -f ${params.kaijuDBD}/kaiju_db_${db_id}.fmi       \
                  -i  ${sgle}                                       \
                  -o ${odir}/${samp_id}_sg.kaiju.${db_id}.out;
                  
          kaiju-addTaxonNames -t ${params.kaijuDBD}/${db_id}/nodes.dmp   \
                -n ${params.kaijuDBD}/${db_id}/names.dmp                 \
                -r superkingdom,phylum,class,order,family,genus,species  \
                -i ${odir}/${samp_id}_sg.kaiju.${db_id}.out              \
                -o ${odir}/${samp_id}_sg.kaiju.${db_id}.names.out;
                  
        #ALL
           #summary table
             kaiju2table -t ${params.kaijuDBD}/${db_id}/nodes.dmp         \
                    -n ${params.kaijuDBD}/${db_id}/names.dmp              \
                    -r species   \
                    -o ${odir}/${samp_id}_all.kaiju.${db_id}.summary.tsv  \
                    ${odir}/${samp_id}_pe.kaiju.${db_id}.out              \
                    ${odir}/${samp_id}_sg.kaiju.${db_id}.out;

            cat  ${odir}/${samp_id}_pe.kaiju.${db_id}.out            \
                 ${odir}/${samp_id}_sg.kaiju.${db_id}.out            \
               > ${odir}/${samp_id}_all.kaiju.${db_id}.out;

            cat  ${odir}/${samp_id}_pe.kaiju.${db_id}.names.out      \
                 ${odir}/${samp_id}_sg.kaiju.${db_id}.names.out      \
               > ${odir}/${samp_id}_all.kaiju.${db_id}.names.out;


            kaiju2krona -u -v                                          \
                    -t ${params.kaijuDBD}/${db_id}/nodes.dmp           \
                    -n ${params.kaijuDBD}/${db_id}/names.dmp           \
                    -i ${odir}/${samp_id}_all.kaiju.${db_id}.out       \
                    -o ${odir}/${samp_id}_all.kaiju.${db_id}.out.krona;

            ktImportText -o ${odir}/${samp_id}_all.kaiju.${db_id}.out.krona.html \
                        ${odir}/${samp_id}_all.kaiju.${db_id}.out.krona;

        """

}



process discard_nonviral {
    input:
     val txn_names_tbl
    
    output:
     val "${params.clnfq_dir}/${samp_id}_pe1.filtered.fastq.gz", emit: PE1out
     val "${params.clnfq_dir}/${samp_id}_pe2.filtered.fastq.gz", emit: PE2out
     val "${params.clnfq_dir}/${samp_id}_sgl.filtered.fastq.gz", emit: SGLout
     
     //when:
       //txn_names_tbl =~ /nr_euk/
     
    script:

     nvids=txn_names_tbl.replaceAll(".names.out", ".nonviral.ids")
     samp_id=txn_names_tbl.split('/')[-1].split('[.]')[0].replaceAll("_all", "")
     
        """
        awk ' (\$1 ~ /C/)  &&  !(\$8 ~ /Viruses|NA/) { print \$2}'    \
          ${txn_names_tbl}   >  ${nvids};

        seqkit grep --invert-match --pattern-file  ${nvids}        \
                    ${params.clnfq_dir}/${samp_id}_pe1.fastq.gz    |\
                    gzip -   > ${params.clnfq_dir}/${samp_id}_pe1.filtered.fastq.gz;
                    
        seqkit grep --invert-match --pattern-file ${nvids}         \
                    ${params.clnfq_dir}/${samp_id}_pe2.fastq.gz    |\
                    gzip -   > ${params.clnfq_dir}/${samp_id}_pe2.filtered.fastq.gz;
                    
        seqkit grep --invert-match --pattern-file  ${nvids}       \
                    ${params.clnfq_dir}/${samp_id}_sgl.fastq.gz   |\
                    gzip -   > ${params.clnfq_dir}/${samp_id}_sgl.filtered.fastq.gz;
        """



}
