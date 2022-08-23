#! /usr/bin/env nextflow

process kaiju_contigs {
    label 'limit_kaiju'
    
    input:
        tuple   val(db_id), val(contigs_fa)

    output:
        //Potser no cal output aqui
        val "${odir}/${samp_id}_all.kaiju.${db_id}.names.out" // ?????
    
    script:
    
    //samp_id=contigs_fa.split('/')[-1].replaceAll(".contigs+singletons.fa","")
    samp_id=contigs_fa.split('/')[-1].split('[.]')[0]
    odir="${params.taxdir}/kaiju/${samp_id}/contigs_taxon"

    """
         [ -d ${odir} ] || mkdir -vp ${odir}
         kaiju -v -z ${task.cpus}                                    \
                  -t  ${params.kaijuDBD}/${db_id}/nodes.dmp         \
                  -f ${params.kaijuDBD}/kaiju_db_${db_id}.fmi       \
                  -i  ${contigs_fa}                                 \
                  -o ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out;
                  
          kaiju-addTaxonNames -t ${params.kaijuDBD}/${db_id}/nodes.dmp   \
                -n ${params.kaijuDBD}/${db_id}/names.dmp                 \
                -r superkingdom,genus,species                            \
                -i ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out                 \
                -o ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.names.out;
        
        kaiju2table -t ${params.kaijuDBD}/${db_id}/nodes.dmp              \
                    -n ${params.kaijuDBD}/${db_id}/names.dmp              \
                    -r species                                            \
                    -o ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.summary.tsv      \
                    ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out;
        
        kaiju2krona -u -v                                              \
                    -t ${params.kaijuDBD}/${db_id}/nodes.dmp           \
                    -n ${params.kaijuDBD}/${db_id}/names.dmp           \
                    -i ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out           \
                    -o ${odir}/${samp_id}.contigs+sgl.kaiju.${db_id}.out.krona;
        
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
    odir="${params.taxdir}/kaiju/${samp_id}/reads_taxon"
        
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
                -r superkingdom,species                                    \
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
                -r superkingdom,species                                  \
                -i ${odir}/${samp_id}_sg.kaiju.${db_id}.out              \
                -o ${odir}/${samp_id}_sg.kaiju.${db_id}.names.out;
                  
        #ALL
           #summary table
             kaiju2table -t ${params.kaijuDBD}/${db_id}/nodes.dmp         \
                    -n ${params.kaijuDBD}/${db_id}/names.dmp              \
                    -r species                                            \
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
     
     when:
       txn_names_tbl =~ /nr_euk/
     
    script:

     celids=txn_names_tbl.replaceAll(".names.out", ".nonviral.ids")
     samp_id=txn_names_tbl.split('/')[-1].split('[.]')[0].replaceAll("_all", "")
     
        """
        awk ' (\$1 ~ /C/)  &&  !(\$8 ~ /Viruses|NA/) { print \$2}'    \
          ${txn_names_tbl}   >  ${celids};

        seqkit grep --invert-match --pattern-file  ${celids}        \
                    ${params.clnfq_dir}/${samp_id}_pe1.fastq.gz    |\
                    gzip -   > ${params.clnfq_dir}/${samp_id}_pe1.filtered.fastq.gz;
                    
        seqkit grep --invert-match --pattern-file ${celids}         \
                    ${params.clnfq_dir}/${samp_id}_pe2.fastq.gz    |\
                    gzip -   > ${params.clnfq_dir}/${samp_id}_pe2.filtered.fastq.gz;
                    
        seqkit grep --invert-match --pattern-file  ${celids}       \
                    ${params.clnfq_dir}/${samp_id}_sgl.fastq.gz   |\
                    gzip -   > ${params.clnfq_dir}/${samp_id}_sgl.filtered.fastq.gz;
        """



}
