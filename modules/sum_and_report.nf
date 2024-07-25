#! /usr/bin/env nextflow

process make_summary_tbl {
    input:
      val sumfls_list
      val asbdir
      val repdir
      val bindir
      val samplestb

    output:
      val ("$outfinal")
    
    script:
 
    // samplestb="${params.samp}";
    conrp="${repdir}/${params.runID}_contigs_counts.tbl";
    // if (params.taxonslow == true){taxfl="${params.taxslowdir}"};
    // if (params.taxonfast == true){ taxfl="${params.taxfastdir}"};
    taxrp="${repdir}/${params.runID}_taxonomy_stats.tbl";
    outfinal="${repdir}/${params.runID}_samples_final_summary.tbl"
    sumflsstr = sumfls_list.join(" ")
    """
   ## Summarize info of contigs and blast:
    echo  "#SAMPLE_ID\tN_CONTIGS\tN_CON+SGL" > ${conrp}
    echo  "#SAMPLE_ID\tN_READS\tN_SEQS\tNSPECS" > ${taxrp}
    assem=\$(echo ${params.assembler} | awk '{print tolower(\$0)}' )
    for i in \$(awk '\$1!~"^#"{print \$1}' $samplestb ); do 
         if [ \$assem == "megahit" ]; 
         then
             ncns=\$(grep -c '^>' "${asbdir}/\${i}/\${i}.contigs.fa");
         elif [ \$assem == "metaspades" ]; 
         then
             ncns=\$(grep -c '^>' "${asbdir}/\${i}/scaffolds.fasta");
         fi;
         
         ncsg=\$(grep -c '^>' "${asbdir}/\${i}/\${i}.contigs+singletons.fa");
         echo \${i}\t\${ncns}\t\${ncsg} >> ${conrp}
    
    done;
    cat  ${sumflsstr} >> $taxrp; 
    # python script
    python3 ${bindir}/virwaste_get_summary_table.py                        \
          ${samplestb}                                                     \
          ${repdir}/${params.runID}_multiqc_raw_data/multiqc_fastqc.txt    \
          ${repdir}/${params.runID}_multiqc_clean_data/multiqc_fastqc.txt  \
          ${repdir}/${params.runID}_multiqc_filt_data/multiqc_fastqc.txt   \
          ${repdir}/${params.runID}_multiqc_align_data/multiqc_bowtie2.txt \
          ${conrp}       \
          ${taxrp}       \
          ${outfinal}    \
          ${params.R1}     ${params.R2} 
     """
}

process fill_html_report {

    input:
       val samples_sum
       val samps
       val bindir
       val figsdir
       val repdir
       val htmldir

    output:
      val final_html

    script:
      
      final_html="${repdir}/${params.runID}_CAPTVRED_final_report.html"
    
      """
       python3 ${bindir}/CAPTVRED_create_report.py                \
            ${samps}                                          \
            ${repdir}                                        \
            ${params.runID}                                              \
            ${samples_sum}                                               \
            ${htmldir}/CAPTVRED_final_report_template.html  \
            ${final_html}     ${params.taxalg}
      """
}







