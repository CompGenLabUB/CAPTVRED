#! /usr/bin/env nextflow

process make_summary_tbl {
    input:
      val x
    output:
      val ("$outfinal")
    
    script:
    just_for_workflow_control=x
    samplestb="${params.sampletbl}";
    conrp="${params.reports_dir}/${params.runID}_contigs_counts.tbl";
    if (params.taxonslow == true){taxfl="${params.taxslowdir}"};
    if (params.taxonfast == true){ taxfl="${params.taxfastdir}"};
    taxrp="${params.reports_dir}/${params.runID}_taxonomy_stats.tbl";
    outfinal="${params.reports_dir}/${params.runID}_samples_final_summary.tbl"
    """
   echo $just_for_workflow_control;
   ## Summarize info of contigs and blast:
    echo  "#SAMPLE_ID\tN_CONTIGS\tN_CON+SGL" > ${conrp}
    echo  "#SAMPLE_ID\tN_READS\tN_SEQS\tNSPECS" > ${taxrp}
    assem=\$(echo ${params.assembler} | awk '{print tolower(\$0)}' )
    for i in \$(awk '\$1!~"^#"{print \$1}' $samplestb ); do 
         if [ \$assem == "megahit" ]; 
         then
             ncns=\$(grep -c '^>' "${params.asbl_dir}/${params.assembler}/\${i}/\${i}.contigs.fa");
         elif [ \$assem == "metaspades" ]; 
         then
             ncns=\$(grep -c '^>' "${params.asbl_dir}/${params.assembler}/\${i}/scaffolds.fasta");
         fi;
         ncsg=\$(grep -c '^>' "${params.asbl_dir}/${params.assembler}/\${i}/\${i}.contigs+singletons.fa");
         echo \${i}\t\${ncns}\t\${ncsg} >> ${conrp}
         
         cat  ${taxfl}/\${i}/\${i}.contigs+singletons_*stats.out >> $taxrp; 
    
    done;

    #python script
    python3 ${params.bindir}/virwaste_get_summary_table.py        \
          ${samplestb}                                    \
          ${params.reports_dir}/${params.runID}_multiqc_raw_data/multiqc_fastqc.txt    \
          ${params.reports_dir}/${params.runID}_multiqc_clean_data/multiqc_fastqc.txt  \
          ${params.reports_dir}/${params.runID}_multiqc_filt_data/multiqc_fastqc.txt   \
          ${params.reports_dir}/${params.runID}_multiqc_bowtie_amplicons_alignment_data/multiqc_bowtie2.txt  \
          ${conrp}  \
          ${taxrp}  \
          ${outfinal}
     """
}

process fill_html_report {

    input:
       val samples_sum

    output:
      val final_html

    script:
      
      final_html="${params.reports_dir}/${params.runID}_CAPTVRED_final_report.html"
    
      """
       python3 ${params.bindir}/CAPTVRED_create_report.py                \
            ${params.sampletbl}                                          \
            ${params.reports_dir}                                        \
            ${params.runID}                                              \
            ${samples_sum}                                               \
            ${params.nfscripts_dir}/CAPTVRED_final_report_template.html  \
            ${final_html}  
      """
}







