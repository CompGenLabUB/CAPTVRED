#! /usr/bin/env nextflow

process coverage_plots {
    input:
        val pebam
        val sgbam
        val blastout
        val brh
        
    //output? no cal
    
    script:
    
    gffblast=blastout.replaceAll(/.tbl$/ ,".gff")
    rdir="${params.amplicon_refseqs_dir}/gff_refgenomes"
    samp_id=pebam.split('/')[-1].toString().replaceAll("_pe.bowtie.sorted.mapped.bam","")
    //idtoname="${params.amplicon_refseqs_dir}/Viral_candidates_zoonosis_refseqs_idrelation.tsv"
    infofl="${params.amplicon_refseqs_dir}/refseqs_ids_relation.tsv"
    
    """
    bash ${params.bindir}/rbbh-blast2gff.sh  ${brh}   ${blastout}   ${gffblast}

    
    Rscript ${params.bindir}/virwaste_coverage_figures.R         \
         ${rdir}          ${rdir}/refseqs_coordinates.tbl        \
         ${samp_id}       ${pebam}         ${sgbam}              \
         ${gffblast}      ${infofl}                              \
         ${params.reports_dir}/coverage_figures                  \
        2> ${params.reports_dir}/coverage_figures/Coverage_${samp_id}.log;
    """

}


process align_counts_plot {
    input:
        val bamlist
    
    script:
    
    bamfiles=bamlist.join(" ")
    tblfl="${params.ampaln_dir}/${params.runID}_bowtie_info.tbl"
    infofl="${params.amplicon_refseqs_dir}/refseqs_ids_relation.tsv"
    figpref="${params.reports_dir}/${params.runID}"
    
    """
    [ -e $tblfl ] && rm $tblfl;
    for i in ${bamfiles}; do \
        sp=\$(echo \$i | rev | cut -d '/' -f1 | rev | cut -d '.' -f1); \
        samtools view \$i |\
          awk -v spid=\$sp 'BEGIN {OFS="\t"}{ print \$1, \$3, \$5, \$9, spid}'  \
              >> $tblfl;\
     done;
     
     awk 'BEGIN{IFS="\t"; OFS="\t"}{print \$4,\$2,\$7 }'       \
          ${params.amplicon_refseqs_dir}/${params.amplicon_refseqs_info}  \
          >  ${infofl};
     
     Rscript ${params.bindir}/aligned_reads_on_amplicons.R     \
             ${tblfl} ${infofl} ${figpref}
    """
}
