#! /usr/bin/env nextflow

process coverage_plots {
    input:
        val pebam
        val sgbam
        val blastout
        val covfl
        val bindir
        val reptdir
        val datab_dir
        val gffdir
    
    output:
        val outdir, emit: ODIR
    
    script:
    
    gffblast=blastout.replaceAll(/.tbl$/ ,".gff")
    samp_id=pebam.split('/')[-1].toString().replaceAll("_pe.bowtie.sorted.mapped.bam","")
    outdir="${reptdir}/coverage_figures"
    
    infofl="${datab_dir}/info_summary.tsv"
    coordfl="${gffdir}/refseqs_coordinates.tsv"

    """
     mkdir -vp ${outdir};

     if [ -s ${covfl} ]; then 
     bash ${bindir}/cov-blast2gff.sh  ${covfl}   ${blastout}   ${gffblast}
     else
          touch  ${gffblast};
     fi;
        
     Rscript ${bindir}/CAPTVRED_coverage_figures.R              \
             ${gffdir}          ${coordfl}                        \
             ${samp_id}       ${pebam}         ${sgbam}         \
             ${gffblast}      ${infofl}                         \
             ${outdir}                                          \
            2> ${reptdir}/coverage_figures/Coverage_${samp_id}.log 1>&2;

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
