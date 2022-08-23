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
    idtoname="${params.amplicon_refseqs_dir}/Viral_candidates_zoonosis_refseqs_idrelation.tsv"
    
    """
    source ${params.bindir}/rbbh-blast2gff.sh    \
        ${brh}                            \
        ${blastout}                       \
        ${gffblast};
    
    
    Rscript ${params.bindir}/virwaste_coverage_figures.R         \
         ${rdir}          ${rdir}/refseqs_coordinates.tbl        \
         ${samp_id}       ${pebam}         ${sgbam}              \
         ${gffblast}      ${idtoname}                            \
         ${params.reports_dir}/coverage_figures                  \
        2> ${params.reports_dir}/coverage_figures/Coverage_${samp_id}.log;
    """

}
