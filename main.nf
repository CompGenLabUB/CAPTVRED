#! /usr/bin/env nextflow
// nextflow.enable.dsl=2

params.wdir = "$BDIR"
params.runID = "$RUNID"
params.sampletbl = "$params.wdir/samples_definition.tbl"
params.rawfq = "$params.wdir/rawseqs_fastq/*R{1,2}*.fastq.gz"
params.amplicon_refseqs_dir="$AMPSQD"
params.amplicon_refseqs="$AMPSQFA"
params.contigs_blast_dir="$CBLASTD"
params.kaijuDBD="$KAIDBD"
params.kaijudbs=["viruses", "rvdb", "refseq", "nr_euk"]

params.NCPUS = 32

params.rawfq_dir="$params.wdir/rawseqs_fastq"
params.reports_dir="$params.wdir/reports"
params.clnfq_dir="$params.wdir/cleanseqs"
params.ampaln_dir="$params.wdir/amplicons_alignment"
params.asbl_dir="$params.wdir/assembly"
params.taxdir="$params.wdir/taxonomy"


// blast parameters //
params.bl_outfmt="6 qseqid qlen sseqid slen qstart qend sstart send length score evalue bitscore pident nident mismatch positive gapopen gaps ppos qframe sframe qcovs qcovhsp qseq sseq";

////
params.assembler="MEGAHIT"
params.blast_approach="BLASTN"

// include { load_sampleids; samplecheck } from './check_inputs.nf'
include { samplecheck } from './check_inputs.nf'

//include { initvars } from './settings.nf'

include { bbduk_clean } from './rawfq_clean.nf'
include { fastQC; multiQC_raw; multiQC_clean; multiQC_bowtie_amp } from './seq_stats.nf'
include { generate_index_bowtie; bowtie_amplicons_alignment; bowtie_amplicons_alignment_sg } from './reads_align.nf'
include { megahit_assembly_all} from './reads_assembly.nf'
include { trinity_assembly_sg; trinity_assembly_pe } from './reads_assembly.nf'
include { make_db_for_blast; do_blastn; do_tblastx; best_reciprocal_hit } from './contigs_align.nf'
include { kaiju_raw; discard_celular } from './taxonomy.nf'
def samplesMap = [:]
SamplesDef = file(params.sampletbl) //(samplestbl_file)
SamplesDef.eachLine {
    line -> {
        def samp = line.split('\t')
        // ignore lines stating with "#"
        if (!samp[0].startsWith("#")) {
            samplesMap.(samp[0]) = (samp[1])
        }
    }
}

log.info """\
 =======================================
 V I R W A S T E - N F   P I P E L I N E
 =======================================
 RUN     : ${params.runID}
 Samples : ${samplesMap}
 SysInfo : ${workflow.userName} SID=${workflow.sessionId} NCPUs=${params.NCPUS} GITcid=${workflow.commitId}
 =======================================
 """
 

 
   
workflow init_samples() {

  main:
    println "# INIT: $samplesMap"

  emit:
    ""

}

workflow fastqc_onrawseqs() {

  take:
    x

  main:
   fastQC(Channel.fromPath(params.rawfq)) | collect | multiQC_raw
    
}


workflow reads_clean() {

  take:

    x

  main:
    
    
    def paths_list=[]
         //paths_list is a LoL: where first item is the raw fastq root 
         //for each sample and the second onw is the new root for the cleanseqs.
    def thysamples = samplesMap
    thysamples.each { sampleID, illuminaID ->
         def newsamp=["$params.rawfq_dir/$illuminaID", "$params.clnfq_dir/$sampleID"]
         paths_list << newsamp
    }
 
    bbduk_clean(x, Channel.from(paths_list)) 

    fastQC( bbduk_clean.out.mix() ) | collect | multiQC_clean

    
  emit:
   PE1=bbduk_clean.out.outPE1
   PE2=bbduk_clean.out.outPE2
   SGL=bbduk_clean.out.outSGL

}

workflow reads_filter_celular() {
    take:
      x
     
     main:
        kaiju_raw(x)
        discard_celular(kaiju_raw.out)
     

    emit:
       discard_celular.out.PE1out
       discard_celular.out.PE2out
       discard_celular.out.SGLout
}


workflow amplicon_sequences_dbinit() {

    take:
      x
      refseqs

    main:
      generate_index_bowtie(x,refseqs)
      
    emit:
      generate_index_bowtie.out


}


workflow amplicon_sequences_align() {

   take:
     cluster_index_path
     pe1
     pe2
     sgl
    
    
   main:
   
    bowtie_amplicons_alignment(cluster_index_path, pe1, pe2, sgl)
    bowtie_amplicons_alignment_sg(cluster_index_path, pe1, pe2, sgl )
    
    multiQC_bowtie_amp(bowtie_amplicons_alignment.out.mix(bowtie_amplicons_alignment_sg.out).collect())
  
  emit:

    bowtie_amplicons_alignment.out.mix(bowtie_amplicons_alignment_sg.out)
    
}

workflow megahit_assembly_flow () {
    take: 
      pe1
      pe2
      sgl      
    
    main:
      ref_fa="${params.amplicon_refseqs}"
      megahit_assembly_all(pe1,pe2, sgl)
      blast_flow( ref_fa, megahit_assembly_all.out.CGSout)
      blast_flow_rev( megahit_assembly_all.out.CGSout, ref_fa)
      if (params.blast_approach ==~ /(?i)blastn/) {
        QOR=blast_flow.out[0]
        ROQ=blast_flow_rev.out[0]
      }else if(params.blast_approach ==~ /(?i)tblastx/ ){
        QOR=blast_flow.out[1]
        ROQ=blast_flow_rev.out[1]
      }
      best_reciprocal_hit(megahit_assembly_all.out.CGSout, ref_fa, QOR, ROQ )

    emit:
      best_reciprocal_hit.out
}


workflow trinity_assembly_flow () {
    take: 
      pe1
      pe2
      sgl      
    
    main:
      ref_fa="${params.amplicon_refseqs}"
      trinity_assembly_pe(pe1, pe2, sgl)
      //trinity_assembly_sg(sgl)
      blast_flow(ref_fa, trinity_assembly_pe.out)
      blast_flow_rev(trinity_assembly_pe.out, ref_fa)
      if (params.blast_approach ==~ /(?i)blastn/) {
        QOR=blast_flow.out[0]
        ROQ=blast_flow_rev.out[0]
      }else if(params.blast_approach ==~ /(?i)tblastx/ ){
        QOR=blast_flow.out[1]
        ROQ=blast_flow_rev.out[1]
      }
      best_reciprocal_hit(trinity_assembly_pe.out, ref_fa, QOR, ROQ )

    emit:
      best_reciprocal_hit.out
}


workflow blast_flow() {

   take:
     ref_fasta
     query_fasta

   main:
    
    make_db_for_blast( ref_fasta)
    do_blastn(query_fasta, make_db_for_blast.out)
    do_tblastx(query_fasta, make_db_for_blast.out)
    
  emit:
    do_blastn.out
    do_tblastx.out
}

workflow blast_flow_rev() {

   take:
     ref_fasta
     query_fasta
     

   main:

    make_db_for_blast(ref_fasta)
    do_blastn(query_fasta, make_db_for_blast.out)
    do_tblastx(query_fasta, make_db_for_blast.out)
  emit:
    do_blastn.out
    do_tblastx.out
}

// // // // // // MAIN // // // // // //  

workflow {

    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Starting  : $workflow.userName $ZERO $workflow.start"
    println "# Reading samples for $params.runID from $params.sampletbl"
    
    
    init_samples()
    fastqc_onrawseqs(init_samples.out)
    reads_clean(init_samples.out)
    
    KDB=Channel.from(params.kaijudbs)
    CLNR=reads_clean.out.merge()
    to_kaiju=KDB.combine(CLNR).merge()
    reads_filter_celular(to_kaiju)
    
    generate_index_bowtie(params.amplicon_refseqs)
    amplicon_sequences_align(generate_index_bowtie.out, reads_filter_celular.out)

    if (params.assembler ==~ /(?i)MEGAHIT/){
         megahit_assembly_flow(reads_filter_celular.out )
    }else if (params.assembler =~ /(?i)TRINITY/ ){
        trinity_assembly_flow(reads_filter_celular.out)
    }

    // workflow for coverage ploting (inputs will be: bowtie out and assembly out)
}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser...\n" : "Oops .. something went wrong: ${workflow.errorMessage}" )
}
