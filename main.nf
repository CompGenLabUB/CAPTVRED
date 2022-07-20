#! /usr/bin/env nextflow
// nextflow.enable.dsl=2

params.wdir = "$BDIR"
params.runID = "$RUNID"
params.sampletbl = "$params.wdir/samples_definition.tbl"
params.rawfq = "$params.wdir/rawseqs_fastq/*R{1,2}*.fastq.gz"
params.amplicon_refseqs="$AMPSQ"
params.amplicon_refseqs_dir="$AMPSQD"

params.NCPUS = 32

params.rawfq_dir="$params.wdir/rawseqs_fastq"
params.reports_dir="$params.wdir/reports"
params.clnfq_dir="$params.wdir/cleanseqs"
params.ampaln_dir="$params.wdir/amplicons_alignment"
params.asbl_dir="$params.wdir/assembly"

// include { load_sampleids; samplecheck } from './check_inputs.nf'
include { samplecheck } from './check_inputs.nf'

//include { initvars } from './settings.nf'

include { bbduk_clean } from './rawfq_clean.nf'
include { fastQC; multiQC_raw; multiQC_clean; multiQC_bowtie_amp } from './seq_stats.nf'
include { generate_index_bowtie; bowtie_amplicons_alignment; bowtie_amplicons_alignment_sg } from './reads_align.nf'
include { megahit_assembly_pe; megahit_assembly_sg; megahit_assembly_all} from './reads_assembly.nf'
include { trinity_assembly_sg; trinity_assembly_pe } from './reads_assembly.nf'
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
 
 // def logFile = new File()
 
 
   
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

    bbduk_clean.out.outPE1
    bbduk_clean.out.outPE2
    bbduk_clean.out.outSGL

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

workflow amplicon_sequences_bowtie() {

    take:
      x

    main:
      bowtie_amplicons_alignment(x) | collect

    emit:
      bowtie_amplicons_alignment.out

}

workflow amplicon_sequences_align() {

   take:
     x

   main:
 
    def cluster_index_path = params.amplicon_refseqs.toString().replaceAll(".fa.gz", "_index")
   
    def samples_lst=[]
    def thysamples = samplesMap
    thysamples.each { sampleID, illuminaID ->
         def newsamp=["${cluster_index_path}",
                      "${sampleID}_pe1.fastq.gz",
                      "${sampleID}_pe2.fastq.gz",
                      "${sampleID}_sgl.fastq.gz"]
         samples_lst << newsamp
    }

   
    bowtie_amplicons_alignment(x, Channel.from(samples_lst))
    bowtie_amplicons_alignment_sg(x, Channel.from(samples_lst))
    
    multiQC_bowtie_amp(bowtie_amplicons_alignment.out.mix(bowtie_amplicons_alignment_sg.out).collect())
  
  emit:

    bowtie_amplicons_alignment.out.mix(bowtie_amplicons_alignment_sg.out)

}


workflow assembly_sequences_flow() {

   take:
    p_end1
    p_end2
    s_end

   main:

     // MEGAHIT //
    //megahit_assembly_pe(p_end1, p_end2)
    //megahit_assembly_sg(s_end)
    megahit_assembly_all(p_end1, p_end2, s_end)
    
    // TRINNITY //
    trinity_assembly_sg(s_end)
    trinity_assembly_pe(p_end1, p_end2)
    
  emit:
    //megahit_assembly_pe.out
    //megahit_assembly_sg.out
    megahit_assembly_all.out
   

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
    generate_index_bowtie(reads_clean.out.mix().collect(),params.amplicon_refseqs)
    amplicon_sequences_align(generate_index_bowtie.out.collect())
    assembly_sequences_flow(reads_clean.out )
    
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser...\n" : "Oops .. something went wrong: ${workflow.errorMessage}" )
}
