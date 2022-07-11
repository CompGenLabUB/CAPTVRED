#! /usr/bin/env nextflow
nextflow.enable.dsl=2

params.wdir = "$BDIR"
params.runID = "$RUNID"
params.sampletbl = "$params.wdir/samples_definition.tbl"
params.rawfq = "$params.wdir/rawseqs_fastq/*R{1,2}*.fastq.gz"

params.NCPUS = 32

params.rawfq_dir="$params.wdir/rawseqs_fastq"
params.reports_dir="$params.wdir/reports"
params.clnfq_dir="$params.wdir/cleanseqs"


// include { load_sampleids; samplecheck } from './check_inputs.nf'
include { samplecheck } from './check_inputs.nf'

//include { initvars } from './settings.nf'

include { bbduk_clean2 } from './rawfq_clean2.nf'
include { fastQC } from './seq_stats.nf'

def samplesMap = [:]

workflow init_samples() {

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
    
        
    println "# Running fastqc on raw fastq sequences ... ..."
    fastQC(Channel.fromPath(params.rawfq))
    println "# Done"
    
}



workflow {

    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Starting  : $workflow.userName $ZERO $workflow.start"
    println "# Reading samples for $params.runID from $params.sampletbl"
    init_samples()
			//println "# Summary of samples to run:"
    
			//paths_list is a LoL: where first item is the raw fastq root for ech sample and the second onw is the new rood for the cleanseqs.
    
   
    println "# Cleaning raw fastq sequences with BBDuk ... ..."
    def paths_list=[]
    samplesMap.each { sampleID, illuminaID ->
         def newsamp=["$params.rawfq_dir/$illuminaID", "$params.clnfq_dir/$sampleID"]
         paths_list << newsamp
    
    }
			// paths_list.each {
			// println ">> $it" // `it` is an implicit parameter corresponding to the current element
			// }
			// println "Total number of samples: $paths_list.size"

    bbduk_clean2(Channel.from(paths_list)) 
    println "# Done"
   
   println "# Running fastqc on raw fastq sequences ... ..."
   fastQC( bbduk_clean2.out.mix() )
   println "# Done"
    
    println "# .........................."
    println "Reads cleaning with BBDuk ..."
    
    
    println "# Finishing : $workflow.userName $workflow.complete"
    
}
