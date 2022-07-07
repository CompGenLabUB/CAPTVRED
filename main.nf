#! /usr/bin/env nextflow
nextflow.enable.dsl=2

params.wdir = "$BDIR"
params.runID = "$RUNID"
params.sampletbl = "$params.wdir/samples_definition.tbl"

//params.rawfq_dir='$launchDir/rawseqs_fastq'
params.rawfq = "$params.wdir/rawseqs_fastq/*_{R1,R2}_001.fastq.gz"

// include { load_sampleids; samplecheck } from './check_inputs.nf'
include { samplecheck } from './check_inputs.nf'

//include { initvars } from './settings.nf'

//include { bbduk; fastqc } from './cleanreads.nf'

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
     
}
    

workflow {

    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Starting  : $workflow.userName $ZERO $workflow.start"

    println "# Reading samples for $params.runID from $params.sampletbl"    
    init_samples()
    println "# Summary of samples to run:"
    samplesMap.each { key, value ->
         println "Sample: $key   Fastq root: $value"
    }
    println "# .........................."

    
   //  Channel
   //      .fromFilePairs( params.rawfq)
   //      .ifEmpty { error "Cannot find any reads matching: ${params.rawfq}" }
   //      .set { read_pairs_ch } 
   //  samplecheck(read_pairs_ch )  // by now it does nothing 

    println "# Finishing : $workflow.userName $workflow.complete"
    
}
