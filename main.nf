#! /usr/bin/env nextflow
nextflow.enable.dsl=2

params.samptbl='$launchDir/samples_definition.tbl'
//params.rawfq_dir='$launchDir/rawseqs_fastq'
params.rawfq="$launchDir/rawseqs_fastq/*_{R1,R2}_001.fastq.gz"

include {  load_sampids; samplecheck } from './check_inputs.nf'

//include { initvars } from './settings.nf'

//include { bbduk; fastqc } from './cleanreads.nf'


workflow  {
	println "entered flow!!"
	ch_tbl=Channel.fromPath(params.samptbl, type: 'file')
	load_sampids(ch_tbl)
	Channel
		.fromFilePairs( params.rawfq)
		.ifEmpty { error "Cannot find any reads matching: ${params.rawfq}" }
		.set { read_pairs_ch } 
	samplecheck(read_pairs_ch )  // by now it does nothing 
}
