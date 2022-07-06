#! /usr/bin/env nextflow


//ch_tbl=Channel.fromPath("${samptbl}", type: 'file')
//Channel
  //.fromPath( "${samptbl}", checkIfExists: true )  		// if empty, complains						                           
  //.set {ch_tbl} 						// make the channel "reads"

process load_sampids {
	input:
	 file x 
	output:
	 stdout
	
	exec:
	SamplesDef = file(x)
	//#SamplesDef = file("${ch_tbl}")
	allLines = SamplesDef.readLines()   //#Do NOT use on large files!
	def sampMap = [:]
	for( line in allLines ) {
		def samp = line.split('\t')
		// ignore lines stating with "#"
		if (!samp[0].startsWith("#")) {
			sampMap.(samp[0]) = (samp[1])
		}
	}
//	println "Sample IDs loaded:"
//	sampMap.each { key, value ->
//		   println "IlluminaID: $key \t SampleID: $value"
// }
}



//ch_fqdir = Channel.fromPath("${rawfq_dir}", type: 'dir')
/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */


process samplecheck {
	println "Entered process samplecheck!!"
    input:
     tuple val(pair_id), path(reads)
     
    output:
     stdout
    
    
    shell:
      '''
      echo "It is availabe!" # $pair_id \t $reads "
      '''
}
