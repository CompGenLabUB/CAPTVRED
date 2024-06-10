#! /usr/bin/env nextflow

include { db_for_kaiju; init_rvdb; merge_rvdb_setref; create_logf } from './modules/init_conf.nf'

workflow prepareblast (){
    take:
      db
      logfl
    main:
      //For blast with kaiju:
      db_for_kaiju(db, logfl)
    
    emit:
      ""
} 

// // // // // // MAIN // // // // // //  
    
    
    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Starting  : $workflow.userName  $workflow.start"
  
workflow {
    def refseqs = "${workflow.projectDir}/references"
    def refseqs_rvdb="$refseqs/${params.rvdb_dir}"
     // create_logf("$workflow.projectDir/references")
        create_logf(refseqs)
     // prepare rvdb + set 
      init_rvdb(refseqs_rvdb, params.rvdb_link, params.db_name, create_logf.out)
      merge_rvdb_setref(refseqs_rvdb, init_rvdb.out,params.db_name, params.set_seqs, create_logf.out )


     // Dbs non-viral discard:
     // db_for_kaiju(params.filtDB, create_logf.out)
    
     // prepareblast(params.blast_kaiju_db, create_logf.out)
}