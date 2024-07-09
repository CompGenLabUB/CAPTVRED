#! /usr/bin/env nextflow

include { db_for_kaiju;   create_logf; stdrd_link; stdrd_link as stdrd_link_otfm } from './modules/init_conf.nf'
include { merge_rvdb_setref; chek_setref_ids; filter_FOI } from './modules/filter_rvdb.nf'
include { get_taxonids; get_taxonids_rvdb} from './modules/db_taxonomy.nf'
include { taxonomizator; taxonomizator as taxonomizator_rvdb} from './modules/db_taxonomy.nf'
include { get_rvdb; get_names_and_nodes } from './modules/update_files.nf' 


workflow check_files () {
   take:
      database
      merged
      subset
   
   main:

      if (params.dbsplit_update==false) {
         if(!new File(subset).exists()){
            params.dbsplit_update=true
            }
      }

      if (params.merge_update==false) {
         if (!new File(merged).exists()){
            params.merge_update=true
            params.dbsplit_update=true
            }
      }
      
      if (params.rvdb_update==false) {
         if (!new File(database).exists()){
            params.rvdb_update=true
            params.merge_update=true
            params.dbsplit_update=true
            }
      }

    println "# RVDB  upd  : $params.rvdb_update"
    println "# Merge upd  : $params.merge_update"
    println "# Split upd  : $params.dbsplit_update"

   emit:
      "FILES CHECKED!"
}

workflow database (){
   
   take:
      dir
      link
      name
      logf
      check

   main:

   check.view()
   get_rvdb(dir, link, name, logf)
   

   emit:
      get_rvdb.out

}

workflow mergefasta () {

   take:
      refseqs_ncbi
      refseqs_rvdb
      db_fa
      set_fa
      merged_fa
      link
      logf

   main:
      get_names_and_nodes(refseqs_ncbi, link, logf)
      chek_setref_ids(refseqs_rvdb, db_fa, set_fa, params.setname)
      merge_rvdb_setref(merged_fa, db_fa, params.db_name, chek_setref_ids.out, logf )
   

   emit:
      MERGED=merge_rvdb_setref.out
      NCBID=get_names_and_nodes.out

}


workflow database_subset (){
   take:
      rvdbdir
      ncbidir
      bindir
      set_fasta
      fulldb_fasta
      fs 
      os
      

   main:
 
      if (params.dbsplit_update==true) {
         namedb=fulldb_fasta.replaceAll("fa.gz", "")

         get_taxonids(rvdbdir, set_fasta, params.setname)
         get_taxonids_rvdb(rvdbdir, fulldb_fasta, namedb )

         taxonomizator(rvdbdir, get_taxonids.out, bindir , ncbidir)
         taxonomizator_rvdb(rvdbdir, get_taxonids_rvdb.out, bindir, ncbidir )

         filter_FOI(rvdbdir, taxonomizator.out, taxonomizator_rvdb.out, fulldb_fasta, bindir )
         F=filter_FOI.out.FOI
         O=filter_FOI.out.OTHER
      }else{
         F=fs 
         O=os
      }

   emit:
      DS_FOI=F
      DS_OTH=O
} 


workflow linkfiles() {
   take:
      outfoi
      outother
      refseqs_rvdb
      log_file
   
   main:
      families_subset_database=Channel.fromPath(outfoi.toString(), checkIfExists: true)
      nonfamilies_subset_database=Channel.fromPath(outother.toString(), checkIfExists: true)
      stdrd_link(families_subset_database, "$refseqs_rvdb/foi_subset_db.fasta.gz", log_file)
      stdrd_link_otfm(nonfamilies_subset_database, "$refseqs_rvdb/other_subset_db.fasta.gz", log_file) 

      /*
      if (new File(outfoi).exists()) {
         stdrd_link(outfoi, "$refseqs_rvdb/foi_subset_db.fasta.gz", logfil )
      } else {
         error "Families subset fasta file not found in \"$refseqs_rvdb/foi_subset_db.fasta.gz\""
      }

      if (new File(outother).exists()) {  
         stdrd_link(outother, "$refseqs_rvdb/other_subset_db.fasta.gz", logfil ) 
      } else {
            error "Other subset fasta file not found in \"$refseqs_rvdb/other_subset_db.fasta.gz\""
      } */
} 

/*
workflow prepareblast (){
    take:
      db
      logfl
    main:
      //For blast with kaiju:
      db_for_kaiju(db, logfl)
    
    emit:
      ""
} */

// // // // // // MAIN // // // // // //  
    
    
    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Starting  : $workflow.userName  $workflow.start"


workflow {
    def refseqs="${workflow.projectDir}/references"
    def refseqs_rvdb="$refseqs/${params.rvdb_dir}"
    def refseqs_ncbi="$refseqs/db/ncbi"

    def bindir="${workflow.projectDir}/bin"

    def rvdb_fa="${refseqs_rvdb}/${params.db_name}"
    def merged_fa="${refseqs_rvdb}/rvdb+${params.setname}.fasta.gz"
    def rvdb_tax="${refseqs_rvdb}/rvdb+${params.setname}.tax"
    def foisubset="${refseqs_rvdb}/rvdb+${params.setname}_foi_subset.fasta.gz"
    def othersubset="${refseqs_rvdb}/rvdb+${params.setname}_other_subset.fasta.gz"
    
     // create_logf("$workflow.projectDir/references")
        create_logf(refseqs)
        def log_file=create_logf.out
        check_files(rvdb_fa, merged_fa, foisubset)
        if(params.rvdb_update==true) {
          database(refseqs_rvdb, params.rvdb_link, params.db_name, log_file, check_files.out)
          dbfasta=database.out
        }else{
          dbfasta=rvdb_fa
        }
        
        if (params.merge_update==true) {
         mergefasta(refseqs_ncbi, refseqs_rvdb, dbfasta, params.set_seqs, merged_fa, params.taxon_link, create_logf.out)
         ncbidir=mergefasta.out.NCBID
         mergfa=mergefasta.out.MERGED
        }else{
         ncbidir=refseqs_ncbi
         mergfa=merged_fa
        }

      database_subset(refseqs_rvdb, ncbidir, bindir, params.set_seqs, mergfa,foisubset,othersubset)

   
      // linkfiles(database_subset.out.DS_FOI, database_subset.out.DS_OTHER, refseqs_rvdb, log_file)

        //stdrd_link(outfoi, "$refseqs_rvdb/foi_subset_db.fasta.gz", create_logf.out )
        //stdrd_link_otfm(outother, "$refseqs_rvdb/otfm_subset_db.fasta.gz", create_logf.out )
  //    families_subset_database=Channel.fromPath(database_subset.out.DS_FOI.toString(), checkIfExists: true)
  //    nonfamilies_subset_database=Channel.fromPath(database_subset.out.DS_OTH.toString(), checkIfExists: true)
  //    stdrd_link(families_subset_database, "$refseqs_rvdb/foi_subset_db.fasta.gz", log_file)
  //    stdrd_link_otfm(nonfamilies_subset_database, "$refseqs_rvdb/other_subset_db.fasta.gz", log_file  ) 


     // Dbs non-viral discard:
        // db_for_kaiju(params.filtDB, create_logf.out)
    
     // prepareblast(params.blast_kaiju_db, create_logf.out)
}