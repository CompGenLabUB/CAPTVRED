#! /usr/bin/env nextflow

// # --- some params --- #
    params.refseqs       = "${workflow.projectDir}/references"
    params.refseqs_fcsgx = "${params.refseqs}/db/gxdb"


include { db_for_kaiju_pred; create_logf; stdrd_link } from './modules/init/init_conf.nf'
include { stdrd_link as stdrd_link_otfm; stdrd_link as stdrd_link_set} from './modules/init/init_conf.nf'
include { stdrd_link as stdrd_link_tax; stdrd_link as stdrd_link_tax_set; stdrd_link as stdrd_link_info} from './modules/init/init_conf.nf'
include { merge_rvdb_setref; chek_setref_ids; filter_FOI } from './modules/init/filter_rvdb.nf'
include { get_taxonids; get_taxonids_rvdb; set_info_files} from './modules/init/db_taxonomy.nf'
include { taxonomizator; taxonomizator as taxonomizator_rvdb} from './modules/init/db_taxonomy.nf'
include { get_rvdb; get_names_and_nodes;  get_accession2taxid } from './modules/init/update_files.nf' 
// include {download_gxdb; download_fcs_image; FCSGX_FETCHDB} from './modules/fcs_tools.nf'
include {FCSGX_FETCHDB} from './modules/fcs_tools.nf'


if ( params.help ) {
    help = """conf/main.nf: Set up the files for the captvred pipeline run.
             |Required arguments:
             |  --set_seqs  Gziped fasta file with genomic sequences of the species targeted in the TES.
             |  --setname   Name of the set of targeted sequences.
             |              [default: "targetset"]
             |
             |Optional arguments:
             |    --db_update      Force to download rvdb database (in the case it is already downloaded). 
             |                     Most recent version will be always downloaded.
             |                      [default: false]
             |    --merge_update   Force to repeat the merge of reference database and viral candidates sequences if it was already done.
             |                      [default: false] Automatically set to "true" when --taxon_update = true.
             |    --taxon_update   Force to download taxonomic classification files from ncbi (in the case it is already downloaded). 
             |                     Most recent version will be always downloaded.
             |                      [default: false]  Automatically set to "true" when --db_update = true or --taxon_update = true.
             |    
             |Database customization:
             |      
             |      --refdb_link  Database link for downloading desired database for blast (it will be automatically downloaded trough wget). 
             |                      [default: rvdb ]
             |      --kaiju_db    Database(s) to be downloaded for kaiju run. To download multiple dbs use coma separated string. 
             |                      [options: nr_euk, refseqs, viruses, rvdb ]][default: nr_euk ] 
             |                    
    """
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

workflow check_files () {
   take:
      database
      merged
      subset
      
   
   main:

      if( new File(subset).exists()){
          println "Subset DB file found: $subset" 
          params.dbsplit_update = params.dbsplit_update ?: false
      } else {
         println "Subset DB file not found: $subset"
         params.dbsplit_update = true
      }


      if( new File(merged).exists()){
            println "Merged DB file found: $merged" 
            params.merge_update = params.merge_update ?: false
      } else {
         println "Merged DB file not found: $merged"
         params.merge_update = true
      }

      if( new File(database).exists()){
            println "Database file found: $database" 
            params.db_update = params.db_update ?: false 
      } else {
         println "Database file not found: $database"
         params.db_update = true
      }

      if (params.db_update)    params.merge_update=true
      if (params.merge_update) params.dbsplit_update=true
      

      println "## FILES CHECKED! ##"
      println "# Updating database   : ${params.db_update}"
      println "# Merging database    : ${params.merge_update}"
      println "# Subsetting database : ${params.dbsplit_update}"

    // dummy output channel to signal completion
    emit:
        completed_ch = Channel.value(true)

}


workflow database (){
   
   take:
      check
      dir
      link
      name
      logf
      fasta


   main:
   if (params.db_update){
      println "Download database is ${params.db_update}... Downloading RVDB database..."
      get_rvdb(dir, link, name, logf)
      outfl=get_rvdb.out
   } else {
      outfl=fasta
   }
   

   emit:
      outfl

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

      if (params.merge_update){
         println "Merge database is ${params.merge_update}... Merging RVDB and targeted set sequences..."
         get_names_and_nodes(refseqs_ncbi, link, logf)
         chek_setref_ids(refseqs_rvdb, db_fa, set_fa, params.setname)
         merge_rvdb_setref(merged_fa, db_fa, params.db_name, chek_setref_ids.out, logf )
         out_m=merge_rvdb_setref.out
         out_ncbi=get_names_and_nodes.out
      } else {
         out_m=refseqs_ncbi
         out_ncbi=merged_fa
      }

   emit:
      MERGED=out_m
      NCBID=out_ncbi

}


workflow database_subset (){
   take:
      rvdbdir
      ncbidir
      bindir
      gffdir
      set_fasta
      fulldb_fasta
      fs 
      os
      link
      logf
      

   main:
      
      stdrd_link_set(set_fasta, "$rvdbdir/setseqs.fasta.gz", logf)
      if (params.dbsplit_update) {
         println "Download database is ${params.dbsplit_update}... Downloading RVDB database..."
         get_accession2taxid(ncbidir, link, logf) 

         get_taxonids(rvdbdir, set_fasta, params.setname, get_accession2taxid.out)
         get_taxonids_rvdb(rvdbdir, fulldb_fasta, get_accession2taxid.out)

         taxonomizator(rvdbdir, get_taxonids.out.TXID, bindir , ncbidir)
         stdrd_link_tax_set(taxonomizator.out, "$rvdbdir/set_tax.tax.gz", logf)
         taxonomizator_rvdb(rvdbdir, get_taxonids_rvdb.out, bindir, ncbidir )
         stdrd_link_tax(taxonomizator_rvdb.out, "$rvdbdir/full_tax.tax.gz", logf)

         set_info_files(get_taxonids.out.GBFL, taxonomizator.out, bindir, gffdir )
         stdrd_link_info(set_info_files.out, "$rvdbdir/info_summary.tsv", logf)
         filter_FOI(rvdbdir, taxonomizator.out, taxonomizator_rvdb.out, fulldb_fasta, bindir )
         
         F=filter_FOI.out.FOI
         O=filter_FOI.out.OTHER
      }else{
         F=fs 
         O=os
      }
      
      
       
      stdrd_link(F, "$rvdbdir/foi_subset_db.fasta.gz", logf)
      stdrd_link_otfm(O, "$rvdbdir/other_subset_db.fasta.gz", logf) 
      

   emit:
      DS_FOI=F
      DS_OTH=O
} 


// // // // // // MAIN // // // // // //  
    
    
    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Starting  : $workflow.userName  $workflow.start"
    

workflow {
//    def refseqs       = "${workflow.projectDir}/references"
    def refseqs_kai   = "${params.refseqs}/db/kaiju"
    def refseqs_rvdb  = "${params.refseqs}/db/${params.refdb_name}"
    def refseqs_gff   = "$refseqs_rvdb/gff_refgenomes"
    def refseqs_ncbi  = "${params.refseqs}/db/ncbi"

   
    def bindir       = "${workflow.projectDir}/bin"

    def rvdb_fa      = "${refseqs_rvdb}/${params.db_name}"
    def merged_fa    = "${refseqs_rvdb}/rvdb+${params.setname}.fasta.gz"
    def rvdb_tax     = "${refseqs_rvdb}/rvdb+${params.setname}.tax"
    def foisubset    = "${refseqs_rvdb}/rvdb+${params.setname}_foi_subset.fasta.gz"
    def othersubset  = "${refseqs_rvdb}/rvdb+${params.setname}_other_subset.fasta.gz"
    
    def gxdb_rm      = "${params.refseqs_fcsgx}/all.README.txt"
    
         log.info """
    ===========================================
    I N I T   P I P E L I N E
    ===========================================
    FCS-GX DB Path : ${params.refseqs_fcsgx}
    Manifest       : ${params.src_mft}
    ===========================================
    """         
        check_files(rvdb_fa, merged_fa, foisubset)
        view(check_files.out)
        chek_completion=check_files.out.completed_ch
        create_logf(params.refseqs)
        def log_file=create_logf.out

        
        database( chek_completion,
                  refseqs_rvdb, 
                  params.rvdb_link, 
                  params.db_name, 
                  log_file, 
                  rvdb_fa
                    )
         
         mergefasta( refseqs_ncbi, 
                        refseqs_rvdb, 
                        database.out,
                        params.set_seqs, 
                        merged_fa, 
                        params.taxon_link, 
                        create_logf.out
                     )

         database_subset( refseqs_rvdb, 
                          mergefasta.out.NCBID, 
                          bindir, 
                          refseqs_gff,
                          params.set_seqs, 
                          mergefasta.out.MERGED,
                          foisubset,
                          othersubset,
                          params.acc2tax_link,
                          create_logf.out
                     )
        
    FCSGX_FETCHDB(params.src_mft)

   
}


      //  check_gxdb(gxdb_rm)
      //  get_gxdb(check_gxdb.out, 
      //           params.src_mft, 
      //           bindir, 
      //           refseqs_fcsgx )
      
      
      
    // manifest_url = "https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"
    // FCSGX_FETCHDB(manifest_url)

   // if (params.customdb) {
   //          db_for_kaiju(params.dbtomake, log_file)
   //          println "AAAAAA"
   //       } else {
   //          println "EEEEEE"
   //          def kaidatabases=Channel.from(params.kaiju_db)
   //             .splitCsv()
   //             .flatten()
   //          println "$refseqs_kai"
   //         // download_kaijudbs(refseqs_kai, Channel.from(params.kaiju_db).splitCsv().flatten() , log_file )
   //         // download_kaijudbs(refseqs_kai, kaidatabases, log_file )
   //       db_for_kaiju_pred(refseqs_kai, kaidatabases, log_file )
   //       }

// 
// workflow check_gxdb () {
//    take:
//       gxdb_txt
//    
//    main:
//       def do_gxdb
//       if( new File(gxdb_txt).exists()){
//           params.gxdb_update = params.gxdb_update ?: false
//       } else {
//          params.gxdb_update=true
//       }
//       
//     
//     println "##   FCS-GX   FILES CHECKED!     ##"
//     println "# Updating database  :   ${params.gxdb_update} "
//    
//     emit:
//         completed_ch = Channel.value(true)
// 
// }
//     
// 
// 

// workflow get_gxdb (){
//    
//    take:
//       check
//       source_mft
//       bindir
//       dbdir
// 
//    main:
//    if (params.gxdb_update){
//       download_fcs_image()
//       download_gxdb(download_fcs_image.out, source_mft, bindir, dbdir)
//    
//    } 
// 
// }


 // workflow download_kaijudbs () {
 //    take:
 //       odir
 //       datab
 //       log_file
 // 
 //    main:
 //       println "DB is: $datab"
 //       println "Link is: $params.K_nr_euk_link"
 //       def link = ":o"
 //       if (val datab == "nr_euk"){
 //          link = params.K_nr_euk_link
 //       } else if (datab == "refseq") {
 //          link = params.K_refseq_link
 //       } else if (datab == "viruses") {
 //          link = params.K_viruses_link
 //       } else if (datab == "rvdb") {
 //          link = params.K_rvdb_link
 //       } else {
 //          println "WARNING!! Unknown kaiju database! \n Available options are: nr_euk, refseqs, viruses, rvdb"
 //       }
 //      println "Link is: $link" 
 //      db_for_kaiju_pred( odir, link, datab ,log_file )
 // 
 // }