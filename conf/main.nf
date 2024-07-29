#! /usr/bin/env nextflow

include { db_for_kaiju; create_logf; stdrd_link } from './modules/init_conf.nf'
include { stdrd_link as stdrd_link_otfm; stdrd_link as stdrd_link_set} from './modules/init_conf.nf'
include { stdrd_link as stdrd_link_tax; stdrd_link as stdrd_link_tax_set; stdrd_link as stdrd_link_info} from './modules/init_conf.nf'
include { merge_rvdb_setref; chek_setref_ids; filter_FOI } from './modules/filter_rvdb.nf'
include { get_taxonids; get_taxonids_rvdb; set_info_files} from './modules/db_taxonomy.nf'
include { taxonomizator; taxonomizator as taxonomizator_rvdb} from './modules/db_taxonomy.nf'
include { get_rvdb; get_names_and_nodes;  get_accession2taxid } from './modules/update_files.nf' 

log.info """\
 =========================================
    C A P T V R E D   -   S E T    U P
 =========================================
 SET NAME    : ${params.runID}
 FASTA FILE  : ${samplesMap}
 SysInfo     : ${workflow.userName} SID=${workflow.sessionId} NCPUs=${params.NCPUS} GITcid=${workflow.commitId}
 =========================================
 """

if ( params.help ) {
    help = """conf/main.nf: Set up the files for the captvred pipeline run.
             |Required arguments:
             |  --set_seqs  Gziped fasta file with genomic sequences of the species targeted in the TES.
             |  --setname   Name of the set of targeted sequences.
             |              [default: "targetset"]
             |
             |Optional arguments:
             |    --rvdb_update    Force to download rvdb database (in the case it is already downloaded). Most recent version will be always downloaded.
             |                      [default: false]
             |    --merge_update   Force to repeat the merging process of reference database and viral candidates sequences (in the case it is already done).
             |                      [default: false] Automatically set to "true" when --taxon_update = true.
             |    --taxon_update   Force to download taxonomic classification files from ncbi (in the case it is already downloaded). Most recent version will be always downloaded.
             |                      [default: false]  Automatically set to "true" when --rvdb_update = true or --taxon_update = true.
             |    
             |    Database customization:
             |      [default: rvdb ]
             |      --refdb_link  Database link for downloading desired database (it will be automatically downloaded trough wget)
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
      gffdir
      set_fasta
      fulldb_fasta
      fs 
      os
      link
      logf
      

   main:
      
      stdrd_link_set(set_fasta, "$rvdbdir/setseqs.fasta.gz", logf)
      if (params.dbsplit_update==true) {
         //namedb=fulldb_fasta.replaceAll("fasta.gz|fa.gz", "")
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
    def refseqs      = "${workflow.projectDir}/references"
    def refseqs_rvdb = "$refseqs/db/${params.refdb_name}"
    def refseqs_gff  = "$refseqs_rvdb/gff_refgenomes"
    def refseqs_ncbi = "$refseqs/db/ncbi"
   
    def bindir       = "${workflow.projectDir}/bin"

    def rvdb_fa      = "${refseqs_rvdb}/${params.db_name}"
    def merged_fa    = "${refseqs_rvdb}/rvdb+${params.setname}.fasta.gz"
    def rvdb_tax     = "${refseqs_rvdb}/rvdb+${params.setname}.tax"
    def foisubset    = "${refseqs_rvdb}/rvdb+${params.setname}_foi_subset.fasta.gz"
    def othersubset  = "${refseqs_rvdb}/rvdb+${params.setname}_other_subset.fasta.gz"
    
     // create_logf("$workflow.projectDir/references")
        create_logf(refseqs)
        def log_file=create_logf.out
        check_files(rvdb_fa, merged_fa, foisubset)
        

        db_for_kaiju(params.raw_kaiju_db, log_file)
        
        if(params.rvdb_update==true) {

            database( refseqs_rvdb, 
                      params.rvdb_link, 
                      params.db_name, 
                      log_file, 
                      check_files.out
                    )
            dbfasta=database.out

        }else{

            dbfasta=rvdb_fa

        }
        
        if (params.merge_update==true) {
            mergefasta( refseqs_ncbi, 
                        refseqs_rvdb, 
                        dbfasta, 
                        params.set_seqs, 
                        merged_fa, 
                        params.taxon_link, 
                        create_logf.out
                     )

            ncbidir=mergefasta.out.NCBID
            mergfa=mergefasta.out.MERGED

        }else{

            ncbidir=refseqs_ncbi
            mergfa=merged_fa
        
        }

      database_subset( refseqs_rvdb, 
                       ncbidir, 
                       bindir, 
                       refseqs_gff,
                       params.set_seqs, 
                       mergfa,
                       foisubset,
                       othersubset,
                       params.acc2tax_link,
                       create_logf.out
                     )


}