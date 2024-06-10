#! /usr/bin/env nextflow

include { bbduk_clean; samps_idtranslate } from './modules/rawfq_clean.nf'
include { fastQC; multiQC_raw; multiQC_clean; multiQC_filt; multiQC_bowtie_amp } from './modules/seq_stats.nf'
include { generate_index_bowtie; bowtie_amplicons_alignment; bowtie_amplicons_alignment_sg } from './modules/reads_align.nf'
include { megahit_assembly_all; metaspades_assembly} from './modules/reads_assembly.nf'
include { make_db_for_blast; do_blastn; do_tblastx; blast_sum_coverage; do_cov_onrefseqs; do_cov_on_viralcandidates } from './modules/contigs_align.nf'
include { kaiju_raw; discard_nonviral; kaiju_contigs; kaiju_summarize; extract_ids  } from './modules/taxonomy.nf'
include { coverage_plots; align_counts_plot } from './modules/plots.nf'
include { handle_contamination_pr } from './modules/contamination.nf'
include { fill_html_report; make_summary_tbl } from './modules/sum_and_report.nf'
include { create_logd; create_filesys } from './modules/init.nf'


 //   params.basedir     =  "$workflow.launchDir"
 //   params.bindir      =  "$workflow.projectDir/bin"
 //   params.refseqs     =  "$workflow.projectDir/references"
//
 //   params.db_dir      =  "$params.refseqs/db"
 //   params.dbset_dir   =  "$params.refseqs/$params.setname"
//
 //   params.tmp_dir     =  "${workflow.launchDir}/tmp"
 //   params.rawqc_dir   =  "${params.basedir}/raw"
 //   params.clnfq_dir   =  "${params.basedir}/clean"
 //   params.ampaln_dir  =  "${params.basedir}/aln"
 //   params.asbl_dir    =  "${params.basedir}/assembly"
 //   params.taxdir      =  "${params.basedir}/taxonomy"
 //   params.reports_dir =  "${params.basedir}/reports"
 //   params.logs_dir    =  "${params.basedir}/logs"
 //   params.html_dir    =  "${params.basedir}/html"


def samplesMap = [:]
SamplesDef = file(params.samp) // (samplestbl_file)
SamplesDef.eachLine {
    line -> {
        def sampl = line.split('\t')
        // ignore lines stating with "#"
        if (!sampl[0].startsWith("#")) {
            samplesMap.(sampl[0]) = (sampl[1])
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
 Parameters 
 ---------------------------------------
    Minimum Contigs Length : ${params.assemblyMINCONLEN}
 =======================================
 """

if ( params.help ) {
    help = """your_script.nf: A description of your script and maybe some examples of how
             |                to run the script
             |Required arguments:
             |  --input_file  Location of the input file file.
             |                [default: ${params.input_file}]
             |
             |Optional arguments:
             |  --use_thing   Do some optional process.
             |                [default: ${params.use_thing}]
             |  -w            The NextFlow work directory. Delete the directory once the process
             |                is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}


workflow init_run() {

  take:
    directories
  main:
    println "# INIT: $samplesMap"

    create_logd(params.logs_dir)
      def fsys = Channel.of( directories)
 //   def fsys = Channel.of( params.tmp_dir,
 //                          params.rawfq_dir,
 //                          params.clnfq_dir,
 //                          params.ampaln_dir,
 //                          params.asbl_dir,
 //                          params.subasb_dir,
 //                          params.taxdir,
 //                          params.subtax_dir,
 //                          params.reports_dir)
    fsys.view() 
    create_filesys(fsys, create_logd.out)

   /* if (params.assembler ==~ /(?i)MEGAHIT/){
        params.subasb_dir="$params.asbl_dir/megahit"
    }else if (params.assembler ==~ /(?i)METASPADES/ ) {
        params.subasb_dir="$params.asbl_dir/metaspades"
    }

    if (params.taxalg ==~ /(?i)KAIJU/ ) {
      params.subtax_dir="$params.taxdir/kaiju"
    } else if (params.taxalg ==~ /(?)BLASTN/ ){
      params.subtax_dir="$params.taxdir/blastn"
    } else if (params.taxalg ==~ /(?)TBLASTX/ ){ 
      params.subtax_dir="$params.taxdir/tblastx"
    }
*/

  emit:
    create_filesys.out
}


workflow fastqc_onrawseqs() {

  take:
    logfl

  main:
   
    def thysamples = samplesMap
    def ids=[]
    thysamples.each { sampleID, illuminaID ->
         ids << illuminaID
    }
    def spstr=ids.join(",")
    def regx="$params.fastq_dir/{$spstr}$params.rawfq_sfx"
    
    println "### $spstr ###"
    spschan=Channel.fromPath("$regx")
    fastQC(spschan) | collect | multiQC_raw
}


workflow reads_clean() {

  take:

    x

  main:
    
    
    def paths_list=[]
         //paths_list is a LoL: where first item is the raw fastq root 
         //for each sample and the second one is the new root for the cleanseqs.
    def thysamples = samplesMap
    thysamples.each { sampleID, illuminaID ->
         def newsamp=["$params.fastq_dir/$illuminaID", "$params.clnfq_dir/$sampleID"]
         paths_list << newsamp
    }samps_idtranslate
    
    
    if (params.trim_adapters == true ) {
       bbduk_clean(x, Channel.from(paths_list)) 
       fastQC( bbduk_clean.out.mix() ) | collect | multiQC_clean
       out1=bbduk_clean.out.outPE1
       out2=bbduk_clean.out.outPE2
       outS=bbduk_clean.out.outSGL
    } else {
       samps_idtranslate(x, Channel.from(paths_list))
       out1=samps_idtranslate.out.outPE1
       out2=samps_idtranslate.out.outPE2
       outS=samps_idtranslate.out.outSGL
    }
    
  emit:
   PE1=out1
   PE2=out2
   SGL=outS

}

workflow idtranslate() {

  take:

    x

  main:
    
    
    def paths_list=[]
         //paths_list is a LoL: where first item is the raw fastq root 
         //for each sample and the second one is the new root for the cleanseqs.
    def thysamples = samplesMap
    thysamples.each { sampleID, illuminaID ->
         def newsamp=["$params.rawfq_dir/$illuminaID", "$params.clnfq_dir/$sampleID"]
         paths_list << newsamp
    }
    
    samps_idtranslate(x, Channel.from(paths_list))
    
  emit:
   PE1=samps_idtranslate.out.outPE1
   PE2=samps_idtranslate.out.outPE2
   SGL=samps_idtranslate.out.outSGL

}

workflow reads_filter_nonviral() {
    take:
      x
     
     main:
        kaiju_raw(x)
        // kaiju_raw.out.view()
        discard_nonviral(kaiju_raw.out)
        // discard_nonviral.out.SGLout.view()
        fastQC( discard_nonviral.out.mix() ) | collect | multiQC_filt

    emit:
       discard_nonviral.out.PE1out
       discard_nonviral.out.PE2out
       discard_nonviral.out.SGLout
}

/*
workflow amplicon_sequences_dbinit() {

    take:
      x
      refseqs

    main:
      generate_index_bowtie(refseqs)
      
    emit:
      generate_index_bowtie.out

}
*/

workflow reads_align_wf() {

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
    ALL=bowtie_amplicons_alignment.out.mix(bowtie_amplicons_alignment_sg.out).collect()
    PE=bowtie_amplicons_alignment.out
    SG=bowtie_amplicons_alignment_sg.out
}

workflow align_summary() {
    take:
        BAMLIST
    
    main:

     multiQC_bowtie_amp(BAMLIST)
     align_counts_plot(BAMLIST)
        
}


workflow trinity_assembly_flow () {
    take: 
      pe1
      pe2
      sgl      
    
    main:
      
      ref_fa="${params.amplicon_refseqs}"
      trinity_assembly_pe(pe1, pe2, sgl)
      blast_flow(ref_fa, trinity_assembly_pe.out)
      blast_flow_rev(trinity_assembly_pe.out, ref_fa)
      
      QOR=blast_flow.out
      ROQ=blast_flow_rev.out
      
      best_reciprocal_hit(trinity_assembly_pe.out.CGSout, ref_fa, QOR, ROQ )

    emit:
      TBL=best_reciprocal_hit.out
      FASTA=trinity_assembly_pe.out.CGSout
      BLOUT=QOR
}


workflow vizualise_results_flow() {
    take:
        pebam
        sgbam
        blasttbl
        brh
        
    main:
        coverage_plots(pebam, sgbam, blasttbl, brh)
        
    emit:
        coverage_plots.out.ODIR

}



workflow direct_blast_n () {
    take:
        ref_fasta
        all_contigs
    main:

        make_db_for_blast( ref_fasta, "FALSE") 
        do_blastn(all_contigs, make_db_for_blast.out.DB, params.taxbndir)
        
        if (params.handle_contamination == true ) {
            handle_contamination_pr( params.cids, 
                                     params.cfaa, 
                                     do_blastn.out.OUT,
                                     do_blastn.out.CONTIGS )
                                     
            blastOut=handle_contamination_pr.out
        } else {
        
            blastOut=do_blastn.out.OUT
        }
        
        blast_sum_coverage(blastOut, "F", "F" )
        
    emit:
        BY_R=blast_sum_coverage.out.BYR
        BY_SQ=blast_sum_coverage.out.BYSQ
        BY_SP=blast_sum_coverage.out.BYSP
        S_SUM=blast_sum_coverage.out.SUM
        CFA=all_contigs
        DONE=blast_sum_coverage.out.SUM2
}

workflow direct_blast_tx () {
    take:
        ref_fasta
        all_contigs
    main:

        make_db_for_blast( ref_fasta, "FALSE") 
        do_tblastx(all_contigs, make_db_for_blast.out.DB, params.taxtbxdir)
        
        if (params.handle_contamination == true ) {
            handle_contamination_pr( params.cids, 
                                     params.cfaa, 
                                     do_tblastx.out.OUT,
                                     do_tblastx.out.CONTIGS )
                                     
            blastOut=handle_contamination_pr.out
        } else {
        
            blastOut=do_tblastx.out.OUT
        }
        
        blast_sum_coverage(blastOut, "F", "F" )
        
    emit:
        BY_R=blast_sum_coverage.out.BYR
        BY_SQ=blast_sum_coverage.out.BYSQ
        BY_SP=blast_sum_coverage.out.BYSP
        S_SUM=blast_sum_coverage.out.SUM
        CFA=all_contigs
        DONE=blast_sum_coverage.out.SUM2
}


workflow coverage_onrefseqs() {
    
    take:
        contigs_fa
        assign_byr

    main:
     
     reffa="${params.amplicon_refseqs_dir}/${params.amplicon_refseqs}"
     make_db_for_blast( reffa, "FALSE")   // FALSE refers to: not redo db if it already exists.
     do_cov_on_viralcandidates("${params.amplicon_refseqs_dir}/${params.amplicon_refseqs_info}", 
                     make_db_for_blast.out.DB,
                     contigs_fa,
                     assign_byr,
                     params.bl_suffix
                     )

   
    emit:
      blastout=do_cov_on_viralcandidates.out.BLOUT
      coverage=do_cov_on_viralcandidates.out.COV

}

workflow set_dep_params () {

    params.basedir     =  workflow.launchDir
    params.bindir      =  workflow.projectDir+'/bin'
    params.refseqs     =  workflow.projectDir+'/references'

    params.db_dir      =  params.refseqs + '/db'
    params.dbset_dir   =  params.refseqs + '/$params.setname'

    if (!params.tmp_dir) { params.tmp_dir     =  "${params.basedir}/tmp" }
    params.rawqc_dir   =  "${params.basedir}/raw"
    params.clnfq_dir   =  "${params.basedir}/clean"
    params.ampaln_dir  =  "${params.basedir}/aln"
    params.asbl_dir    =  "${params.basedir}/assembly"
    params.taxdir      =  "${params.basedir}/taxonomy"
    params.reports_dir =  "${params.basedir}/reports'"    
    params.logs_dir    =  "${params.bdir}/logs"
    params.html_dir    =  "${params.basedir}/html"
  
  emit:
    [params.tmp_dir, params.rawfq_dir, params.clnfq_dir, params.ampaln_dir]
 //   params.asbl_dir
 //   params.subasb_dir
 //   params.taxdir
 //   params.subtax_dir
 //   params.reports_dir

}

// // // // // // MAIN // // // // // //  

workflow {

    
    //check_params()
    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Bdir      : $workflow.launchDir"
    println "# Starting  : $workflow.userName $ZERO $workflow.start"
    println "# Reading samples for $params.runID from $params.samp"
    
    set_dep_params()

    println " ### $workflow.launchDir ## $params.basedir ## $params.tmp_dir ## $params.html_dir"
    //set_dep_params() // | collect | init_run()
    //Channel.of(set_dep_params.out).view()
    // def d = set_dep_params.out.collect()
 /*   init_run("ll") 
    
    // 1 // Clean reads // //
    
        fastqc_onrawseqs(init_run.out)
        reads_clean(init_run.out)
    
    // 2 // Discard reads identified as nonviral // //
    
        //KDB=Channel.from(params.kaijudbs)
        KDB=Channel.from(params.filtDB)
        CLNR=reads_clean.out.merge()
        to_kaiju=KDB.combine(CLNR).merge()
        reads_filter_nonviral(to_kaiju)

    // 3 // Map reads on database   //  //
        dbtobowtie=init_refs.out.FULLDB
        generate_index_bowtie(dbtobowtie)
        reads_align_wf(generate_index_bowtie.out, reads_filter_nonviral.out)
*//*    // 3 // Map reads on reference genomes fasta // //  (used to design the probes) 
    
        tobowtie="${params.amplicon_refseqs_dir}/${params.amplicon_refseqs}";
        generate_index_bowtie(tobowtie)
        amplicon_sequences_align(generate_index_bowtie.out, reads_filter_nonviral.out)
    //    align_summary(amplicon_sequences_align.out.ALL)
    
    // 4 // Assembly reads into contigs // //
    
        if (params.assembler ==~ /(?i)MEGAHIT/){
             megahit_assembly_all(reads_filter_nonviral.out)
             CNFA=megahit_assembly_all.out.CGSout

             
        }else if (params.assembler ==~ /(?i)METASPADES/ ){
             metaspades_assembly(reads_filter_nonviral.out)
             CNFA=metaspades_assembly.out.CGSout

        }
    
    // 5 // Taxonomic classification of contigs // //
        
        if (params.taxalg ==~ /(?i)KAIJU/ ) {
            // 5.1. FAST APPROACH// 
              // 5.1.1. Protein-level classification (KAIJU) //
            KCDB=Channel.from(params.kaijudb)
            // KCDB.view()
            to_kaiju_contigs=KCDB.combine(CNFA.merge()).merge()
            kaiju_contigs(to_kaiju_contigs)
            kaiju_summarize(kaiju_contigs.out.NM)
        //    kaiout=kaiju_contigs.out.filter{it.contains(params.blastfastdb)}.flatten()
        //    kaiout.view()
        //    all_contigs_blast(kaiout)
                     
                       
            blOUT=""                   
            covOUT=""
            byreadfl=kaiju_summarize.out.BYR
            byseqfl=kaiju_summarize.out.BYSQ
            byspecfl=""
            statssum=kaiju_summarize.out.SUM
            kronaPT=kaiju_contigs.out.KP

            FIN=kaiju_summarize.out.DONE.collect()
          
        } else if (params.taxalg ==~ /(?)BLASTN/ ){
            //Directly blast into database //
            ref_fa="${params.blast_refseqs_dir}/${params.blast_ref_db_name}";
            // blast_flow( ref_fa, CNFA.merge())

            // direct_blast(ref_fa, CNFA.merge())
            direct_blast_n(ref_fa, CNFA)
            
             if ( params.general_only == false) {
                  coverage_onrefseqs(direct_blast_n.out.CFA, 
                                     direct_blast_n.out.BY_R
                                    )
                  blOUT=coverage_onrefseqs.out.blastout                   
                  covOUT=coverage_onrefseqs.out.coverage
             } else {
                  blOUT=""                   
                  covOUT=""
             }
            
            byreadfl=direct_blast_n.out.BY_R
            byseqfl=direct_blast_n.out.BY_SQ
            byspecfl=direct_blast_n.out.BY_SP
            statssum=direct_blast_n.out.S_SUM
            kronaPT=""     

            FIN=direct_blast_n.out.DONE.collect()           
           
         } else if (params.taxalg ==~ /(?)TBLASTX/ ){

            // Directly blast into database //
            ref_fa="${params.blast_refseqs_dir}/${params.blast_ref_db_name}";
            // blast_flow( ref_fa, CNFA.merge())

            // direct_blast(ref_fa, CNFA.merge())
            direct_blast_tx(ref_fa, CNFA)
            
            if ( params.general_only == false) {
                  coverage_onrefseqs(direct_blast_tx.out.CFA, 
                                     direct_blast_tx.out.BY_R
                     )
                  blOUT=coverage_onrefseqs.out.blastout                   
                  covOUT=coverage_onrefseqs.out.coverage
            } else {
                  Channel.from("TAXONOMY APPROACH IS NOT AVAILABLE... please check your spelling!").view()
                  blOUT=""                   
                  covOUT=""
            }
            
            byreadfl=direct_blast_tx.out.BY_R
            byseqfl=direct_blast_tx.out.BY_SQ
            byspecfl=direct_blast_tx.out.BY_SP
            statssum=direct_blast_tx.out.S_SUM
            kronaPT=""

            FIN=direct_blast_tx.out.DONE.collectFile(name: "${params.tmp_dir}/allstats.txt", newLine: true)
            // collect()
         }


    // 6 // Reporting results
    
        // 6.1. // Summary tables // //
        
        make_summary_tbl(FIN)
        
        // 6.2. // Plot coverage by genome (reads and contigs)
        
        vizualise_results_flow(
                       reads_align_wf.out.PE,
                       reads_align_wf.out.SG,
                       blOUT,
                       covOUT
                       )
        figsDIR=vizualise_results_flow.out
        
        // 6.3. // HTML report. 
        fill_html_report(make_summary_tbl.out, figsDIR)
        
                       
        // 6.2 // Taxonomy summary
            // taxon_outkaiju  (*.rvdb.names.out file)
            // blast_unc_x_cl.RPT
            // blast_unc_x_cl.out.RPT
*/         
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser...\n" : "Oops .. something went wrong: ${workflow.errorMessage}" )
}