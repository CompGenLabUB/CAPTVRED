#! /usr/bin/env nextflow

include { bbduk_clean; samps_idtranslate } from './rawfq_clean.nf'
include { fastQC; multiQC_raw; multiQC_clean; multiQC_filt; multiQC_bowtie_amp } from './seq_stats.nf'
include { generate_index_bowtie; bowtie_amplicons_alignment; bowtie_amplicons_alignment_sg } from './reads_align.nf'
include { megahit_assembly_all; metaspades_assembly} from './reads_assembly.nf'
include { index_seqs; index_seqs as index_refs; make_db_for_blast; do_blastn; do_tblastx; blast_sum_coverage; do_cov_onrefseqs; do_cov_on_viralcandidates } from './contigs_align.nf'
include { kaiju_raw; discard_nonviral; kaiju_contigs; kaiju_summarize; extract_ids  } from './taxonomy.nf'
include { coverage_plots; align_counts_plot } from './plots.nf'
include { handle_contamination_pr } from './contamination.nf'
include { fill_html_report; make_summary_tbl } from './sum_and_report.nf'

/* THE DISCARDED:
include { do_blast_kaiju; merge_blast_outs } from './contigs_align.nf'
include { taxonid_to_fasta; readid_to_fasta }from './taxonomy.nf'
include { best_reciprocal_hit } from './contgis_align'
*/

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
 Parameters 
 ---------------------------------------
    Minimum Contigs Length : ${params.assemblyMINCONLEN}
 =======================================
 """

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
   
    def thysamples = samplesMap
    def ids=[]
    thysamples.each { sampleID, illuminaID ->
         ids << illuminaID
    }
    def spstr=ids.join(",")
    def regx="$params.rawfq_dir/{$spstr}$params.rawfq_sfx"
    
    println "### $spstr ###"
    //def regexp= "^($spstr)" << "(${R1}|${R2})*.gz" 
    // Channel.fromPath(params.rawfq).view()
    spschan=Channel.fromPath("$regx")
    //spschan=spsch.filter({ it =~ /"$spstr"/ })
    //spschan=spsch.grep({ it.toString() =~ /"$spstr"/ })
    // spschan.view()
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
         def newsamp=["$params.rawfq_dir/$illuminaID", "$params.clnfq_dir/$sampleID"]
         paths_list << newsamp
    }
    
    
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
  //  bbduk_clean(x, Channel.from(paths_list)) 

  //  fastQC( bbduk_clean.out.mix() ) | collect | multiQC_clean

    
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


workflow amplicon_sequences_dbinit() {

    take:
      x
      refseqs

    main:
      //generate_index_bowtie(x,refseqs)
      generate_index_bowtie(refseqs)
      
    emit:
      generate_index_bowtie.out


}


workflow amplicon_sequences_align() {

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



/*
workflow megahit_assembly_flow () {
    take:
      pe1
      pe2
      sgl
    
    main:
      ref_fa="${params.amplicon_refseqs_dir}/${params.amplicon_refseqs}";
      megahit_assembly_all(pe1,pe2, sgl)
      
      blast_flow( ref_fa, megahit_assembly_all.out.CGSout)
      blast_flow_rev( megahit_assembly_all.out.CGSout, ref_fa)

      QOR=blast_flow.out
      ROQ=blast_flow_rev.out
          
      best_reciprocal_hit(megahit_assembly_all.out.CGSout, ref_fa, QOR, ROQ )
      brh=best_reciprocal_hit.out

    emit:
      TBL=brh
      FASTA=megahit_assembly_all.out.CGSout
      BLOUT=QOR
}
*/

workflow trinity_assembly_flow () {
    take: 
      pe1
      pe2
      sgl      
    
    main:
      
      ref_fa="${params.amplicon_refseqs}"
      trinity_assembly_pe(pe1, pe2, sgl)
      //trinity_assembly_sg(sgl)
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


//workflow blast_flow() {
//   take:
//    ref_fasta
//    query_fasta
     
//   main:
//        make_db_for_blast( ref_fasta, "FALSE") 
//        if (params.blast_approach ==~ /(?i)blastn/) {
//            do_blastn(query_fasta, make_db_for_blast.out.DB, params.taxslowdir)
//            outf=do_blastn.out
//        }
        
//        if (params.blast_approach ==~ /(?i)tblastx/) {
//            do_tblastx(query_fasta, make_db_for_blast.out.DB)
//            outf=do_tblastx.out
//        }
//  emit:
//    outf
//}

//workflow blast_flow_rev() {
//   take:
//     ref_fasta
//     query_fasta
//   main:
//        make_db_for_blast(ref_fasta, "FALSE")
//        if ( make_db_for_blast.out.CTRL.toString() == 1){
//            if (params.blast_approach ==~ /(?i)blastn/) {
//                do_blastn(query_fasta, make_db_for_blast.out.DB.toString(), params.taxfastdir)
//                outf=do_blastn.out
//            }
//            
//            if (params.blast_approach ==~ /(?i)tblastx/) {
//                do_tblastx(query_fasta, make_db_for_blast.out.DB)
//                outf=do_tblastx.out
//            }
//        }else{
//            blast_q=query_fasta.toString().split('/')[-1].replaceAll(".gz", "").replaceAll(".fa", "")
//            blast_r=ref_fasta.toString().split('/')[-1]
//            outf="${params.contigs_blast_dir}/${blast_q}_ON_${blast_r}.${params.blast_approach}.tbl"
//            write_out= new File("${outf}")
//            write_out.write ""  
//        }
//    emit:
//    outf
//}


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
        
        //readid_to_fasta(unaligned_ids)
        //ref_database="${params.blast_refseqs_dir}/${params.blast_ref_db_name}"

        make_db_for_blast( ref_fasta, "FALSE") 
        // all_contigs.view()
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
        // blast_byread=blast_sum_coverage.out.BYR
        
    emit:
        // REP=blast_byread
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
        
        //readid_to_fasta(unaligned_ids)
        //ref_database="${params.blast_refseqs_dir}/${params.blast_ref_db_name}"

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
        // blast_byread=blast_sum_coverage.out.BYR
        
    emit:
        // REP=blast_byread
        BY_R=blast_sum_coverage.out.BYR
        BY_SQ=blast_sum_coverage.out.BYSQ
        BY_SP=blast_sum_coverage.out.BYSP
        S_SUM=blast_sum_coverage.out.SUM
        CFA=all_contigs
        DONE=blast_sum_coverage.out.SUM2
}

/*
workflow unc_contigs_blast (){
    take:
        uncids
    
    main:
        sp=uncids.toString().split('/')[-1].split('[.]')[0]

        // Queries fasta:
        readid_to_fasta(uncids)
        
        // Reference fasta:
        ref_fa="${params.blast_refseqs_dir}/${params.blast_ref_db_name}";
             
        // Blast and output processing:
            
        make_db_for_blast(ref_fa, "FALSE")
        //  make_db_for_blast.out.DB.view()
        do_tblastx(readid_to_fasta.out.FA, make_db_for_blast.out.DB, params.taxfastdir)
       // blast_sum_coverage(do_blastn.out, uncids )
       // blast_rpt=blast_sum_coverage.out.TBL
       // stillunc=blast_sum_coverage.out.UNIDS
       
    emit:
          SP=sp
          BLOUT=do_tblastx.out
         // tuple sp, val(do_blastn.out)
         
}
*/

/* workflow class_contigs_blast (){
    take:
        clids
        taxids
        
    main:
       // sp=clids.toString().split('/')[-1].split('[.]')[0]
        
        
            // Queries fasta:
            readid_to_fasta(clids)
            
            // refseqs fasta:
            taxonid_to_fasta(taxids)
            ref_fa=taxonid_to_fasta.out.REF
            
            // Blast and output processing:
            
            make_db_for_blast(ref_fa, "TRUE")
            do_blast_kaiju(readid_to_fasta.out.FA, params.taxfastdir, make_db_for_blast.out.DB)

    emit:
        BLOUT=do_blast_kaiju.out

} */

/* workflow all_contigs_blast (){
    take:
      kaijuout
        
    main:
        
            extract_ids(kaijuout)
            
           // a. // Blast of unclass contigs into reference db. 
                unc_contigs_blast(extract_ids.out.UNIDS)
                
           // b. // Blast of class contigs into reference sequences.
                class_contigs_blast(extract_ids.out.CLIDS, extract_ids.out.CLT)
                
           // c. // Summarize results
                merge_blast_outs(class_contigs_blast.out.BLOUT, unc_contigs_blast.out.BLOUT )
                // merge_blast_outs.out.BALL.view()
                blast_sum_coverage(merge_blast_outs.out.BALL, extract_ids.out.CLIDS, extract_ids.out.UNIDS)

    
    emit:
        ""

} */

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

// // // // // // MAIN // // // // // //  

workflow {

    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Starting  : $workflow.userName $ZERO $workflow.start"
    println "# Reading samples for $params.runID from $params.sampletbl"
    
    // 1 // Clean reads // //
   
        init_samples()
        fastqc_onrawseqs(init_samples.out)
        reads_clean(init_samples.out)
   
    // 2 // Discard reads identified as nonviral // //
   
        //KDB=Channel.from(params.kaijudbs)
        KDB=Channel.from(params.kaijuDBRAW)
        CLNR=reads_clean.out.merge()
        to_kaiju=KDB.combine(CLNR).merge()
        reads_filter_nonviral(to_kaiju)
        
    // 3 // Map reads on reference genomes fasta // //  (used to design the probes) 
    
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
                       amplicon_sequences_align.out.PE,
                       amplicon_sequences_align.out.SG,
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
            


}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser...\n" : "Oops .. something went wrong: ${workflow.errorMessage}" )
}

///// DISCARDED //////
 // 5.2. // Blast of unclass contigs into classified set. 
        
       // taxon_outkaiju=kaiju_contigs.out.filter{ it.contains("rvdb")}
       // blast_unc_x_cl(taxon_outkaiju)     // Outputs blast report if exists OR empty channel in nothing were found
