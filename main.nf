#! /usr/bin/env nextflow

include { bbduk_clean } from './rawfq_clean.nf'
include { fastQC; multiQC_raw; multiQC_clean; multiQC_filt; multiQC_bowtie_amp } from './seq_stats.nf'
include { generate_index_bowtie; bowtie_amplicons_alignment; bowtie_amplicons_alignment_sg } from './reads_align.nf'
include { megahit_assembly_all} from './reads_assembly.nf'
include { trinity_assembly_sg; trinity_assembly_pe } from './reads_assembly.nf'
include { make_db_for_blast; do_blastn; do_tblastx; best_reciprocal_hit; blast_sum_coverage } from './contigs_align.nf'
include { kaiju_raw; discard_nonviral; kaiju_contigs; extract_ids; taxonid_to_fasta;  readid_to_fasta} from './taxonomy.nf'
include { coverage_plots; align_counts_plot } from './plots.nf'
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
   fastQC(Channel.fromPath(params.rawfq)) | collect | multiQC_raw
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
 
    bbduk_clean(x, Channel.from(paths_list)) 

    fastQC( bbduk_clean.out.mix() ) | collect | multiQC_clean

    
  emit:
   PE1=bbduk_clean.out.outPE1
   PE2=bbduk_clean.out.outPE2
   SGL=bbduk_clean.out.outSGL

}

workflow reads_filter_nonviral() {
    take:
      x
     
     main:
        kaiju_raw(x)
        discard_nonviral(kaiju_raw.out)
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
      generate_index_bowtie(x,refseqs)
      
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
    
    //multiQC_bowtie_amp(bowtie_amplicons_alignment.out.mix(bowtie_amplicons_alignment_sg.out).collect())
  
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


workflow blast_flow() {

   take:
    ref_fasta
    query_fasta
     
   main:
        make_db_for_blast( ref_fasta) 
        if (params.blast_approach ==~ /(?i)blastn/) {
            do_blastn(query_fasta, make_db_for_blast.out.DB)
            outf=do_blastn.out
        }
        
        if (params.blast_approach ==~ /(?i)tblastx/) {
            do_tblastx(query_fasta, make_db_for_blast.out.DB)
            outf=do_tblastx.out
        }

  emit:
    outf
    
}

workflow blast_flow_rev() {

   take:
     ref_fasta
     query_fasta


   main:

        make_db_for_blast(ref_fasta)
        if ( make_db_for_blast.out.CTRL.toString() == 1){
            if (params.blast_approach ==~ /(?i)blastn/) {
                do_blastn(query_fasta, make_db_for_blast.out.DB.toString())
                outf=do_blastn.out
            }
            
            if (params.blast_approach ==~ /(?i)tblastx/) {
                do_tblastx(query_fasta, make_db_for_blast.out.DB)
                outf=do_tblastx.out
            }
        }else{
            blast_q=query_fasta.toString().split('/')[-1].replaceAll(".gz", "").replaceAll(".fa", "")
            blast_r=ref_fasta.toString().split('/')[-1]
            outf="${params.contigs_blast_dir}/${blast_q}_ON_${blast_r}.${params.blast_approach}.tbl"
            write_out= new File("${outf}")
            write_out.write ""  
            
        }
    
    
  emit:
    outf
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
     ""

}

workflow blast_unc_x_cl (){
    take:
        tax_names
    
    main:
    extract_ids(tax_names)
    extract_ids.out.CLT.view()
    taxonid_to_fasta(extract_ids.out.CLT)
    readid_to_fasta(extract_ids.out.UN)
    
    rfasta=new File(taxonid_to_fasta.out.REF.toString())
    if (rfasta.size() > 0 ) {
        make_db_for_blast(taxonid_to_fasta.out.REF)
        do_blastn(taxonid_to_fasta.out.UNC, make_db_for_blast.out.DB)
        blast_sum_coverage(do_blastn.out, extract_ids.out.UN )
        blast_rpt=blast_sum_coverage.out.TBL
        uncids=blast_sum_coverage.out.UNIDS
    } else {
        blast_rpt=""
        uncids=extract_ids.out.UN
    }
    emit:
      REP=blast_rpt
      UN=uncids
}

workflow blast_unc_x_rvdb () {
    take:
        unaligned_ids
    main:
        
        readid_to_fasta(unaligned_ids)
        ref_database="${params.blast_refseqs_dir}/${params.blast_ref_db_name}"
        do_blastn(readid_to_fasta.out, ref_database)
        blast_sum_coverage(do_blastn.out, unaligned_ids )
        blast_rpt=blast_sum_coverage.out.TBL
        
    emit:
        REP=blast_rpt
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
        align_summary(amplicon_sequences_align.out.ALL)
    
    // 4 // Assembly reads into contigs // //
    
        if (params.assembler ==~ /(?i)MEGAHIT/){
             megahit_assembly_flow(reads_filter_nonviral.out )
             CNFA=megahit_assembly_flow.out.FASTA  //.merge()
             
             brhtbl=megahit_assembly_flow.out.TBL
             Blastout=megahit_assembly_flow.out.BLOUT
             
        }else if (params.assembler =~ /(?i)TRINITY/ ){
             trinity_assembly_flow(reads_filter_nonviral.out)
             CNFA=trinity_assembly_flow.out.FASTA
             brhtbl=trinity_assembly_flow.out.TBL
             Blastout=trinity_assembly_flow.out.BLOUT
        }
    
    // 5 // Taxonomic classification of contigs // //
    
        // 5.1. // Protein-level classification (KAIJU)

            KCDB=Channel.from(params.kaijudbs)
            to_kaiju_contigs=KCDB.combine(CNFA.merge()).merge()
            kaiju_contigs(to_kaiju_contigs)
        
        // 5.2. // Blast of unclass contigs into classified set. 
        
            taxon_outkaiju=kaiju_contigs.out.filter{ it.contains("rvdb")}
            blast_unc_x_cl(taxon_outkaiju)     // Outputs blast report if exists OR empty channel in fothing were found
        
        // 5.3. // Blast of unclass contigs into classified set. 
            blast_unc_x_rvdb(blast_unc_x_cl.out.UN) 
    
    
    // // Blast all contigs into database // //
    
    //ref_database="${params.blast_refseqs_dir}/${params.blast_ref_db_name}"
    //do_blastn(CNFA, ref_database)
    
    
    // 6 // Reporting results
    
        // 6.1 // Plot coverage by genome (reads and contigs) // //
        
        vizualise_results_flow(
                       amplicon_sequences_align.out.PE,
                       amplicon_sequences_align.out.SG,
                       Blastout,
                       brhtbl)
                       
        // 6.2 // Taxonomy summary
            // taxon_outkaiju  (*.rvdb.names.out file)
            // blast_unc_x_cl.RPT
            // blast_unc_x_cl.out.RPT
            
}


workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser...\n" : "Oops .. something went wrong: ${workflow.errorMessage}" )
}
