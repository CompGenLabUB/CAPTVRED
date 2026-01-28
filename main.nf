#! /usr/bin/env nextflow
nextflow.enable.dsl=2
params.run_fcsgx = false 

include { bbduk_clean; samps_idtranslate } from './modules/rawfq_clean.nf'
//include { fastQC; multiQC_raw; multiQC_clean; multiQC_filt; multiQC_bowtie_amp } from './modules/seq_stats.nf'
// include { fastQC; multiQC_bowtie } from './modules/seq_stats.nf'
include { fastQC as fastQC_raw; fastQC as fastQC_clean; fastQC as fastQC_filt } from './modules/seq_stats.nf'
include { multiQC as multiQC_raw; multiQC as multiQC_clean; multiQC as multiQC_filt } from './modules/seq_stats.nf'
include { fq2fasta; discard_contaminants } from './modules/seq_stats.nf'
include { FCSGX_RUNGX } from './modules/nf-core/fcsgx/main.nf'
include { generate_index_bowtie; bowtie_amplicons_alignment; bowtie_amplicons_alignment_sg } from './modules/reads_align.nf'
include { megahit_assembly_new; refine_assembly; metaspades_assembly_new} from './modules/reads_assembly.nf'
//include { make_db_for_blast; do_blastn; do_tblastx; blast_sum_coverage; do_cov_on_viralcandidates } from './modules/contigs_align.nf'
include { blast_makedb; blast_search; blast_sum_coverage } from './modules/contigs_align.nf'
include { kaiju_raw; discard_nonviral; kaiju_contigs; kaiju_summarize; extract_ids  } from './modules/taxonomy.nf'
include { coverage_plots; align_counts_plot } from './modules/plots.nf'
include { handle_contamination_pr } from './modules/contamination.nf'
include { fill_html_report; make_summary_tbl } from './modules/sum_and_report.nf'
include { create_logd; create_filesys } from './modules/init.nf'



log.info """\
 ========================================
  C A P T V R E D - N F   P I P E L I N E
 ========================================
 RUN     : ${params.runID}
 Samples : 
 SysInfo : ${workflow.userName} SID=${workflow.sessionId} NCPUs=${params.cpus} GITcid=${workflow.commitId}
 ========================================

 """

if ( params.help ) {
    help = """main.nf: Start CAPTVRED pipeline for TES datasets sequences analyses.
             |Required arguments:
             |    --samp       Tabulated file with samples information.Templete named "samples_Definition_template.tbl" is provided with this repository.
             |    --fastq_dir  Root directory for read fastq files.
             |    --runID      Identification (name or code) of the run.
             | 
             |
             |Optional arguments:
             |  
             |    --help          Print this help message
             |    --NCPUS         Maximum number of threads available for the demanding steps of the pipeline. [default:32]
             |    --R1            Suffix for read 1 fastq files. [default: _R1_001]
             |    --R2            Suffix for read 2 fastq files. [default: _R2_001]
             |    --assembler     megahit(default) or metaspades
             |    --taxalg        Taxonomy algorithm. blastn(default), tblastx or kaiju
             |    --handle_contamination [default: false]. If true, fasta file must be provided.
             |    --do_cov_figures       [default: true ]
    """
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}


workflow reads_align_wf() {

   take:
     cluster_index_path
     pe1
     pe2
     sgl
    
   main:
   
    bowtie_amplicons_alignment(cluster_index_path, pe1, pe2, sgl, params.ampaln_dir)
    bowtie_amplicons_alignment_sg(cluster_index_path, pe1, pe2, sgl, params.ampaln_dir )
    // samps_align=multiQC_bowtie_amp(bowtie_amplicons_alignment.out.mix(bowtie_amplicons_alignment_sg.out).collect())
    peout=bowtie_amplicons_alignment.out.LOG.mix(bowtie_amplicons_alignment.out.STS)
    sgout=bowtie_amplicons_alignment_sg.out.LOG.mix(bowtie_amplicons_alignment_sg.out.STS)
    samps_align=peout.mix(sgout).collect()
   // samps_align.view()
    multiQC_bowtie( samps_align, params.reports_dir, params.logs_dir, "align")
  
  emit:
    ALL=samps_align
    PE=bowtie_amplicons_alignment.out.bamPE
    SG=bowtie_amplicons_alignment_sg.out.bamSG
}




workflow vizualise_results_flow() {
    take:
        pebam
        sgbam
        blasttbl
        brh
        bindir
        repdir
        dbdir
        gffdir
        
    main:
        coverage_plots(pebam, 
                       sgbam, 
                       blasttbl, 
                       brh, 
                       bindir, 
                       repdir, 
                       dbdir, 
                       gffdir
                       )
        
    emit:
        coverage_plots.out.ODIR

}

workflow direct_blastn_new () {
    take:
        ref_fasta
        all_contigs
    
    main:
        blast_makedb( ref_fasta, params.blastdbname)


}

 // workflow direct_blast_n () {
 //     take:
 //         ref_fasta
 //         all_contigs
 //     main:
 // 
 //         make_db_for_blast( ref_fasta, "TRUE") 
 //         do_blastn(all_contigs, make_db_for_blast.out.DB, params.subtax_dir)
 //         
 //         if (params.handle_contamination == true ) {
 //             handle_contamination_pr( params.cids, 
 //                                      params.cfaa, 
 //                                      do_blastn.out.OUT,
 //                                      do_blastn.out.CONTIGS )
 //                                      
 //             blastOut=handle_contamination_pr.out
 //         } else {
 //         
 //             blastOut=do_blastn.out.OUT
 //         }
 //         
 //         blast_sum_coverage(blastOut, 
 //                            "F", 
 //                            "F", 
 //                            params.refdb_dir, 
 //                            params.subasb_dir, 
 //                            params.bindir,
 //                            params.reports_dir )
 //         
 //     emit:
 //         BY_R=blast_sum_coverage.out.BYR
 //         BY_SQ=blast_sum_coverage.out.BYSQ
 //         BY_SP=blast_sum_coverage.out.BYSP
 //         S_SUM=blast_sum_coverage.out.SUM
 //         CFA=all_contigs
 //         DONE=blast_sum_coverage.out.STA
 // }

 // workflow direct_blast_tx () {
 //     take:
 //         ref_fasta
 //         all_contigs
 //     main:
 // 
 //         make_db_for_blast( ref_fasta, "FALSE") 
 //         do_tblastx(all_contigs, make_db_for_blast.out.DB, params.taxtbxdir)
 //         
 //         if (params.handle_contamination == true ) {
 //             handle_contamination_pr( params.cids, 
 //                                      params.cfaa, 
 //                                      do_tblastx.out.OUT,
 //                                      do_tblastx.out.CONTIGS )
 //                                      
 //             blastOut=handle_contamination_pr.out
 //         } else {
 //         
 //             blastOut=do_tblastx.out.OUT
 //         }
 //         
 //         blast_sum_coverage(blastOut, "F", "F" )
 //         
 //     emit:
 //         BY_R=blast_sum_coverage.out.BYR
 //         BY_SQ=blast_sum_coverage.out.BYSQ
 //         BY_SP=blast_sum_coverage.out.BYSP
 //         S_SUM=blast_sum_coverage.out.SUM
 //         CFA=all_contigs
 //         DONE=blast_sum_coverage.out.STA
 // }


workflow coverage_compute() {
    
    take:
        contigs_fa
        assign_byr
        reffa

    main:
     
     // reffa="${params.amplicon_refseqs_dir}/${params.amplicon_refseqs}"
     make_db_for_blast( reffa, "FALSE")   // FALSE refers to: not redo db if it already exists.
     do_cov_on_viralcandidates("${params.refdb_dir}/${params.set_tax}", 
                                make_db_for_blast.out.DB,
                                contigs_fa,
                                assign_byr,
                                params.bl_suffix,
                                params.bindir, 
                                params.tmp_dir, 
                                params.cov_dir
                                )

   
    emit:
      blastout=do_cov_on_viralcandidates.out.BLOUT
      coverage=do_cov_on_viralcandidates.out.COV

}

// // // // // // MAIN // // // // // //  


params.tmp_dir      =  "${params.results_dir}/tmp"      
// params.rawqc_dir    =  "${params.basedir}/raw"     
// params.clnfq_dir    =  "${params.basedir}/clean"   
params.ampaln_dir   =  "${params.results_dir}/aln"     
params.asbl_dir     =  "${params.results_dir}/assembly"
params.taxdir       =  "${params.results_dir}/taxonomy"
params.cov_dir      =   "${params.results_dir}/coverage"
// params.reports_dir  =  "${params.results_dir}/reports" 
params.logs_dir     =  "${params.results_dir}/logs"    
params.html_dir     =  "${params.ctvdir}/html"    

workflow () {
  
  
    //check_params()
    println "# Running   : $workflow.scriptId - $workflow.scriptName"
    println "# Project   : $workflow.projectDir"
    println "# Bdir      : $workflow.launchDir"
    println "# Starting  : $workflow.userName $ZERO $workflow.start"
    println "# Reading samples for $params.runID from $params.samp"

    println " ### $workflow.launchDir ## $params.bindir ## $params.refseqs ## $params.tmp_dir ## $params.html_dir ## $params.logs_dir"

    
    ch_samples = Channel
        .fromPath(params.samp)
        .splitCsv(header: false, sep: '\t')
        .filter { row -> !row[0].startsWith('#') }
        .map { row -> 
            def meta = [ id: row[0], illumina: row[1] ]

            def path1  = "${params.fastq_dir}/${row[1]}${params.R1}.${params.rawfq_sfx}"
            def f1_list     =  file(path1)  // creates a list with all combinations: fq.gz and fastq.gz.
            if( f1_list.isEmpty() ) exit 1, "MISSING FILE: ${f1}\nCheck if ILLUMINA_ID matches the filename!"
        
            def path2  = "${params.fastq_dir}/${row[1]}${params.R2}.${params.rawfq_sfx}"
            def f2_list     = file(path2)
            if( f2_list.isEmpty() ) exit 1, "MISSING FILE: ${f2}\nCheck if ILLUMINA_ID matches the filename!"
        
            return [ meta, [f1_list.first(), f2_list.first()]  ]
        }
        .view { meta, files -> "ID: ${meta.id} | Files: ${files.collect { it.name }}" }


    // Initial Setup: Infrastructure Folders
    [params.logs_dir, params.tmp_dir].each { dir ->
      if (dir) {
        def d = file(dir)
        if( !d.exists() ) {
            d.mkdirs()
            log.info "Created Infrastructure Dir: $dir"
          }
         }
      }

    // init_run(filesystem, create_logd.out) 
    
    // 1 // Clean reads // //
    
        //fastqc_onrawseqs(ch_samples)
        fastQC_raw(ch_samples) 
        sampsqual=fastQC_raw.out.zip
                 .map {it[1]}
                 .collect()
                 .view { "Sending to MultiQC: " + it.collect { file -> file.name }.join(', ') }
        multiQC_raw(sampsqual, "raw")

        if (params.trim_adapters == true ) {
            bbduk_clean(ch_samples) 
            fastQC_clean(bbduk_clean.out.cleanReads)
            clean_sampsqual=fastQC_raw.out.zip
                 .map {it[1]}
                 .collect()
                 .view { "Sending to MultiQC: " + it.collect { file -> file.name }.join(', ') }
            multiQC_clean(clean_sampsqual, "clean")

          }

        // 2 // Filter contamination (Discard reads identified as nonviral) // //

        if (params.run_fcsgx) {
              ch_fasta = fq2fasta (bbduk_clean.out.cleanReads)
              ch_fasta_for_fcs = ch_fasta.fasta.map { meta, fasta ->  
                      return [ meta, '10239',fasta]  
                    }

              FCSGX_RUNGX ( ch_fasta_for_fcs,  file(params.fcsgx_db),  []  )
              ch_to_filter = bbduk_clean.out.cleanReads.join( FCSGX_RUNGX.out.fcsgx_report )
              // ch_to_filter.view { item -> 
              //   meta  = item[0]
              //   reads = item[1]
              //   conta = item[2]
              //   return ">>> CHECK: Sample [${meta.id}] \n    - Reads: ${reads} \n    - Contaminants: ${conta}\n"
              // }
              discard_contaminants(ch_to_filter)
      
              thy_clean_reads=discard_contaminants.out.cleaned_reads
             // ch_to_assembly = discard_contaminants.out.cleaned_reads.map { meta, reads -> 
             //                   return [ meta, [reads[0], reads[1]], reads[2] ] 
             //                 }
        }else {
              thy_clean_reads=bbduk_clean.out.cleanReads
            //  ch_to_assembly = bbduk_clean.out.cleanReads.map { meta, reads -> 
            //                    return [ meta, [reads[0], reads[1]], reads[2] ] 
            //                  }
             
        }
        
        ch_to_assembly = thy_clean_reads.map { meta, reads -> 
                                return [ meta, [reads[0], reads[1]], reads[2] ] 
                              }
        ch_to_assembly.view { item -> 
                meta  = item[0]
                reads = item[1]
                return ">>> CHECK: Sample [${meta.id}] \n    - Reads: ${reads}\n"
              }
     
     if (params.assembler ==~ /(?i)MEGAHIT/){
              megahit_assembly_new(ch_to_assembly)
              ch_raw_assembly = megahit_assembly_new.out.contigs
      } else if (params.assembler ==~ /(?i)METASPADES/ ){
              metaspades_assembly_new(ch_to_assembly)
              ch_raw_assembly = metaspades_assembly_new.out.contigs
      }
      ch_to_refine=ch_raw_assembly
                    .join( thy_clean_reads )
                    .map{ meta, contigs, reads ->
                            def r1  = reads.find { it.name.contains('pe1') }
                            def r2  = reads.find { it.name.contains('pe2') }
                            def sgl = reads.find { it.name.contains('sgl') } ?: file("NO_SGL_FILE")
                         return [ meta, contigs, [r1, r2], sgl ]
                       }
     refine_assembly(ch_to_refine)

     ref_fasta=file("${params.refdb_dir}/${params.fams_subset}");
     blast_makedb( ref_fasta, params.blastdbname)
     blast_search(refine_assembly.out.final_assbly, blast_makedb.out.db.collect())
     taxonomy_database=file("${params.refdb_dir}/${params.full_tax}")
     ch_to_sum=blast_search.out.blOUT
                    .join(refine_assembly.out.idxstats)
                    .map{ meta, blastout, readstats ->
                            def pe=readstats.find { it.name.contains('pe') }
                            def se=readstats.find { it.name.contains('se') } ?: file("NO_SGL_FILE")
                            return [ meta, blastout, pe, se ]
                         }
                      
     blast_sum_coverage( ch_to_sum , taxonomy_database )
}

workflow future (){

    // 3 // Map reads on database   // //
        dbtobowtie="$params.refdb_dir/$params.fams_subset"
        generate_index_bowtie(dbtobowtie)
        reads_align_wf(generate_index_bowtie.out, reads_filter_nonviral.out)

    // 4 // Assembly reads into contigs // //
    
        if (params.assembler ==~ /(?i)MEGAHIT/){
             megahit_assembly_all(reads_filter_nonviral.out, 
                                  params.subasb_dir, 
                                  params.logs_dir  )
             CNFA=megahit_assembly_all.out.CGSout
        }else if (params.assembler ==~ /(?i)METASPADES/ ){
             metaspades_assembly(reads_filter_nonviral.out, 
                                 params.subasb_dir
                                 )
             CNFA=metaspades_assembly.out.CGSout
        }
    
    // 5 // Taxonomic classification of contigs // //
        
        if (params.taxalg ==~ /(?i)KAIJU/ ) {
            // 5.1. FAST APPROACH// 
              // 5.1.1. Protein-level classification (KAIJU) //
            KCDB=Channel.from(params.kaijudb)
            to_kaiju_contigs=KCDB.combine(CNFA.merge()).merge()
            kaiju_contigs(to_kaiju_contigs)
            kaiju_summarize(kaiju_contigs.out.NM)
                       
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
            ref_fa="${params.refdb_dir}/${params.fams_subset}";
            // blast_flow( ref_fa, CNFA.merge())

            // direct_blast(ref_fa, CNFA.merge())
            direct_blast_n(ref_fa, CNFA)
            
             if ( params.do_cov_figures == true) {
                  
                  ref_fasta="$params.refdb_dir/${params.set_seqs}"

                  coverage_compute(direct_blast_n.out.CFA, 
                                     direct_blast_n.out.BY_R, 
                                     ref_fasta
                                    )

                  blOUT=coverage_compute.out.blastout                   
                  covOUT=coverage_compute.out.coverage
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

            direct_blast_tx(ref_fa, CNFA)
            
            if ( params.do_cov_figures == true) {
                  coverage_compute(direct_blast_tx.out.CFA, 
                                     direct_blast_tx.out.BY_R
                     )
                  blOUT=coverage_compute.out.blastout                   
                  covOUT=coverage_compute.out.coverage
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
        
        
        sampdef="$params.samp"
        make_summary_tbl(FIN, params.subasb_dir, params.reports_dir, params.bindir, sampdef)
        
        // 6.2. // Plot coverage by genome (reads and contigs)
        
        vizualise_results_flow(
                       reads_align_wf.out.PE,
                       reads_align_wf.out.SG,
                       blOUT,
                       covOUT,
                       params.bindir,
                       params.reports_dir,
                       params.refdb_dir, 
                       params.gff_dir
                       )
        figsDIR=vizualise_results_flow.out
        
        // 6.3. // HTML report. 
        fill_html_report (make_summary_tbl.out, 
                          sampdef,
                          params.bindir,
                          figsDIR, 
                          params.reports_dir,
                          params.html_dir)
        
                       
        // 6.2 // Taxonomy summary
            // taxon_outkaiju  (*.rvdb.names.out file)
            // blast_unc_x_cl.RPT
            // blast_unc_x_cl.out.RPT
    
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the reports in your browser...\n" : "Oops .. something went wrong: ${workflow.errorMessage}" )
}
