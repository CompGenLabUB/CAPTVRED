#! /usr/bin/env nextflow

//process megahit_assembly_pe {
//    input:
//     val(pe1)
//     val(pe2)
//    
//    output:
//     val pe_root
//      
//    script:
//        // idir="$params.clnfq_dir"
//        odir="$params.asbl_dir/megahit"
//        pe_root=pe1.split('/')[-1].toString().replaceAll("_pe1.fastq.gz", "_pe")
//        //wholepath=odir.concat("/".concat(pe_root)) 
//        """
//        megahit  -t ${params.NCPUS}  --presets meta-sensitive\
//            -1 ${pe1} -2 ${pe2}  \
//            --out-dir ${odir}/${pe_root}  --out-prefix ${pe_root};
//        """
//}


//process megahit_assembly_sg {
//
//    input:
//     val(sgle)
//    
//    output:
//    val sgle_root
//      
//    script:
//        //idir="$params.clnfq_dir"
//        sgle_root=sgle.split('/')[-1].replaceAll("_sgl.fastq.gz", "_sg")
//        odir="$params.asbl_dir/megahit/${sgle_root}"
//        """
//        #[-d ${odir} ] && rm -r ${odir}
//        [-d ${odir} ] && echo " !!! ---- ${odir} ALREADY EXISTS ---- !!! "
//        #mkdir ${odir}
//        megahit  -t ${params.NCPUS}  --presets meta-sensitive \
//            -r ${sgle}   \
//            --out-dir ${odir} --out-prefix ${sgle_root}; 
//        """
//}



process megahit_assembly_all {
    
    
    
    input:
      val(pe1)
      val(pe2)
      val(sgle)
    
    output:
      val  outfl , emit : CGSout

      
    script:
    
        sp_root=sgle.split('/')[-1].replaceAll("_sgl.filtered.fastq.gz", "")
        odir="$params.asbl_dir/megahit/${sp_root}"
        outfl="${odir}/${sp_root}.contigs+singletons.fa"
        
        """
        [ -d ${odir} ] && rm -r ${odir};
        ## 1 ## Assembly reads into contigs:
        megahit  -t ${params.NCPUS}  --presets meta-large \
            -1 ${pe1} -2 ${pe2}  -r ${sgle}                   \
            --min-contig-len ${params.assemblyMINCONLEN}      \
            --out-dir ${odir} --out-prefix ${sp_root}         \
            2> ${params.logs_dir}/${sp_root}.megahit_assembly.log;
        
         cat  ${odir}/${sp_root}.contigs.fa  >  ${odir}/${sp_root}.contigs+singletons.fa;
        
        ## 2A ## If at least 1 contig has been assembled map reads on the contigs.
        if [ -s ${odir}/${sp_root}.contigs.fa ]; then
               
           ## 2.1. ## Create index (contigs as db):
                 bowtie2-build --threads $params.NCPUS    -f     \
                               ${odir}/${sp_root}.contigs.fa     \
                               ${odir}/${sp_root}_index          \
                                  2> ${odir}/${sp_root}_index.bowtiedb.log;
           
           ## 2.2. ## Signle-End reads:
                ## a. map reads:
                  bowtie2 --threads $params.NCPUS -q  -x ${odir}/${sp_root}_index   \
                            --end-to-end --sensitive    -N 1             \
                             --met-file ${odir}/${sp_root}.bowtie_onto_contigs.metrics     \
                            -U  ${sgle}                          \
                            -S ${odir}/${sp_root}_se.bowtie_onto_contigs.sam           \
                            2> ${odir}/${sp_root}_se.bowtie_onto_contigs.log;
                
                ## b. Get fasta of unmapped:
                  samtools fasta -f4  ${odir}/${sp_root}_se.bowtie_onto_contigs.sam  |\
                    seqkit seq -m ${params.assemblyMINCONLEN} >> ${odir}/${sp_root}.contigs+singletons.fa;
                
                ## c. Get idxstats of mapped reads:
                 ( 
                  samtools view -Sb -G 4  \
                                     ${odir}/${sp_root}_se.bowtie_onto_contigs.sam         \
                                   > ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.bam;
                  samtools sort      ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.bam    \
                                   > ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.bam;
                  samtools index     ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.bam;
                  samtools idxstats  ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.bam   \
                                   > ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.stats;
                ) > ${odir}/${sp_root}_se.get_maped_reads_info.log 2>&1;
                
            ## 2.3. ## Paired-End reads:
                ## a. map reads:
                  bowtie2 --threads $params.NCPUS -q  -x ${odir}/${sp_root}_index     \
                         --end-to-end --sensitive   -N 1                             \
                         --met-file ${odir}/${sp_root}.bowtie_onto_contigs.metrics   \
                         -1  ${pe1}  -2  ${pe2}                                      \
                         -S ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam            \
                         2> ${odir}/${sp_root}_pe.bowtie_onto_contigs.log;
                         
                ## b. Get fasta of unmapped:
                  samtools fasta -f4  ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam  |\
                    seqkit seq -m ${params.assemblyMINCONLEN}                        |\
                    sed 's/^>\\([^ ]*\\) \\([12]\\)/>\\1:R\2 \\2/'                   |\
                    >>  ${odir}/${sp_root}.contigs+singletons.fa;
                ## c. Get idxstats of mapped reads:
                 ( 
                  samtools view -Sb -G 4  \
                                     ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam         \
                                   > ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.bam;
                  samtools sort      ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.bam    \
                                   > ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.bam;
                  samtools index     ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.bam;
                  samtools idxstats  ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.bam   \
                                   > ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.stats;
                ) > ${odir}/${sp_root}_pe.get_maped_reads_info.log 2>&1;
                
        else
        ## 2B ## if no contigs obtained -> keep al reads > MinimumContigsLength nt
            for fl in $pe1 $pe2 $sgle; do
               seqkit seq -g -m ${params.assemblyMINCONLEN}  \${fl}  |\
                  seqkit fq2fa - \
                  >> ${odir}/${sp_root}.contigs+singletons.fa;
            done;
            touch ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.stats;
            touch ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.stats;
            
        fi
           
        """
        
        
        
}





process trinity_assembly_pe {

    label 'trinity_crash'
    
    input:
     val(pe1)
     val(pe2)
     val(sgle)
     
    output:
     val  "${odir}/${pe_root}.contigs+singletons.fa"
      
    script:
        pe_root=pe1.split('/')[-1].toString().replaceAll("_pe1.fastq.gz", "_pe")
        odir="$params.asbl_dir/trinity/${pe_root}_trinity"
        
        
        """
        [ -d ${odir} ] && rm -r ${odir}
        mkdir -vp ${odir}
        Trinity --seqType fq --CPU ${params.NCPUS} \
          --max_memory ${params.trinityMAXM} \
          --min_contig_length ${params.assemblyMINCONLEN}  \
          --left ${pe1}  --right ${pe2}  \
          --output ${odir}    \
          2> ${odir}/${pe_root}_trinity.log 1>&2;
          
          
          bowtie2-build --threads $params.NCPUS    -f     \
                       ${odir}/Trinity.fasta     \
                       ${odir}/Trinity_fa_index          \
                          2> ${odir}/Trinity_fa_index.bowtiedb.log;
        
         bowtie2 --threads $params.NCPUS -q  -x ${odir}/Trinity_fa_index  \
                    --end-to-end --fast    -N 1             \
                     --met-file ${odir}/${pe_root}.bowtie_pe_onto_contigs.metrics     \
                    -U  ${sgle}                          \
                    -S ${odir}/${pe_root}.bowtie_pe_onto_contigs.sam           \
                    2> ${odir}/${pe_root}.bowtie_pe_onto_contigs.log;

         awk 'length(\$10) >= ${params.assemblyMINCONLEN}' ${odir}/${pe_root}.bowtie_pe_onto_contigs.sam |\
           samtools fasta -f4 >>  ${odir}/${pe_root}.contigs+singletons.fa
        
        
         bowtie2 --threads $params.NCPUS -q  -x ${odir}/Trinity_fa_index    \
                 --end-to-end --fast    -N 1             \
                 -met-file ${odir}/${pe_root}.bowtie_onto_contigs.metrics     \
                 -m1  ${pe1}  -m2  ${pe2}                  \
                 -S ${odir}/${pe_root}.bowtie_se_onto_contigs.sam           \
                 2> ${odir}/${pe_root}.bowtie_se_onto_contigs.log;
          
          awk 'length(\$10) >= ${params.assemblyMINCONLEN}' ${odir}/${pe_root}.bowtie_se_onto_contigs.sam |\
           samtools fasta -f4 >>  ${odir}/${pe_root}.contigs+singletons.fa
          
        """
}


process trinity_assembly_sg {
    label 'trinity_crash'
     
    input:
      val(sgle)
    
    output:
      val sgle_root
      
    script:
    
        sgle_root=sgle.split('/')[-1].replaceAll("_sgl.fastq.gz", "_sg")
        odir="$params.asbl_dir/trinity/${sgle_root}_trinity"
        
        
        """
        [ -d ${odir} ] && rm -r ${odir}
        mkdir -vp ${odir}
        Trinity --seqType fq --CPU ${params.NCPUS} \
          --max_memory ${params.trinityMAXM} --NO_SEQTK \
          --min_contig_length 100  \
          --single ${sgle}  \
          --output ${odir}  \
          2> ${odir}/${sgle_root}_trinity.log 1>&2;
        """
}



