#! /usr/bin/env nextflow


process megahit_assembly_all {

    input:
      val (pe1)
      val (pe2)
      val (sgle)
      val (mhdir)
      val (logd)
    
    output:
      val  outfl , emit : CGSout

      
    script:
    
        sp_root=sgle.split('/')[-1].replaceAll("_sgl.filtered.fastq.gz", "")
        odir="$mhdir/${sp_root}"
        outfl="${odir}/${sp_root}.contigs+singletons.fa"
        
        """
        [ -d ${odir} ] && rm -r ${odir};
        ## 1 ## Assembly reads into contigs:
        megahit  -t ${params.NCPUS}  --presets meta-large     \
            -1 ${pe1} -2 ${pe2}  -r ${sgle}                   \
            --min-contig-len ${params.assemblyMINCONLEN}      \
            --out-dir ${odir} --out-prefix ${sp_root}         \
            2> ${logd}/${sp_root}.megahit_assembly.log;
        
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
                  samtools fasta -f4 -F2048  ${odir}/${sp_root}_se.bowtie_onto_contigs.sam  |\
                    seqkit seq -m ${params.assemblyMINCONLEN} >> ${odir}/${sp_root}.contigs+singletons.fa;
                
                ## c. Get idxstats of mapped reads:
                 ( 
                  samtools view -Sb -f2 -F2052  \
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
                  samtools fasta -f4 -F2048  ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam  |\
                    seqkit seq -m ${params.assemblyMINCONLEN}                        |\
                    sed 's/^>\\([^ ]*\\) \\([12]\\)/>\\1:R\2 \\2/'                   |\
                    >>  ${odir}/${sp_root}.contigs+singletons.fa;
                ## c. Get idxstats of mapped reads:
                 ( 
                  samtools view -Sb -f2 -F2052 \
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



process metaspades_assembly {

    input:
      val (pe1)
      val (pe2)
      val (sgle)
      val (spdir)
    
    output:
      val  outfl , emit : CGSout

      
    script:
        sp_root=sgle.split('/')[-1].replaceAll("_sgl.filtered.fastq.gz", "")
        odir="$spdir/${sp_root}"
        outfl="${odir}/${sp_root}.contigs+singletons.fa"

       """
          mkdir -vp ${odir};
          
         ## 1 ## Assembly reads into contigs:
         spades    --meta                         \
         -t ${params.NCPUS}  -m ${params.MAXMEM}  \
         -1 ${pe1}  -2 ${pe2}                     \
         --phred-offset ${params.phred}           \
         -o $odir  2> ${odir}.metaspades_assembly.log 1>&2;
         
        cat ${odir}/scaffolds.fasta >  ${outfl};
         
          ## 2A ## If at least 1 scaffold has been assembled map PE reads on the scaffolds.
          
    if [ -s ${odir}/scaffolds.fasta ]; then 
          
          ## 2.1. Create index: 
          echo "$sp_root building idx" >> /data/virpand/pandemies/TEST_SET/SIMDATA/kkk;
              bowtie2-build --threads $params.NCPUS    -f            \
                               ${odir}/scaffolds.fasta     \
                               ${odir}/${sp_root}_index              \
                                  2> ${odir}/${sp_root}_index.bowtiedb.log;
          ## 2.2. ## Signle-End reads:
            echo "$sp_root let'sgo SG" >> /data/virpand/pandemies/TEST_SET/SIMDATA/kkk;
            sz=\$(zcat ${sgle} | wc -l);
            echo "\$sz"  >> /data/virpand/pandemies/TEST_SET/SIMDATA/kkk;
           if [ \$sz -ne  0 ]; then
                ## a. map reads:
                 echo "$sp_root SG  COND IS TRUE" >> /data/virpand/pandemies/TEST_SET/SIMDATA/kkk;
                  bowtie2 --threads $params.NCPUS -q  -x ${odir}/${sp_root}_index   \
                            --end-to-end --sensitive    -N 1             \
                             --met-file ${odir}/${sp_root}.bowtie_onto_scaffolds.metrics     \
                            -U  ${sgle}                          \
                            -S ${odir}/${sp_root}_se.bowtie_onto_scaffolds.sam           \
                            2> ${odir}/${sp_root}_se.bowtie_onto_scaffolds.log;
                
            
                ## b. Get fasta of unmapped:
                  samtools fasta -f4  ${odir}/${sp_root}_se.bowtie_onto_scaffolds.sam  |\
                    seqkit seq -m ${params.assemblyMINCONLEN} >> ${odir}/${sp_root}.contigs+singletons.fa;
                
                ## c. Get idxstats of mapped reads:
                 ( 
                  samtools view -Sb -G 4  \
                                     ${odir}/${sp_root}_se.bowtie_onto_scaffolds.sam         \
                                   > ${odir}/${sp_root}_se.bowtie_onto_scaffolds.maped.bam;
                  samtools sort      ${odir}/${sp_root}_se.bowtie_onto_scaffolds.maped.bam    \
                                   > ${odir}/${sp_root}_se.bowtie_onto_scaffolds.maped.sorted.bam;
                  samtools index     ${odir}/${sp_root}_se.bowtie_onto_scaffolds.maped.sorted.bam;
                  samtools idxstats  ${odir}/${sp_root}_se.bowtie_onto_scaffolds.maped.sorted.bam   \
                                   > ${odir}/${sp_root}_se.bowtie_onto_scaffolds.maped.sorted.stats;
                ) > ${odir}/${sp_root}_se.get_maped_reads_info.log 2>&1;
           else 
                 echo "$sp_root block ignoredddd" >> kkk;
           fi;
            echo "$sp_root end of SG" >> kkk;
            
             ## 2.3. ## Paired-End reads:
             echo "$sp_root let'sgo PE" >> kkk;
             
                    echo "$sp_root PE  COND IS TRUE" >> kkk;
                    ## a. map reads:
                      bowtie2 --threads $params.NCPUS -q  -x ${odir}/${sp_root}_index      \
                             --end-to-end --sensitive   -N 1                               \
                             --met-file ${odir}/${sp_root}.bowtie_onto_scaffolds.metrics   \
                             -1  ${pe1}  -2  ${pe2}                                        \
                             -S ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.sam            \
                             2> ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.log;
                             
                    ## b. Get fasta of unmapped:
                      samtools fasta -f4  ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.sam  |\
                        seqkit seq -m ${params.assemblyMINCONLEN}                          |\
                        sed 's/^>\\([^ ]*\\) \\([12]\\)/>\\1:R\2 \\2/'                      \
                        >>  ${outfl} 2> ${odir}/${sp_root}_pe.getunmapedfa.log ;
                    
                    ## c. Get idxstats of mapped reads:
                     ( 
                      samtools view -Sb -G 4  \
                                         ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.sam         \
                                       > ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.maped.bam;
                      samtools sort      ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.maped.bam    \
                                       > ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.maped.sorted.bam;
                      samtools index     ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.maped.sorted.bam;
                      samtools idxstats  ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.maped.sorted.bam   \
                                       > ${odir}/${sp_root}_pe.bowtie_onto_scaffolds.maped.sorted.stats;
                    ) > ${odir}/${sp_root}_pe.get_maped_reads_info.log 2>&1;
                
        else
                    ## 2B ## if no contigs obtained -> keep all reads > MinimumContigsLength nt
                        for fl in $pe1 $pe2 $sgle; do
                           seqkit seq -g -m ${params.assemblyMINCONLEN}  \${fl}  |\
                              seqkit fq2fa -                                      \
                              >> ${outfl};
             done;
             touch ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.stats;
             touch ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.stats;
            
          fi;
     
      """
}
         
         
         
         
         
