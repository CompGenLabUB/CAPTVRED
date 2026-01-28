#! /usr/bin/env nextflow


process megahit_assembly_new{

    input:
      tuple val(meta), path(pe), path(single)
    
    output:
    tuple val(meta), path("megahit_out/${meta.id}.contigs.fa"), emit: contigs
    path "*.log"
      
    script:
    
      //  sp_root=sgle.split('/')[-1].replaceAll("_sgl.filtered.fastq.gz", "")
      // odir="$mhdir/${sp_root}"
      //  outfl="${odir}/${sp_root}.contigs+singletons.fa"
        out_dir="megahit_out"
      """
        rm -rf ${out_dir}
        megahit  -t ${task.cpus}  --presets meta-large     \
            -1 ${pe[0]} -2 ${pe[1]}  -r ${single}          \
            --min-contig-len ${params.assemblyMINCONLEN}   \
            --out-dir ${out_dir}  --out-prefix ${meta.id}  \
            2> ${meta.id}.megahit_assembly.log;
      """
}

process metaspades_assembly_new {
    tag "$meta.id"
    input:
      tuple val(meta), path(pe), path(single)

    output:
    tuple val(meta), path("metaspades_out/${meta.id}.scaffolds.fasta"), emit: contigs
    path "*.log"
      
      
    script:
        odir="metaspades_out"
       """
         mkdir -vp ${odir};
          
         ## 1 ## Assembly reads into contigs:
         spades    --meta                            \
         -t ${task.NCPUSNCPUS}  -m ${task.memory.toGiga()}  \
         -1 ${pe[0]}  -2 ${pe[1]}                    \
         --phred-offset ${params.phred}              \
         -o ${odir}  2> ${meta.id}.metaspades_assembly.log 1>&2;
         
      """
}


process refine_assembly {
    tag "$meta.id"
  
    input:
    tuple val(meta), path(contigs), path(pe), path(single)

    output:
    tuple val(meta), path("${meta.id}.contigs+singletons.fa"), emit: final_assbly
    tuple val(meta), path("*.stats.txt"), emit: idxstats
    path "*.log"

    script:
    // def prefix = meta.id
    """
    # 1. Verificamos si hay contigs válidos (Megahit a veces saca archivos vacíos)
    if [ -s "${contigs}" ]; then

        cat ${contigs} > ${meta.id}.contigs+singletons.fa
        ## 1. ## Create index (contigs as db):
        bowtie2-build --threads $task.cpus ${contigs} ${meta.id}_index

        ## 2. Map reads back to contigs to extract unmapped reads ##
        ### 2.1. Single-End reads:  ###
        if [ -s  "${single}" ]; then
           ( bowtie2 -p $task.cpus         \
                    -x ${meta.id}_index   \
                    -U ${single}          \
                    --end-to-end --sensitive -N 1 \
                    --met-file ${meta.id}.bowtie_onto_contigs.se.metrics.txt |\
              samtools view -bS - | samtools sort -o ${meta.id}_se.bam -
            
            # Estadísticas
            samtools index ${meta.id}_se.bam
            samtools idxstats ${meta.id}_se.bam > ${meta.id}_se.stats.txt
            
            # Extraer unmapped y añadir al assembly final
            samtools fasta -f4 -F2048 ${meta.id}_se.bam | \
            seqkit seq -m ${params.assemblyMINCONLEN} >> ${meta.id}.contigs+singletons.fa 
            ) 2> ${meta.id}_map_se_to_contigs.log 
        fi

        ### 2.2. Paired-End reads:  ###
         ( bowtie2 -p $task.cpus           \
                  -x ${meta.id}_index     \
                  -1 ${pe[0]} -2 ${pe[1]} \
                  --end-to-end --sensitive -N 1 \
                  --met-file ${meta.id}.bowtie_onto_contigs.pe.metrics.txt  |\
            samtools view -bS - | samtools sort -o ${meta.id}_pe.bam -

        # Estadísticas
        samtools index ${meta.id}_pe.bam
        samtools idxstats ${meta.id}_pe.bam > ${meta.id}_pe.stats.txt

        # Extract unmapped and add to final assembly
        samtools fasta -f4 -F2048 ${meta.id}_pe.bam | \
        seqkit seq -m ${params.assemblyMINCONLEN} | \
        sed -E 's/^>([^ ]+)[/ ]([12]).*/>\\1:R\\2 \\2/' >> ${meta.id}.contigs+singletons.fa 
        ) 2> ${meta.id}_map_pe_to_contigs.log 

    else
        ## if no contigs obtained -> keep al reads > MinimumContigsLength nt
        # Procesar SE
        if [ -s "${single}" ]; then
            seqkit seq -m ${params.assemblyMINCONLEN} ${single} >> ${meta.id}.contigs+singletons.fa
        fi
        
        # Procesar PE
        seqkit seq -m ${params.assemblyMINCONLEN} ${pe[0]} >> ${meta.id}.contigs+singletons.fa
        seqkit seq -m ${params.assemblyMINCONLEN} ${pe[1]} >> ${meta.id}.contigs+singletons.fa
        
        #mock files to output
        cat "no contigs available to align with" > ${meta.id}.metrics.txt 
        cat "no contigs available to align with" > ${meta.id}_non_mapped.log 
    fi
    """
}

 // process megahit_assembly_all {
 // 
 //     input:
 //       tuple val(meta), path(pe), path(single)
 //     
 //     output:
 //       tuple val(meta), val(outfl)  , emit : CGSout
 // 
 //       
 //     script:
 //     
 //       //  sp_root=sgle.split('/')[-1].replaceAll("_sgl.filtered.fastq.gz", "")
 //       // odir="$mhdir/${sp_root}"
 //       //  outfl="${odir}/${sp_root}.contigs+singletons.fa"
 //         
 //         """
 //         [ -d ${odir} ] && rm -r ${odir};
 //         ## 1 ## Assembly reads into contigs:
 //         megahit  -t ${params.NCPUS}  --presets meta-large     \
 //             -1 ${pe[0]} -2 ${pe[1]}  -r ${single}                   \
 //             --min-contig-len ${params.assemblyMINCONLEN}      \
 //             --out-dir ${odir} --out-prefix ${sp_root}         \
 //             2> ${logd}/${sp_root}.megahit_assembly.log;
 //         
 //          cat  ${odir}/${sp_root}.contigs.fa  >  ${odir}/${sp_root}.contigs+singletons.fa;
 //         
 //         ## 2A ## If at least 1 contig has been assembled map reads on the contigs.
 //         if [ -s ${odir}/${sp_root}.contigs.fa ]; then
 //                
 //            ## 2.1. ## Create index (contigs as db):
 //                  bowtie2-build --threads $params.NCPUS    -f     \
 //                                ${odir}/${sp_root}.contigs.fa     \
 //                                ${odir}/${sp_root}_index          \
 //                                   2> ${odir}/${sp_root}_index.bowtiedb.log;
 //            
 //            ## 2.2. ## Signle-End reads:
 //                 ## a. map reads:
 //                   bowtie2 --threads $params.NCPUS -q  -x ${odir}/${sp_root}_index   \
 //                             --end-to-end --sensitive    -N 1             \
 //                              --met-file ${odir}/${sp_root}.bowtie_onto_contigs.metrics   \
 //                             -U  ${sgle}    \
 //                             -S ${odir}/${sp_root}_se.bowtie_onto_contigs.sam  \
 //                             2> ${odir}/${sp_root}_se.bowtie_onto_contigs.log;
 //                 
 //                 ## b. Get fasta of unmapped:
 //                   samtools fasta -f4 -F2048  ${odir}/${sp_root}_se.bowtie_onto_contigs.sam  |\
 //                     seqkit seq -m ${params.assemblyMINCONLEN} >> ${odir}/${sp_root}.contigs+singletons.fa;
 //                 
 //                 ## c. Get idxstats of mapped reads:
 //                  ( 
 //                   samtools view -Sb -f2 -F2052  \
 //                                      ${odir}/${sp_root}_se.bowtie_onto_contigs.sam         \
 //                                    > ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.bam;
 //                   samtools sort      ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.bam    \
 //                                    > ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.bam;
 //                   samtools index     ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.bam;
 //                   samtools idxstats  ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.bam   \
 //                                    > ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.stats;
 //                 ) > ${odir}/${sp_root}_se.get_maped_reads_info.log 2>&1;
 //             
 //             ## 2.3. ## Paired-End reads:
 //                 ## a. map reads:
 //                   bowtie2 --threads $params.NCPUS -q  -x ${odir}/${sp_root}_index     \
 //                          --end-to-end --sensitive   -N 1                             \
 //                          --met-file ${odir}/${sp_root}.bowtie_onto_contigs.metrics   \
 //                          -1  ${pe1}  -2  ${pe2}                                      \
 //                          -S ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam            \
 //                          2> ${odir}/${sp_root}_pe.bowtie_onto_contigs.log;
 //                          
 //                 ## b. Get fasta of unmapped:
 //                   samtools fasta -f4 -F2048  ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam  |\
 //                     seqkit seq -m ${params.assemblyMINCONLEN}                        |\
 //                     sed 's/^>\\([^ ]*\\) \\([12]\\)/>\\1:R\2 \\2/'                   |\
 //                     >>  ${odir}/${sp_root}.contigs+singletons.fa;
 //                 ## c. Get idxstats of mapped reads:
 //                  ( 
 //                   samtools view -Sb -f2 -F2052 \
 //                                      ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam         \
 //                                    > ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.bam;
 //                   samtools sort      ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.bam    \
 //                                    > ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.bam;
 //                   samtools index     ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.bam;
 //                   samtools idxstats  ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.bam   \
 //                                    > ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.stats;
 //                 ) > ${odir}/${sp_root}_pe.get_maped_reads_info.log 2>&1;
 //                 
 //         else
 //         ## 2B ## if no contigs obtained -> keep al reads > MinimumContigsLength nt
 //             for fl in $pe1 $pe2 $sgle; do
 //                seqkit seq -g -m ${params.assemblyMINCONLEN}  \${fl}  |\
 //                   seqkit fq2fa - \
 //                   >> ${odir}/${sp_root}.contigs+singletons.fa;
 //             done;
 //             touch ${odir}/${sp_root}_se.bowtie_onto_contigs.maped.sorted.stats;
 //             touch ${odir}/${sp_root}_pe.bowtie_onto_contigs.maped.sorted.stats;
 //             
 //         fi
 //            
 //         """
 // 
 // }


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
         
         
         
         
         
