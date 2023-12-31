#! /usr/bin/env nextflow

process generate_index_bowtie {
    input:

         //val dep
         val x

    output:
    
        val idx_path, emit: IDX

    script:

        idx_path=x.toString().replaceAll(".fa.gz", "_index")
        idx_files=file("${idx_path}.*.bt2")

       
        //bowtie2-buid generates AT LEAST 6 files, if less than 6 are found, db might be corrupted
        if( idx_files.size() >= 6 )
        
            """
            cat ${params.amplicon_refseqs_dir}/bowtie_genome_indexed.cdate 1>&2;
            
            """
            
        else
        
            """
             bowtie2-build --threads $params.NCPUS -f     \
                          $x     ${idx_path}              \
                          2> ${idx_path}.bowtiedb.log;
             date +"DB created on %Y/%m/%d %T %Z %s"      \
                  > ${params.amplicon_refseqs_dir}/bowtie_genome_indexed.cdate;
            """
}


process bowtie_amplicons_alignment {

    input:
      val(idx_path)
      val(pe1)
      val(pe2)
      val(sgle)
    
    output:
    
      val "${odir}/${pe_root}.bowtie.sorted.mapped.bam",   emit: bamPE
      
    script:

                      
        //idir=params.clnfq_dir
        odir=params.ampaln_dir

        pe_root=pe1.split('/')[-1].toString().replaceAll("_pe1.filtered.fastq.gz", "_pe")

        """
        
        bowtie2 --threads $params.NCPUS -q  -x ${idx_path} \
                --end-to-end --sensitive                   \
                -N ${params.bowtie_nmismatch}              \
                -L ${params.bowtie_seedlen}                \
                --met-file ${odir}/${pe_root}.metrics      \
                -1 ${pe1} -2 ${pe2}                        \
                -S ${odir}/${pe_root}.bowtie.sam           \
                2> ${odir}/${pe_root}.bowtie.log;
        
        samtools view -Sb -o ${odir}/${pe_root}.bowtie.bam    \
                        ${odir}/${pe_root}.bowtie.sam;
        samtools sort ${odir}/${pe_root}.bowtie.bam           \
                       > ${odir}/${pe_root}.bowtie.sorted.bam;
        rm -vf ${odir}/${pe_root}.bowtie.bam;
        samtools index  ${odir}/${pe_root}.bowtie.sorted.bam;
        samtools view -F2052  -q ${params.alignMINQ} -h ${odir}/${pe_root}.bowtie.sorted.bam \
                    -o ${odir}/${pe_root}.bowtie.sorted.mapped.bam
        touch ${odir}/${pe_root}.ok
        
        """
    
}

process bowtie_amplicons_alignment_sg {

    input:
      val(idx_path)
      val(pe1)
      val(pe2)
      val(sgle)
    
    output:

      val "${odir}/${sgle_root}.bowtie.sorted.mapped.bam", emit: bamSG
      
    script:

        //idir=params.clnfq_dir
        odir=params.ampaln_dir
        
        sgle_file=file("${sgle}")
        sgle_root=sgle.split('/')[-1].toString().replaceAll("_sgl.filtered.fastq.gz", "_sg")
        
        if( sgle_file.size() > 0 )
            
            """
            bowtie2 --threads $params.NCPUS -q  -x ${idx_path}   \
                    --end-to-end --sensitive   -N 1              \
                    --met-file ${odir}/${sgle_root}.metrics      \
                    -U  ${sgle}                                  \
                    -S ${odir}/${sgle_root}.bowtie.sam           \
                    2> ${odir}/${sgle_root}.bowtie.log;
                    
            samtools view -Sb  -o ${odir}/${sgle_root}.bowtie.bam    \
                          ${odir}/${sgle_root}.bowtie.sam;
            samtools sort ${odir}/${sgle_root}.bowtie.bam           \
                          > ${odir}/${sgle_root}.bowtie.sorted.bam;
            rm -vf ${odir}/${sgle_root}.bowtie.bam;
            samtools index  ${odir}/${sgle_root}.bowtie.sorted.bam;
            samtools view -F2052 -h -q ${params.alignMINQ} ${odir}/${sgle_root}.bowtie.sorted.bam \
                  -o ${odir}/${sgle_root}.bowtie.sorted.mapped.bam;
            touch ${odir}/${sgle_root}.ok
            """

        else
            """
            touch ${odir}/${sgle_root}.bowtie.sam ${odir}/${sgle_root}.bowtie.bam ${odir}/${sgle_root}.bowtie.sorted.bam;
            """

}
