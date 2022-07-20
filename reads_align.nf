#! /usr/bin/env nextflow

process generate_index_bowtie {
    input:

         val dep
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
      val dep
      tuple val(idx_path), val(pe1), val(pe2), val(sgle)
    
    output:
    
      val "${odir}/${pe_root}.bowtie.sorted.bam",   emit: bamPE
      
    script:

                      
        idir=params.clnfq_dir
        odir=params.ampaln_dir

        pe_root=pe1.toString().replaceAll("_pe1.fastq.gz", "_pe")

        
        
        """
        
        bowtie2 --threads $params.NCPUS -q  -x ${idx_path} \
                --local --sensitive-local  -N 1            \
                --met-file ${odir}/${pe_root}.metrics      \
                -1 ${idir}/${pe1} -2 ${idir}/${pe2}        \
                -S ${odir}/${pe_root}.bowtie.sam           \
                2> ${odir}/${pe_root}.bowtie.log;
        
        samtools view -Sb -o ${odir}/${pe_root}.bowtie.bam    \
                        ${odir}/${pe_root}.bowtie.sam;
        samtools sort ${odir}/${pe_root}.bowtie.bam           \
                       > ${odir}/${pe_root}.bowtie.sorted.bam;
        rm -vf ${odir}/${pe_root}.bowtie.bam;
        samtools index  ${odir}/${pe_root}.bowtie.sorted.bam;
        samtools view -F4 -f67 -h ${odir}/${pe_root}.bowtie.sorted.bam \
                    -o ${odir}/${pe_root}.bowtie.sorted.mapped.bam
        touch ${odir}/${pe_root}.ok
        
        """
    
}

process bowtie_amplicons_alignment_sg {

    input:
      val dep
      tuple val(idx_path), val(pe1), val(pe2), val(sgle)
    
    output:

      val "${odir}/${sgle_root}.bowtie.sorted.bam", emit: bamSG
      
    script:

        idir=params.clnfq_dir
        odir=params.ampaln_dir
        
        sgle_file=file("${idir}/${sgle}")
        sgle_root=sgle.toString().replaceAll("_sgl.fastq.gz", "_sg")
        
        if( sgle_file.size() > 0 )
            
            """
            bowtie2 --threads $params.NCPUS -q  -x ${idx_path}   \
                    --local --sensitive-local   -N 1             \
                     --met-file ${odir}/${sgle_root}.metrics     \
                    -U  ${idir}/${sgle}                          \
                    -S ${odir}/${sgle_root}.bowtie.sam           \
                    2> ${odir}/${sgle_root}.bowtie.log;
                    
            samtools view -Sb -o ${odir}/${sgle_root}.bowtie.bam    \
                          ${odir}/${sgle_root}.bowtie.sam;
            samtools sort ${odir}/${sgle_root}.bowtie.bam           \
                          > ${odir}/${sgle_root}.bowtie.sorted.bam;
            rm -vf ${odir}/${sgle_root}.bowtie.bam;
            samtools index  ${odir}/${sgle_root}.bowtie.sorted.bam;
            samtools view -F2052 -h  ${odir}/${sgle_root}.bowtie.sorted.bam \
                  -o ${odir}/${sgle_root}.bowtie.sorted.mapped.bam;
            touch ${odir}/${sgle_root}.ok
            """

        else
            """
            touch ${odir}/${sgle_root}.bowtie.sam ${odir}/${sgle_root}.bowtie.bam ${odir}/${sgle_root}.bowtie.sorted.bam;
            """

}
