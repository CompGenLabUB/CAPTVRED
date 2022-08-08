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
      val "" , emit: DEP
      
    script:
    
        sp_root=sgle.split('/')[-1].replaceAll("_sgl.filtered.fastq.gz", "")
        odir="$params.asbl_dir/megahit/${sp_root}"
        outfl="${odir}/${sp_root}.contigs+singletons.fa"
        """
        [ -d ${odir} ] && rm -r ${odir};

        megahit  -t ${params.NCPUS}  --presets meta-sensitive \
            -1 ${pe1} -2 ${pe2}  -r ${sgle}                   \
            --min-contig-len 100                              \
            --out-dir ${odir} --out-prefix ${sp_root}; 
        
        cat  ${odir}/${sp_root}.contigs.fa > ${odir}/${sp_root}.contigs+singletons.fa;
        
         bowtie2-build --threads $params.NCPUS    -f     \
                       ${odir}/${sp_root}.contigs.fa     \
                       ${odir}/${sp_root}_index          \
                          2> ${odir}/${sp_root}_index.bowtiedb.log;
        
         bowtie2 --threads $params.NCPUS -q  -x ${odir}/${sp_root}_index   \
                    --local --sensitive-local   -N 1             \
                     --met-file ${odir}/${sp_root}.bowtie_onto_contigs.metrics     \
                    -U  ${sgle}                          \
                    -S ${odir}/${sp_root}_se.bowtie_onto_contigs.sam           \
                    2> ${odir}/${sp_root}_se.bowtie_onto_contigs.log;

         awk 'NF < 4 || length(\$10) > 99' ${odir}/${sp_root}_se.bowtie_onto_contigs.sam |\
           samtools fasta -f4 - >>  ${odir}/${sp_root}.contigs+singletons.fa;
        
        
         bowtie2 --threads $params.NCPUS -q  -x ${odir}/${sp_root}_index   \
                 --local --sensitive-local   -N 1             \
                 --met-file ${odir}/${sp_root}.bowtie_onto_contigs.metrics     \
                 -1  ${pe1}  -2  ${pe2}                  \
                 -S ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam           \
                 2> ${odir}/${sp_root}_pe.bowtie_onto_contigs.log;
          
          awk 'NF < 4 || length(\$10) > 99' ${odir}/${sp_root}_pe.bowtie_onto_contigs.sam |\
           samtools fasta -f4 - >>  ${odir}/${sp_root}.contigs+singletons.fa;
           
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
          --max_memory 250G \
          --min_contig_length 100  \
          --left ${pe1}  --right ${pe2}  \
          --output ${odir}    \
          2> ${odir}/${pe_root}_trinity.log 1>&2;
          
          bowtie2-build --threads $params.NCPUS    -f     \
                       ${odir}/Trinity.fasta     \
                       ${odir}/Trinity_fa_index          \
                          2> ${odir}/Trinity_fa_index.bowtiedb.log;
        
         bowtie2 --threads $params.NCPUS -q  -x ${odir}/Trinity_fa_index  \
                    --local --sensitive-local   -N 1             \
                     --met-file ${odir}/${pe_root}.bowtie_pe_onto_contigs.metrics     \
                    -U  ${sgle}                          \
                    -S ${odir}/${pe_root}.bowtie_pe_onto_contigs.sam           \
                    2> ${odir}/${pe_root}.bowtie_pe_onto_contigs.log;

         awk 'length(\$10) > 99' ${odir}/${pe_root}.bowtie_pe_onto_contigs.sam |\
           samtools fasta -f4 >>  ${odir}/${pe_root}.contigs+singletons.fa
        
        
         bowtie2 --threads $params.NCPUS -q  -x ${odir}/Trinity_fa_index    \
                 --local --sensitive-local   -N 1             \
                 -met-file ${odir}/${pe_root}.bowtie_onto_contigs.metrics     \
                 -m1  ${pe1}  -m2  ${pe2}                  \
                 -S ${odir}/${pe_root}.bowtie_se_onto_contigs.sam           \
                 2> ${odir}/${pe_root}.bowtie_se_onto_contigs.log;
          
          awk 'length(\$10) > 99' ${odir}/${pe_root}.bowtie_se_onto_contigs.sam |\
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
          --max_memory 250G --NO_SEQTK \
          --min_contig_length 100  \
          --single ${sgle}  \
          --output ${odir}  \
          2> ${odir}/${sgle_root}_trinity.log 1>&2;
        """
}



