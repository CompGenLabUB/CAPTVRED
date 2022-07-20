#! /usr/bin/env nextflow

process megahit_assembly_pe {

    input:
     val(pe1)
     val(pe2)
    
    output:
     val pe_root
      
    script:
        // idir="$params.clnfq_dir"
        odir="$params.asbl_dir/megahit"
        pe_root=pe1.split('/')[-1].toString().replaceAll("_pe1.fastq.gz", "_pe")
        //wholepath=odir.concat("/".concat(pe_root)) 
        """
        megahit  -t ${params.NCPUS}  --presets meta-sensitive\
            -1 ${pe1} -2 ${pe2}  \
            --out-dir ${odir}/${pe_root}  --out-prefix ${pe_root};
        """
}


process megahit_assembly_sg {

    input:
     val(sgle)
    
    output:
    val sgle_root
      
    script:
        //idir="$params.clnfq_dir"
        sgle_root=sgle.split('/')[-1].replaceAll("_sgl.fastq.gz", "_sg")
        odir="$params.asbl_dir/megahit/${sgle_root}"

        """
        #[-d ${odir} ] && rm -r ${odir}
        [-d ${odir} ] && echo " !!! ---- ${odir} ALREADY EXISTS ---- !!! "
        #mkdir ${odir}
        megahit  -t ${params.NCPUS}  --presets meta-sensitive \
            -r ${sgle}   \
            --out-dir ${odir} --out-prefix ${sgle_root}; 
        """
}



process megahit_assembly_all {

    input:
    val(pe1)
    val(pe2)
    val(sgle)
    
    output:
    val sp_root
      
    script:
    
        sp_root=sgle.split('/')[-1].replaceAll("_sgl.fastq.gz", "")
        odir="$params.asbl_dir/megahit/${sp_root}"
        """
        [ -d ${odir} ] && rm -r ${odir}
        megahit  -t ${params.NCPUS}  --presets meta-sensitive \
            -1 ${pe1} -2 ${pe2}  -r ${sgle}   \
            --out-dir ${odir} --out-prefix ${sp_root}; 
        """
}





process trinity_assembly_pe {

    label 'trinity_crash'
    
    input:
     val(pe1)
     val(pe2)
    
    output:
     val pe_root
      
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



