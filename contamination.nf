#! /usr/bin/env nextflow

process handle_contamination_pr() {
    input:
	val (rstd)  // List of refseqs to discard (NCBI IDS)
	val (fatod) // Protein fasta of seqs to discard (potential contamination).
	val (BLTOF) // Output of the BLAST that will be filtered
	val (CONTIGSFA) //contigfs fasta
	
    output:
	env (OUT)
	
    
    script:
	idstof=BLTOF.toString().replaceAll(".tbl", ".potentialcontamin.ids")
	fatof=BLTOF.toString().replaceAll(".tbl", ".potentialcontamin.fa")
	fatoflog=BLTOF.toString().replaceAll(".tbl", ".potentialcontamin.log")
	
	bldbo=fatod.toString().replaceAll(".faa", "_blastdb") //blast fb name
	
	bxo=BLTOF.toString().replaceAll(".tbl", ".blastxoncont.tbl")
	bxolog=BLTOF.toString().replaceAll(".tbl", ".blastxoncont.log")
	bxcov=bxo.toString().replaceAll(".tbl", ".cov.tbl")
	bxcovlog=bxo.toString().replaceAll(".tbl", ".cov.log")
	
	cotod=bxo.toString().replaceAll(".tbl", ".contigstodiscard.ids")
	BLFLT=BLTOF.toString().replaceAll(".tbl", ".filt.tbl")
        """
        #1. Filter contig ids that mapped refseqs potentially discarded:
        grep -f ${rstd}  ${BLTOF} | awk '{print \$1}'  >  ${idstof}
        if [ -s ${idstof} ]; then 
	    #2.Extract fasta of this sequences:
	    seqkit grep  --pattern-file  ${idstof} ${CONTIGSFA}   \
		> ${fatof}  2> ${fatoflog};
	    
	    #3. Create blastdb of the fasta to discard:
	    makeblastdb -in ${fatod}  -dbtype prot  -out ${bldbo};
	    
	    #4. Blastx:
	    blastx -query ${fatof} -db ${bldbo}  -out  ${bxo}            \
		   -outfmt \"${params.bl_outfmt}\"  -num_threads ${params.NCPUS} \
		   2> ${bxolog} 1>&2;
	    
	    #5. Compute QueryCoverage:
	     ${params.bindir}/coverage_blastshorttbl.pl \
		-prog=BLASTX  ${bxo} > ${bxcov}  2> ${bxcovlog};
	    
	    #6. Discard cov > threshold.
	     awk '(\$9> ${params.cont_min_cov}) && (\$1=="Q"){print \$2}' ${bxcov} > ${cotod};
	    grep -vwf ${cotod} ${BLTOF} > ${BLFLT};
	    OUT=${BLFLT};
	else
	    OUT=${BLTOF};
	fi;
        """
    
}
