#! /usr/bin/env nextflow

process bbduk_clean {
	input:
	  //file x
	  tuple val(old), val(new) from samplesMap
	 
	output:
	 
	 
	shell:
	'''
	odir="${PATHTO}/cleanseqs";
	mkdir -vp $odir;
	infile="${PATHFROM}/!{old}";
	outfile="${odir}/!{new}";
	# Adapters Trimming:
	cat <<"EOF" 
	$BBDUK_PATH/bbduk.sh  \
	         in=${infile}_R1_001.fastq.gz \
		    in2=${infile}_R2_001.fastq.gz \
		    out=${outfile}_pe1.fastq.gz   \
		   out2=${outfile}_pe2.fastq.gz   \
		   outs=${outfile}_sgl.fastq.gz   \
            ref=$BBDUK_PATH/resources/adapters.fa \
              k=13 ktrim=r useshortkmers=t mink=5 \
          qtrim=t trimq=20 minlength=32           \
        threads=!{NCPUS} overwrite=true maq=10    \
             2> $outfile.log 1>&2 ;
    touch $outfile.bbduk_clean.ok;
EOF    
	'''

}
