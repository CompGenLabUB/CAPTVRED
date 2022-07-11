#! /usr/bin/env nextflow
//mychan=Chanel.fromPath("params.rawfq/samplesMap.keySet()*.fastq.gz")
//params.rawfq = "$params.wdir/rawseqs_fastq/samplesMap.keySet()_{R1,R2}*.fastq.gz"
//params.rawfqD = "$params.wdir/rawseqs_fastq/"
//params.clnseqsD="$params.wdir/cleanseqs"

// process bbduk_clean2   {
//
//     samplesMap.each { key, value ->
//          if (reads[0].contains("${key}") && reads[1].contains("$key")) {               //if infile contains $key 
//          //def $out1 = reads[0].replace("$key", "$value");   //set $value to outfile
//          //def $out2 = reads[1].replace("$key", "$value");   //set $value to outfile
//          def inR1 = reads[0].split("$key")[0] + value + '_pe1.fastq.gz'
//          def inR1 = reads[1].split("$key")[0] + value + '_pe1.fastq.gz'
//          def out1 = reads[0].split("$key")[0] + value + '_pe1.fastq.gz'
//          def out2 = reads[1].split("$key")[0] + value + '_pe2.fastq.gz'
//          //def $outSary= reads[0].split("$key")
//          def outS = reads[0].split("$key")[0] + value + '_se.fastq.gz'
//          
//          println "Sample: $key   Fastq root: $value"
//          println "1: $out1 \n 2: $out2 \n S: $outS"
//          
//          }
//     }
//
//     input:
//      //tuple val(old), val(new) from samplesMap
//     //path raw from params.rawfqD
//     //path x
//     // tuple val(pair_id), file(reads)
//      path inR1
//      path inR2
//      
//    output:
//     //tuple val("${out1}") val("${out2}") val("${outS}")
//     path out1
//     path out2
//     path outS
//
//
//     script:
//
//     //shell:
//     '''
//     # Adapters Trimming:
//     cat <<"EOF" 
//     $BBDUK_PATH/bbduk.sh  \
//             in=$inR1 \
//         in2=$inR2 \
//         out=$out1   \
//         out2=$out2   \
//         outs=$outS   \
//            ref=$BBDUK_PATH/resources/adapters.fa \
//            k=13 ktrim=r useshortkmers=t mink=5 \
//            qtrim=t trimq=20 minlength=32           \
//        threads=!{NCPUS} overwrite=true maq=10    \
//             2> $outfile.log 1>&2 ;
//    touch $outfile.bbduk_clean.ok;
//EOF    
//     '''
//
//}

process bbduk_clean {
	input:
	  val idir
	  val odir
      val samples
	 
	output:
	  // path "${odir}/${id_sample}_pe1.fastq.gz", emit: outPE1
	  // path "${odir}/${id_sample}_pe2.fastq.gz", emit: outPE2
	  // path "${odir}/${id_sample}_sgl.fastq.gz", emit: outSGL
	  stdout
	
	script:
	samples.each { id_sample, id_illumina -> 

		println "# BBDUK on $idir/$id_sample (from $id_illumina)"
				"""
		mkdir -vp ${odir};
		touch ${odir}/prova2.txt
		cat "# Adapters Trimming :";
#		cat <<EOF 
#		$BBDUK_PATH/bbduk.sh  \
#				 in=!{idir}/!{id_illumina}_R1_001.fastq.gz \
#				in2=!{idir}/!{id_illumina}_R2_001.fastq.gz \
#				out=!{odir}/!{id_sample}_pe1.fastq.gz \
#			   out2=!{odir}/!{id_sample}_pe2.fastq.gz \
#			   outs=!{odir}/!{id_sample}_sgl.fastq.gz \
#				ref=$BBDUK_PATH/resources/adapters.fa \
#				  k=13 ktrim=r useshortkmers=t mink=5 \
#			  qtrim=t trimq=20 minlength=32           \
#			threads=!{param.NCPUS} overwrite=true maq=10 \
#				 2> !{odir}/!{id_sample}.bbduk_clean.log 1>&2 ;
#		touch !{odir}/!{id_sample}.bbduk_clean.ok;
#	EOF    
		"""

     println "   outfile: $odir/$id_sample ... "
   } 
   // samples.each

}

process bbduk_clean2 {
	input:
	  tuple val(sampid), val(illuid)
	 
	output:
	  val "${illuid}_pe1.fastq.gz", emit: outPE1
	  val "${illuid}_pe2.fastq.gz", emit: outPE2
	  val "${illuid}_sgl.fastq.gz", emit: outSGL
	  
	  //tuple path ("${illuid}_pe1.fastq.gz"),
	    //    path ("${illuid}_pe2.fastq.gz"),
	      //  path ("${illuid}_sgl.fastq.gz)"
	
	  
	
	script:
	
	//println "# BBDUK on $sampid"
	
	"""
	mkdir -vp $params.clnfq_dir;
	# Adapters Trimming:
	
	\$BBDUK_PATH/bbduk.sh  \
	         in=${sampid}_R1_001.fastq.gz \
		    in2=${sampid}_R2_001.fastq.gz \
		    out=${illuid}_pe1.fastq.gz \
		   out2=${illuid}_pe2.fastq.gz \
		   outs=${illuid}_sgl.fastq.gz \
            ref=\$BBDUK_PATH/resources/adapters.fa \
              k=13 ktrim=r useshortkmers=t mink=5 \
          qtrim=t trimq=20 minlength=32           \
        threads=$params.NCPUS overwrite=true maq=10 \
             2> ${illuid}.bbduk_clean.log 1>&2 ;
    touch ${illuid}.bbduk_clean.ok; 
	"""

}


