#!/bin/bash

if [ $# -lt 1 ] ; then
    echo "";
    echo "usage: source rbbh-blast2gff.sh  [BRHfile] [OUTBLASTfilename]  [outGFF]";
    echo "Blast outfile expected fields: qseqid, qlen, sseqid, slen,
          qstart, qend, sstart, send, length, score, evalue, 
          bitscore, pident, nident, mismatch, positive, gapopen,
          gaps, ppos, qframe, sframe, qcovs, qcovhsp, qseq, sseq";
          
    echo "Best reciprocal hit tabular file (obtained from rhhb.py)"
    exit 0;
fi;

BRH=$1;
echo "BRH file: $BRH";
BLO=$2;
echo "blast outfile: $BLO";
#GFFout=$(echo $BLO | sed 's/.tbl$/.gff/' - );
GFFout=$3;
echo "Output gff file $GFFout";

[[ -f $GFFout ]] && rm $GFFout;
touch  $GFFout;

idslist=$(awk '!( $1 ~ /^#/) { print $2}' $BRH | sort | uniq);
refgnlist=$(awk '!( $1 ~ /^#/) { print $1}' $BRH | sort | uniq);


for i in $idslist; 
    do for j in $refgnlist; 
       do
          awk -v id=$i -v rg=$j ' 
            ($1 ~ id && $3 ~ rg) {
                    if($7>$8) { 
                        S=$8; E=$7; T="-"; 
                    } else { 
                        S=$7; E=$8; T="+"; 
                    };
                    print $3"\tblastn\ttranscript\t"S"\t"E"\t"$13"\t"T"\t.\tID="$1
                }' $BLO >> $GFFout;
           
       done;
    done;


