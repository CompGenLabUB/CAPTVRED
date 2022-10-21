import sys
import numpy as np
import re
import pandas as pd


'''
USAGE: 

python3 /dir/to/virwaste_get_summary_tables.py     \
        /dir/to/sample_definitions.tbl \
        /dir/to/raw/multiqc_fastqc.txt \
        /dir/to/clean/multiqc_fastqc.txt \
        /dir/to/filter/multiqc_fastqc.txt  \
         /dir/to/align/multiqc_bowtie2.txt \
         /dir/to/outfile.csv

OUTPUT FILE:
Fields are:  IDSample, N.-Raw-Reads, N.-Clean-Reads-PE, N.-Clean-Reads-SG, N.-Clean-Reads-TOT, 
             N.-Filt-Reads-PE, N.-Filt-Reads-SG, N.-Filt-Reads-TOT, N.-Aligned-SG, N.-Aligned-PE-Both,
             N.-Aligned-PE-OneMate, N.-Aligned-Total
'''

samples={}
decoder={}
try:
        with open(sys.argv[1], 'r') as file:
                for line in file.readlines():
                        if not line.startswith("#"):
                                ln=line.rstrip().split("\t")
                                samples[ln[0]]= [int(0)]*11 #[ int(0), int(0), int(0), int(0), int(0), int(0), int(0), int(0), int(0), int(0)]
                                decoder[ln[1]]=ln[0]
                samples["TOTAL"] = [int(0)]*11
            
except ValueError:
        print("No samples file found... Please give samlpes tbl file.")
print(decoder)
print( samples )
print( "\n\n")
try:
        with open(sys.argv[2], 'r') as file:
                for line in file.readlines():
                        ln=line.rstrip().split("\t")
                        idf=re.sub('_R[1,2]_001.fastq.gz', "", ln[0])
                        if not idf=="Sample":
                                samples[decoder[idf]][0] += int(float(ln[4]))
                                samples["TOTAL"][0] += int(float(ln[4]))
except ValueError:
        print("No multiqc found for RAW READS... Table will show zeroes for all related fields.")
print( samples )
print( "\n\n")


try:
        with open(sys.argv[3], 'r') as file:
                for line in file.readlines():
                        ln=line.rstrip().split("\t")
                        idf=ln[0]
                        if not idf=="Sample":
                                if bool(re.search('pe[1,2]', idf)):
                                        idf=re.sub('_pe[1,2].fastq.gz', '', idf)
                                        samples[idf][1] += int(float(ln[4]))
                                        samples[idf][3] += int(float(ln[4]))
                                        samples["TOTAL"][1] += int(float(ln[4]))
                                        samples["TOTAL"][3] += int(float(ln[4]))
                                elif bool(re.search('sgl', idf)):
                                        idf=re.sub('_sgl.fastq.gz', "", idf)
                                        samples[idf][2] += int(float(ln[4]))
                                        samples[idf][3] += int(float(ln[4]))
                                        samples["TOTAL"][2] += int(float(ln[4]))
                                        samples["TOTAL"][3] += int(float(ln[4]))
except ValueError:
        print("No multiqc found for CLEAN READS... Table will show zeroes for all related fields.")
print( samples )
print( "\n\n")
                                
try:
        with open(sys.argv[4], 'r') as file:
                for line in file.readlines():
                        ln=line.rstrip().split("\t")
                        idf=ln[0]
                        if not idf=="Sample":
                                if bool(re.search('pe[1,2]', idf)):
                                        idf=re.sub('_pe[1,2].filtered.fastq.gz', '', idf)
                                        samples[idf][4] += int(float(ln[4]))
                                        samples[idf][6] += int(float(ln[4]))
                                        samples["TOTAL"][4] += int(float(ln[4]))
                                        samples["TOTAL"][6] += int(float(ln[4]))
                                elif bool(re.search('sgl', idf)):
                                        idf=re.sub('_sgl.filtered.fastq.gz', "", idf)
                                        samples[idf][5] += int(float(ln[4]))
                                        samples[idf][6] += int(float(ln[4]))
                                        samples["TOTAL"][5] += int(float(ln[4]))
                                        samples["TOTAL"][6] += int(float(ln[4]))
except ValueError:
        print("No multiqc found for FILTERED READS... Table will show zeroes for all related fields.")
print( samples )
print( "\n\n")

try:
        with open(sys.argv[5], 'r') as file:
                for line in file.readlines():
                        ln=line.rstrip().split("\t")
                        idf=ln[0]
                        if not idf=="Sample":
                                if bool(re.search('sg', idf)):
                                         idf=re.sub('_sg.bowtie.log', "", idf)
                                         samples[idf][7] += int(float(ln[16]) + float(ln[17]))    ## Sgl end -> one + multiple
                                         samples[idf][10] += int(float(ln[16]) + float(ln[17])) 
                                elif bool(re.search('pe', idf)):
                                        idf=re.sub('_pe.bowtie.log', '', idf)
                                        samples[idf][8] += int(float(ln[4]) + float(ln[5]))        ## pend both mates -> one + multiple hits
                                        samples[idf][10]  += int(float(ln[4]) + float(ln[5]))
                                        
                                        samples[idf][9] += int(float(ln[8]) + float(ln[9]))    ## pend one mate -> one + multiple hits
                                        samples[idf][10] += int(float(ln[8]) + float(ln[9]))
                                 
except ValueError:
        print("No multiqc found for BAM ALIGNED READS... Table will show zeroes for all related fields.")
print( samples )
print( "\n\n")

try:
        with open(sys.argv[6], 'r') as file:
                for line in file.readlines():
                        ln=line.rstrip().split("\t")
                        idf=ln[0]
                        samples[idf][11] += int(float(ln[1])
                        samples[idf][12] += int(float(ln[4])
except ValueError:
        print("No CONTIGS SUMMARY file... Table will show zeroes for all related fields.")
print( samples )
print( "\n\n")


#print(samples)
outd=open(sys.argv[7], 'w')
outd.write("#Sample_ID,Sample_Type,Sequencing_Method,N_RawReads,N_cleanReads_PE,N_CleanReads_SGE,")
outd.write("N_AlignedReads_PE_both,N_AlignedReads_PE_onemate,N_AlignedReads_SGLE,N_AlignedReads_Total\n")
for s in sorted(samples):
        #strg=[s, [i for i in samples[s]]]
        #print (strg)
        outd.write("{}".format(s))
        [outd.write(",{}".format(i)) for i in samples[s]]
        outd.write("\n")
outd.close()

                                
exit(0);

