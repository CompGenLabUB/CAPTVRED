#Please, note that this is a TAB-separated file, where:
##
#SAMPLE_ID(*) - desired prefix to be used for all the files created in the downstream analyses (it can be the same as sequencing_id). 
#SEQUENCING_ID(*) - prefix that appears in the raw fastq files outputed from the basecalling.
#SAMPLE_NAME(*) - description of the biological content of the sample. Please try to aviod spaces.
#PAIRED (boolean) -  use YES or NO. (Note/This version only considers paired data).
#SEQUENCING_METHOD - technology, targetted aproach,... Please try to aviod spaces.
#SAMPLE_ORIGIN - Where the sample was collected,... Please try to aviod spaces.
#METADATA - any other information of interest. Please try to aviod spaces.
# Fields marked with "(*)" are required for the CAPTVRED analysis. 
#Other fields are not necessary (but might be used in future versions od CAPTVRED protocol)
##
#SAMPLE_ID	SEQUENCING_ID	SAMPLE_NAME	PAIRED	SEQUENCING_METHOD	SAMPLE_ORIGIN   METADATA
example_00  EX00    An_example_sample   YES Illumina    River_water_locationX   