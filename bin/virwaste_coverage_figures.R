

args <- commandArgs(TRUE);

GFF.refs.dir  <- args[1];  # Directory where gff of regerence sequences are stored
refs.coord.tbl <- args[2];  # tbl containing fields: ID,  START_POS,  END_POS
samp.id <- args[3];
BOWTIEpe.filename <- args[4];  
BOWTIEsg.filename <- args[5];
outplots.dir <- args[6];



#GFF.refs.dir  <- "refseqs"; #args[4];
#refs.coot.tbl <- "refseqs/refseqs_coordinates.tbl"
#BOWTIEpe.filename <- "prova_C10/R01_C10_G02_pe.bowtie.sorted.mapped.bam"
#BOWTIEsg.filename <- "prova_C10/R01_C10_G02_sg.bowtie.sorted.mapped.bam"
#GFF.seqid <- "NC_002645" # !! sense la versio (.1, .2, .3, ...)
#Xmin  <- as.integer(2000);
#Xmax  <- as.integer(27317); # genome length
#   Xmax <- round_any(Xmax, 500, f = ceiling)
#BOWTIEpe.seq <- "NC_002645"
#BOWTIEsg.seq <- "NC_002645"



#library(tidyverse)
library(rtracklayer);
library(GenomicFeatures);
library(Rsamtools);
library(GenomicAlignments);
library(VariantAnnotation);
library(ggplot2);
library(ggbio);
library(plyr);

# taxdb <- loadTaxonomyDb()
# taxdb[2697049, ] <- c(2697049, 'Betacoronavirus', 'SARS-CoV-2')
# 1928434 -> Betacoronavirus sp.
colors.set <- c("grey",
                "#d96d58",  # red
                "#619cdb",  # blue
                "#d295db",  # purple
                "#84d579",  # green
                "#cec75b"); # yellow
                
pointSize   <-  4;
textSize    <- 10;
spaceLegend <-  0.5;


#GFF.fields <- c("SeqID", "Source", "Feat", "Start", "End", "Score", "Strand", "Frame", "Info");
#thygap <- read.table(GFF.filename);
#colnames(thygap) <- GFF.fields

#txdb <- makeTxDbFromGFF(file=GFF.filename,
#                        format="gff3"); 
             
                        

 #OrYel  /  #

#Read genomes coordinates:
ref.coor<-read.table(refs.coord.tbl, row.names=1, col.names=c("ID", "START", "END"))

#Read PE bam
pe.bam <- readGAlignments(BOWTIEpe.filename, # PE
                        use.names=TRUE)
                        
                # param=ScanBamParam(which=GRanges(BOWTIEpe.seq, IRanges(Xmin, Xmax))));

#Reads SE bam
sg.bam <- readGAlignments(BOWTIEsg.filename, # SG
                        use.names=TRUE)
                        
                # param=ScanBamParam(which=GRanges(BOWTIEsg.seq, IRanges(Xmin, Xmax))));

#Get genomes to which there are maped reads
pe.refs<-as.data.frame(table(seqnames(pe.bam)))
sg.refs<-as.data.frame(table(seqnames(sg.bam)))
reflst <- union(sg.refs[sg.refs$Freq!=0, 1], pe.refs[pe.refs$Freq!=0, 1])

#NA Plots:
NA_plot <- function(xmax, ymax, message) {
                theplot <- ggplot() +
                        theme_void() +
                        geom_text(aes(xmax/2, ymax/2, label=message),
                                    size=10, color='grey', alpha=0.6) +
                        xlab(NULL) +
                        xlim(0, xmax) +
                        ylim(0,ymax) +
                        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6)) 
                return(theplot)
          }


for (rgn in reflst) {  ##REPENSAR, VAL LA PENA... mehhh
    #Reference genome
        
        #filename
        
        GFF.filename  <- paste0(GFF.refs.dir, "/", rgn ,".gff")
        rgn.v <- rgn
        rgn.id <- unlist(strsplit(rgn, '[.]'))[1];
        
        if (!file.exists(GFF.filename)) {
            rgn <- rgn.id;
            GFF.filename  <- paste0(GFF.refs.dir, "/", rgn ,".gff")
        }
        
        txdb <- makeTxDbFromGFF(file=GFF.filename,
                            format="gff3");
                            
        print(paste0("rgn: ", rgn, "  //  rgn.v: ", rgn.v, "  //  rgn.id: ", rgn.id))
        
        

        # TxDb object
        
        colors.genes <- hcl.colors(length(levels(as.factor(as.list(txdb)$transcripts$tx_id))), palette="Spectral")
        # parameters
        Xmin  <- as.integer(ref.coor[rgn.id, 1,]);
        Xmax  <- as.integer(ref.coor[rgn.id, 2,]); # WE NEED genome length
        Xmax <- round_any(Xmax, 500, f = ceiling)
        
        #plot
        print("Ploting refgen")
        refgen <- autoplot(txdb, which=GRanges(rgn.id, IRanges(Xmin, Xmax)),
               names.expr = "gene_id", fill=colors.genes) +
        theme_bw() +
        theme(panel.border = element_blank(),
            axis.text.y=element_blank(),
            text = element_text(size = 10),
            plot.title = element_text(size = 15)) +
        scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0));
    
    # Paired ends coverage
        
        #subset bam
        pesub.bam<-pe.bam[seqnames(pe.bam) == rgn.v ]     ## SI ESTA EMPTY FER UN PLOT EN BLANC
        pe.avgcov <- mean(coverage(pesub.bam));
        
        # plot
        print("Ploting PE_COV")
        if (length( pesub.bam) > 0) {
            pe.coverage <- autoplot(pesub.bam , geom = "line", stat = "coverage") +
                theme_bw() + theme(text = element_text(size = 10)) + ylab("") +
                scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0)) #+
                #geom_hline(yintercept=pe.avgcov, color="red", linetype="dotdash"); 
        }else{
            pe.coverage <- NA_plot( Xmax, 1, 'No alignments found')
        }
    
    
    #  Single ends coverage
        
        #subset bam
        sgsub.bam<-sg.bam[seqnames(sg.bam) == rgn.v ]        ## SI ESTA EMPTY FER UN PLOT EN BLANC
        sg.avgcov <- mean(coverage(sgsub.bam));
        
        #plot
        print("Ploting SG_COV")
        if (length( sgsub.bam) > 0) {    
            sg.coverage <- autoplot(sgsub.bam , geom = "line", stat = "coverage") +
              theme_bw() + theme(text = element_text(size = 10)) + ylab("") +
              scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0)) #+
              #geom_hline(yintercept=sg.avgcov, color="red", linetype="dotdash"); 
        }else{
            sg.coverage <- NA_plot( Xmax, 1, 'No alignments found')
        }
        
    #Pending: Contigs
    
    #  Merge all plots
        print("Ready to merge the plots")
        title=paste0( "Reads and amplicons coverage from sample ", samp.id, " onto ", rgn.v);
        wholeplot <- tracks(CDS=refgen, PEcovg=pe.coverage, SEcovg=sg.coverage ,
               heights = c(0.15, 0.3, 0.3), xlim = c(Xmin, Xmax), main=title) +
               scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0));
    PNG.filename<- paste0(outplots.dir, "/Coverage_", samp.id, "_onto_", rgn.v, ".png") 
    #  Save plot
        ggsave(PNG.filename, plot=wholeplot,
               width = 25, height = 15, units = "cm", dpi = 600);
}
