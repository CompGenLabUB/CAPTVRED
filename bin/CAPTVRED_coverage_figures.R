args <- commandArgs(TRUE);


## READ CL inputs
GFF.refs.dir  <- args[1];
refs.coord.tbl <- args[2];
samp.id <- args[3];
BOWTIEpe.filename <- args[4];
BOWTIEsg.filename <- args[5];
contigs.gff.filename <- args[6]
id.to.name <- args[7]
outplots.dir <- args[8];

#GFF.refs.dir  <- "/data/virpand/pandemies/refseqs/ampliconseqs/gff_refgenomes"
#refs.coord.tbl <- "/data/virpand/pandemies/refseqs/ampliconseqs/gff_refgenomes/refseqs_coordinates.tbl"  
#samp.id <- "TS_SD_00"
#BOWTIEpe.filename <- "/data/virpand/pandemies/TEST_SET/DEV_TESTSET/amplicons_alignment/TS_SD_00_pe.bowtie.sorted.mapped.bam"
#BOWTIEsg.filename <- "/data/virpand/pandemies/TEST_SET/DEV_TESTSET/amplicons_alignment/TS_SD_00_sg.bowtie.sorted.mapped.bam" 
#contigs.gff.filename <- "/data/virpand/pandemies/TEST_SET/DEV_TESTSET/taxonomy/taxon_viral_candidates/TS_SD_00/TS_SD_00_blastn_on_viralcandidates.gff"
#id.to.name <- "/data/virpand/pandemies/refseqs/ampliconseqs/Viral_candidates_zoonosis.ids"
#outplots.dir <-  "/data/virpand/pandemies/TEST_SET/DEV_TESTSET/reports/coverage_figures"


## Upload required libraries
library(rtracklayer);
library(GenomicFeatures);
library(Rsamtools);
library(GenomicAlignments);
library(VariantAnnotation);
library(ggplot2);
library(ggbio);
library(plyr);
library(gridExtra);

## Viz variables
colors.set <- c("grey",
                "#d96d58",  # red
                "#619cdb",  # blue
                "#d295db",  # purple
                "#84d579",  # green
                "#cec75b"); # yellow
                
pointSize   <-  4;
textSize    <- 10;
spaceLegend <-  0.5;

ref.coor<-read.table(refs.coord.tbl, row.names=1, col.names=c("ID", "START", "END"))

## Read data info:
  #Read PE bam
pe.bam <- readGAlignments(BOWTIEpe.filename, # PE
                        use.names=TRUE)     # param=ScanBamParam(which=GRanges(BOWTIEpe.seq, IRanges(Xmin, Xmax))));
  #Reads SE bam
sg.bam <- readGAlignments(BOWTIEsg.filename, # SG
                        use.names=TRUE)    # param=ScanBamParam(which=GRanges(BOWTIEsg.seq, IRanges(Xmin, Xmax))));
  #Contigs from gff:

print(file.info(contigs.gff.filename)$size)
ContigsFound=FALSE;
if (file.info(contigs.gff.filename)$size > 0) {
        contigs.txdb <- makeTxDbFromGFF(file=contigs.gff.filename,
                                    format="gff3");
        ContigsFound=TRUE;
}
print(contigs.txdb$user_seqlevels)

#NA Plots:
NA_plot <- function(xmax, ymax, message) {
                theplot <- ggplot() +
                        theme_bw() +
                        geom_text(aes(xmax/2, ymax/2, label=message),
                                    size=10, color='grey', alpha=0.6) +
                        xlab(NULL) + ylab(NULL) +
                        xlim(0, xmax) +
                        ylim(0,ymax) +
                        scale_x_continuous(expand=c(0,0)) +
                        theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.6),
                              #panel.border = element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    text = element_text(size = 10),
                                    plot.title = element_text(size = 15)) 
                return(theplot)
          }
    
## Convert genomes ID to Human Readable name
genomes_rel <- read.table(id.to.name, header=TRUE, sep="\t", comment.char="")

print(genomes_rel$SpecTaxonId)
for (txn in genomes_rel$SpecTaxonId) {
     print("taxon is")
     print(txn)
     txvec<-as.vector(genomes_rel[genomes_rel$SpecTaxonId==txn,]$SeqID)
     plot_list <- list()
     for (rgn in txvec) {
        print("...")
        print(rgn)
        GFF.filename  <- paste0(GFF.refs.dir, "/", rgn ,".gff")
        rgn.v <- rgn
        rgn.id <- unlist(strsplit(rgn, '[.]'))[1];
        rgn.name <-  genomes_rel[genomes_rel$SeqID==rgn,]$Name
        
      if (!file.exists(GFF.filename)) {
            rgn <- rgn.id;
            GFF.filename  <- paste0(GFF.refs.dir, "/", rgn ,".gff")
        }
        
        print("preparint txdb object")
        txdb <- makeTxDbFromGFF(file=GFF.filename,
                            format="gff3");
        # TxDb object
        ntx <- length(levels(as.factor(as.list(txdb)$transcripts$tx_id)));
        colors.genes <- ifelse(ntx > 1, hcl.colors(ntx, palette="Spectral"), "blue");
        # parameters
        Xmin  <- as.integer(ref.coor[rgn.id, 1,]);
        Xmax  <- as.integer(ref.coor[rgn.id, 2,]); # WE NEED genome length
        Xmax <- round_any(Xmax, 500, f = ceiling)
        
        #plot
        print("Ploting refgen")
        refgen <- autoplot(txdb,
                           which=GRanges(rgn.id, IRanges(Xmin, Xmax)),
                                         names.expr = "gene_id", fill=colors.genes) +
                      theme_bw() +
                      theme(panel.border = element_blank(),
                            axis.text.y=element_blank(),
                            text = element_text(size = 10),
                            plot.title = element_text(size = 15)) +
                      scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0));
    # Paired ends coverage
        
        #subset bam
        pesub.bam<-pe.bam[seqnames(pe.bam) == rgn.v, ]     ## SI ESTA EMPTY FER UN PLOT EN BLANC
        pe.avgcov <- mean(coverage(pesub.bam));
        
        # plot
        print("Ploting PE_COV")
        if (length( pesub.bam) > 0) {
            pe.coverage <- autoplot(pesub.bam , geom = "line", stat = "coverage") +
                               theme_bw() + 
                               theme(text = element_text(size = 10)) + ylab("") +
                               scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0)) 
        }else{
            pe.coverage <- NA_plot( Xmax, 1, 'No alignments found') 
        }
    #  Single ends coverage
        
        #subset bam
        sgsub.bam<-sg.bam[seqnames(sg.bam) == rgn.v, ]        ## SI ESTA EMPTY FER UN PLOT EN BLANC
        sg.avgcov <- mean(coverage(sgsub.bam));
        
        #plot
        print("Ploting SG_COV")
        if (length( sgsub.bam) > 0) {
            sg.coverage <- autoplot(sgsub.bam , geom = "line", stat = "coverage") +
              theme_bw() + theme(text = element_text(size = 10)) + ylab("") +
              scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0)) 
        }else{
            sg.coverage <- NA_plot( Xmax, 1, 'No alignments found')
        }

        #plot
        print("Ploting Assembled contigs")
        print(ContigsFound)
        if ( ContigsFound) {
                print(" --- AAA --- ")
	        print(rgn.v)
                print(contigs.txdb$user_seqlevels)
                if (any(contigs.txdb$user_seqlevels == rgn.v)) {
#          if ( ContigsFound & any(contigs.txdb$user_seqlevels == rgn.v) ) {
                    print(" --- BBB --- ")
                    colors.genes <- "lightseagreen" #hcl.colors(1, palette="Spectral");
                    contigs <- autoplot(contigs.txdb, which=GRanges(rgn.v, IRanges(Xmin, Xmax)),
                                names.expr = "transcript", fill=colors.genes) +
                                theme_bw() +
                                theme(panel.border = element_blank(),
                                    axis.text.y=element_blank(),
                                    text = element_text(size = 10),
                                    plot.title = element_text(size = 15)) +
                                scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0));
                } else {
                    print(" --- bbb --- ")
                    contigs <- NA_plot( Xmax, 1, 'No contigs > 100nt assembled')+
                                scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0)) +
                                theme(axis.text.x=element_text())
                }
        } else {
                    print(" --- CCC --- ")
                    contigs <- NA_plot( Xmax, 1, 'No contigs > 100nt assembled')+
                                scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0)) +
                                theme(axis.text.x=element_text())
        }
       #  Merge all plots
          print("Ready to merge the plots")
          title=paste0( rgn.name, " (", rgn.v, ") \n Sample=", samp.id);
          wholeplot <- tracks(
                        CDS=refgen,
                        PEcovg=pe.coverage,
                        SEcovg=sg.coverage,
                        Assembly=contigs,
                        heights = c(0.2, 0.2, 0.2, 0.2),
                        main.height=3,
                        xlim = c(Xmin, Xmax), 
                        main=title
                    )  +
                 scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0));
      PNG.filename<- paste0(outplots.dir, "/Coverage_", samp.id, "_onto_", rgn.v, ".png") 
      
      #  Save plot
      ggsave(PNG.filename, plot=wholeplot,
              width = 25, height = 15, units = "cm", dpi = 600);
    };
};






