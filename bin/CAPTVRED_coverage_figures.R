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

if (FALSE) {
GFF.refs.dir  <- "/data/virpand/pandemies/CAPTVRED/references/db/rvdb_nt/gff_refgenomes"
refs.coord.tbl <- "/data/virpand/pandemies/CAPTVRED/references/db/rvdb_nt/gff_refgenomes/refseqs_coordinates.tsv"  
samp.id <- "TS_SD_00"
BOWTIEpe.filename <- "/data/virpand/pandemies/CAPTVRED/aln/TS_SD_00_pe.bowtie.sorted.mapped.bam"
BOWTIEsg.filename <- "/data/virpand/pandemies/CAPTVRED/aln/TS_SD_00_sg.bowtie.sorted.mapped.bam" 
contigs.gff.filename <- "/data/virpand/pandemies/CAPTVRED/coverage/TS_SD_00/TS_SD_00_blastn_on_viralcandidates.gff"
id.to.name <- "/data/virpand/pandemies/CAPTVRED/references/db/rvdb_nt/info_summary.tsv"
outplots.dir <-  "/data/virpand/pandemies/CAPTVRED/references/db/rvdb_nt/info_summary.tsv"
}
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
library(stringr);

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


## Read data info:
  #Read PE bam
  print ("Reading PE aligned reads file...")
pe.bam <- readGAlignments(BOWTIEpe.filename, # PE
                          use.names=TRUE)     # param=ScanBamParam(which=GRanges(BOWTIEpe.seq, IRanges(Xmin, Xmax))));
seq_lvs <- gsub("^[^|]+\\|[^|]+\\|([^|]+)\\|.*$", "\\1", seqlevels(pe.bam), perl=TRUE)
seqlevels(pe.bam) <- seq_lvs
seq_names <- gsub("^[^|]+\\|[^|]+\\|([^|]+)\\|.*$", "\\1", seqnames(pe.bam), perl=TRUE)
seqnames(pe.bam) <- seq_names
  
  #Reads SE bam
  print ("Reading SE aligned reads file...")
sg.bam <- readGAlignments(BOWTIEsg.filename, # SG
                        use.names=TRUE)    # param=ScanBamParam(which=GRanges(BOWTIEsg.seq, IRanges(Xmin, Xmax))));
seq_lvs2 <- gsub("^[^|]+\\|[^|]+\\|([^|]+)\\|.*$", "\\1", seqlevels(sg.bam), perl=TRUE)
seqlevels(sg.bam) <- seq_lvs2
seq_names2 <- gsub("^[^|]+\\|[^|]+\\|([^|]+)\\|.*$", "\\1", seqnames(sg.bam), perl=TRUE)
seqnames(sg.bam) <- seq_names2
  
  #Contigs from gff:
 print ("Reading contigs gff...")


ContigsFound=FALSE;
if (file.info(contigs.gff.filename)$size > 0) {
        print ("   Read gff info...")
        
         contigstbl <- read.table(contigs.gff.filename,sep="\t",header=FALSE)
         colnames(contigstbl) <- c("SeqID", "Source", "Feat", "Start", "End", "Score", "Strand", "Frame", "Info")
         contigstbl$name=gsub("ID=","",contigstbl$Info)
  #      contigs.txdb <- makeTxDbFromGFF(file=contigs.gff.filename,
  #                                  format="gff3");
        print ("   done!...")
        ContigsFound=TRUE;
}
#print ("   printing!...")
#print(contigs.txdb$user_seqlevels)

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
                        theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6),
                              #panel.border = element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    text = element_text(size = 10),
                                    plot.title = element_text(size = 15)) 
                return(theplot)
          }
    
print ("Reading coordinates table...")
#ref.coor<-read.table(refs.coord.tbl, row.names=1, col.names=c("ID", "START", "END"))
#dim(ref.coor) 

## Convert genomes ID to Human Readable name
 print ("Reading info table...")
genomes_rel <- read.table(id.to.name, header=TRUE, sep="\t", comment.char="", row.names=NULL)
colnames(genomes_rel) <- c("SeqId", "Name", "Tax", "Start", "End", "Description", "SpecId", "SpecTax", "Genus", "GenTax", "Family", "FamTax")
head(genomes_rel, 4)
for (txn in unique(genomes_rel$Tax)) { ## CHECK
     print(paste0("taxon is ",txn))
     txvec<-as.vector(genomes_rel[genomes_rel$Tax==txn,]$SeqId)
  #   print(length(txvec))
     plot_list <- list()
     for (rgn in txvec) {
        print("...")
        print(rgn)
        GFF.filename  <- paste0(GFF.refs.dir, "/", rgn ,".gff")
        rgn.v <- rgn
        rgn.id <- unlist(strsplit(rgn, '[.]'))[1];
        rgn.name <-  genomes_rel[genomes_rel$SeqId==rgn,]$Name
        
      if (!file.exists(GFF.filename)) {
            print(" ..")
            rgn <- rgn.id;
            GFF.filename  <- paste0(GFF.refs.dir, "/", rgn ,".gff")
        }
        print(paste0(">>>",GFF.filename))
        print("preparing txdb object")
        txdb <- makeTxDbFromGFF(file=GFF.filename,
                            format="gff3");
        # TxDb object
        ntx <- length(levels(as.factor(as.list(txdb)$transcripts$tx_id)));
        colors.genes <- ifelse(ntx > 1, hcl.colors(ntx, palette="Spectral"), "blue");
        # parameters
        Xmin  <- as.integer(genomes_rel[genomes_rel$SeqId==rgn.v,]$Start);
        Xmax  <- as.integer(genomes_rel[genomes_rel$SeqId==rgn.v,]$End); # WE NEED genome length
        Xmax <- round_any(Xmax, 500, f = ceiling)
      #  Xmin  <- as.integer(ref.coor[rgn.id, 1,]);
      #  Xmax  <- as.integer(ref.coor[rgn.id, 2,]); # WE NEED genome length
      #  Xmax <- round_any(Xmax, 500, f = ceiling)
        
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
   #     pesub.bam<-pe.bam[seqnames(pe.bam) == rgn.v,]     ## SI ESTA EMPTY FER UN PLOT EN BLANC
        pesub.bam<-pe.bam[seqnames(pe.bam) == rgn.v]
        pe.avgcov <- mean(coverage(pesub.bam));
        
        # plot
        print("Ploting PE_COV")
        if (length( pesub.bam) > 0) {
            print("there are reads!")
            pe.coverage <- autoplot(pesub.bam , geom = "line", stat = "coverage") +
                               theme_bw() + 
                               theme(text = element_text(size = 10)) + ylab("") +
                               scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0)) 
        }else{
            print("No reads found!")
            pe.coverage <- NA_plot( Xmax, 1, 'No alignments found')
        }
    #  Single ends coverage
        
        #subset bam
        sgsub.bam<-sg.bam[seqnames(sg.bam) == rgn.v]        ## SI ESTA EMPTY FER UN PLOT EN BLANC
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
        #print(contigs.txdb$user_seqlevels)
       ## FER SUBSET PER CONTIG
    subs.contigs.tbl<-contigstbl[contigstbl$SeqID==rgn.v,]
    if (nrow(subs.contigs.tbl) > 0) {
        
        #    if (any(contigs.txdb$user_seqlevels == rgn.v)) {
       #          if ( ContigsFound & any(contigs.txdb$user_seqlevels == rgn.v) ) {
                    print(" --- BBB --- ")
                    colors.genes <- "lightseagreen" #hcl.colors(1, palette="Spectral");
                    gr <- GRanges(seqnames=subs.contigs.tbl$SeqID, 
                                  ranges=IRanges(start=subs.contigs.tbl$Start, end=subs.contigs.tbl$End),
                                  gene_id=subs.contigs.tbl$name,
                                  score=subs.contigs.tbl$Score
                                 )
                    
                    contigs <- autoplot(gr, geom="rect", fill=colors.genes) +
                                theme_bw() +
                                theme(panel.border = element_blank(),
                                    axis.text.y=element_blank(),
                                    text = element_text(size = 10),
                                    plot.title = element_text(size = 15)) +
                                scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0));
                    ## CARREGAR EN GRANGES
     #               gr <- GRanges(seqnames=GFF.seq, # Rle(rep(GFF.seq, length(blst2$SeqId))),  # blst2$SeqId),
     #                           IRanges(start=blst2$Start,
     #                           end=blst2$End,
     #                           names=substr(blst2$Group,4,length(blst2$Group))),
     #                           Rle(strand(blst2$Strand)),
     #                           score=blst2$Score, seqinfo=GNM)

                       
                   ## AUTOPLOT
  #                 p3 <- autoplot(gr, geom="rect", aes(fill=score)) + labs(fill = "%Identity") +
  #    scale_fill_gradient(low=colors.set[2], high=colors.set[4],
  #                        breaks=seq(90,100,2.5), limits=c(90,100)) +
  #    theme_bw() + theme(axis.text.y=element_blank(), text = element_text(size = 10),
  #                       legend.justification = c(0,1), legend.position = c(0, 1),
  #                       legend.title = element_text(size = textSize/1.25),
  #                       legend.text  = element_text(size = textSize/1.5),
  #                       legend.key.size = unit(spaceLegend, "lines")) +
  #    scale_y_discrete(breaks=NULL) + ylab("") +
  #    scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0));

  #                  contigs <- autoplot(gr, geom="rect")  #which=GRanges(rgn.v, IRanges(Xmin, Xmax)),
  #                              names.expr = "transcript", fill=colors.genes) +
  #                              theme_bw() +
  #                              theme(panel.border = element_blank(),
  #                                  axis.text.y=element_blank(),
  #                                  text = element_text(size = 10),
  #                                  plot.title = element_text(size = 15)) +
  #                              scale_x_continuous(limits = c(Xmin, Xmax), expand = c(0, 0));
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
      print("RSaving plot...")
      ggsave(PNG.filename, plot=wholeplot,
              width = 25, height = 15, units = "cm", dpi = 600);
      print("DONE!")
    };
};