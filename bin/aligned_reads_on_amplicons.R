library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
#outfl<-"amplicons_alignment/figures/RUN01_WGSvsCAP"
outfl <- args[3] 

#DT <- read.table("RUN01_bowtie_info.tbl", header=FALSE, sep="\t")
DT <- read.table(args[1], header=FALSE, sep="\t")
colnames(DT)=c("readID", "REFgenome", "QUAL", "seqLEN", "sampID")

genomes_rel<- read.table(args[2], header=FALSE, sep="\t")
#genomes_rel <- read.table("Viral_candidates_zoonosis_refseqs_idrelation.tsv", header=TRUE, sep="\t")
colnames(genomes_rel)<-c("NC", "HRsimp", "HRcomp" )
NCSimp <- unlist(str_split(genomes_rel$NC, "[.]", simplify=T))[,1]
genomes_rel$NCsimp <- NCSimp

HRids<-genomes_rel$HRcomp
names(HRids)<-genomes_rel$NCsimp
HRids["UNC"]="UNC"

HRids2<-genomes_rel$HRsimp
names(HRids2)<-genomes_rel$NCsimp
HRids2["UNC"]="UNC"


DT <- DT %>% 
  mutate(
    SEQAP = if_else(
      condition = str_detect(sampID, "_M\\d\\d_" ), 
      true      = "WGS", 
      false     = if_else(
	      condition = str_detect(sampID, "_C\\d\\d_" ), 
	      true      = "CAP", 
	      false     = "UNC"
	          )
             )
        ) %>% 
      mutate(
        SORI = if_else(
          condition = str_detect(sampID, "_G\\d\\d_" ), 
          true      = "BATGNO", 
          false     = if_else(
              condition = str_detect(sampID, "_S\\d\\d_" ), 
              true      = "SEWAGE", 
              false     = "UNCLAS"
                         )
                       )
           ) %>% 
          mutate(
            ENDS = if_else(
              condition = str_detect(sampID, "_pe" ), 
              true      = "PEND", 
              false     = if_else(
                  condition = str_detect(sampID, "_sg" ), 
                  true      = "SEND", 
                  false     = "UNCLAS"
                             )
                           )
               )
               
               
SP <- unlist(str_split(DT$sampID, "_", simplify=T))[,3]
DT$SP<-SP
sprefg<-  unlist(str_split(DT$REFgenome, '[.]', simplify=T))[,1]
Name <- HRids[sprefg]
NameSimp <- HRids2[sprefg]
DT$REFgenSimp <- sprefg
DT$REFName <- Name
DT$REFNameSimp <- NameSimp
levels(DT$REFName) <- unique(unname(HRids))
levels(DT$REFNameSimp) <- unique(unname(HRids2))
DT$ENDS <- factor( DT$ENDS, levels=c("SEND", "PEND") )



#group.colors <- c(CAP = "#0453a2" , WGS = "#ffb84e" , UNC ="#CCCCCC")
#group.colors <- c(PEND ="#de7966", SEND ="#f4a292" , UNC ="#CCCCCC") 
group.colors <- c(PEND ="#de7966", SEND ="#ffb7aa" , UNC ="#CCCCCC") 



#### Aligned reads ####
bysamp_all_v <- ggplot(DT, aes(x=factor(REFNameSimp), fill=ENDS)) + 
	        geom_bar() +
            scale_x_discrete(limits=factor(levels(DT$REFNameSimp)), drop=FALSE) +
            scale_y_continuous(labels = function(x) format( log10(x), scientific = FALSE), trans='log10') +
	        coord_flip() +
	        scale_fill_manual(values=group.colors) +
            theme_bw() + 
            facet_wrap(~SP, nrow=1) +
            ggtitle( "Mapped sequences into reference genomes by sample" ) +
            ylab("Log10 Counts of mapped reads") + 
            xlab("Reference genomes" ) +
            theme( 
               axis.title = element_text( size = 15, color="#222222", hjust=0.5),
               plot.title = element_text( size = 18, face = "bold", hjust=0.5),
               strip.text = element_text(size = 9),
               legend.title = element_text(size=10), #change legend title font size
               legend.text = element_text(size=9)
            )
          
          
ggsave(plot=bysamp_all_v, filename=paste0(outfl,"_maped_reads_by_sample_allvirus.ok.png"),width=25, height=10, dpi=300) 

quit()

####   other explorative analysis - Discarded for the pipeline   ####

####  Quality distribution  ####

 q1 <- ggplot(DT, aes(x=QUAL)) +
   geom_histogram( bins=27, fill="#994d69") +
   #geom_density(alpha=.3, color="#d79aaa", fill="#d79aaa") +
   geom_bar(color="#994d69", fill="#994d69") +
   facet_wrap(~SP, ncol=3, scales="free") +
   theme_bw()
 ggsave(plot=q1, filename=paste0(outfl,"quality_distribution.png"),width=25, height=10, dpi=300) 
 
###

## By specie ## discarded
byspec_all_v <- ggplot(DT, aes(x=factor(SP), fill=ENDS)) + 
	        geom_bar() +
            #scale_x_discrete(limits=factor(levels(DT$SP)), drop=TRUE) +
	        #coord_flip() +
	        scale_fill_manual(values=group.colors) +
            theme_bw() + 
            facet_wrap(~REFNameSimp) +
            scale_y_log10() +
            ggtitle( "Mapped sequences into reference genomes by sample" ) +
            ylab("Log10 Counts of mapped reads") + 
            xlab("Reference genomes" ) +
            theme( 
               axis.title = element_text( size = 15, color="#222222", hjust=0.5),
               plot.title = element_text( size = 18, face = "bold", hjust=0.5),
               strip.text = element_text(size = 9),
               legend.title = element_text(size=10), #change legend title font size
               legend.text = element_text(size=9)
            )
          
          
ggsave(plot=bysamp_all_v, filename=paste0(outfl,"_aligned_reads_by_sample_allvirus.ok.png"),width=25, height=10, dpi=300) 
###
