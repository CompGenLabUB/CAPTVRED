args = commandArgs(trailingOnly=TRUE)


#This script generates two static UpSet plots to visualize
#the interaction of kaiju outputs in 5 different databases
#in terms of species and number of reads respectively


library(gplots)
library(UpSetR)
library(stringr)
library(grid)
## DEFINE INPUT FILES ##
print(file.exists(args[1]))
if(file.exists(args[1]) & file.size(args[1])>0){
  euk<-read.table(args[1], sep=",") 
} else {
  euk<- data.frame(matrix("NAN", nrow=1, ncol=3)); 
  colnames(euk)=c("V1", "V2", "V3")
}

if(file.exists(args[2]) & file.size(args[2])>0){
  refs<-read.table(args[2], sep=",")
} else {
  refs<- data.frame(matrix("NAN", nrow=1, ncol=3)); 
  colnames(refs)=c("V1", "V2", "V3")
}
if(file.exists(args[3]) & file.size(args[3])>0){
  rvdbP<-read.table(args[3], sep=",")
} else {
  rvdbP<-data.frame(matrix("NAN", nrow=1, ncol=3)); 
  colnames(rvdbP)=c("V1", "V2", "V3")
  }
  
if(file.exists(args[4]) && file.size(args[4])>0){
  vir<-read.table(args[4], sep=",")
} else {
  vir<-data.frame(matrix("NAN", nrow=1, ncol=3)); 
  colnames(vir)=c("V1", "V2", "V3")
  }


if(file.exists(args[5]) & file.size(args[5])>0){
  cap<-read.table(args[5], sep=",")
} else {
  cap<-data.frame(matrix("NAN", nrow=1, ncol=3)); 
  colnames(cap)=c("V1", "V2", "V3")
}

sampid=rev(str_split(args[6], "/")[[1]])[1] 
print(sampid)
# files args 1 to 5 are csv consisting con 3 columns:
        # column 1: read id
        # column 2: kingdom
        # column 3: specie
 #this csv files are obtained from kaiju output. (nee nextflow protocol)

#####   1. UpSetR plot  ######

## INTERSECTION BY IDENTIFIER ##

lt_byid <- list( nr_euk=as.vector(euk$V1),  refseqs=as.vector(refs$V1),  rvdb=as.vector(rvdbP$V1), viruses=as.vector(vir$V1), captureSeqs=as.vector(cap$V1) )
  
  # Plot  -> ## falta descartar si id=NA
print("PLOTTTING")

png(filename=paste0(args[6], "_reads_found_by_database.png"), width = 1200, height = 700)
upset(fromList(lt_byid), order.by = "freq",
        point.size = 4.5, line.size = 1, number.angles = 20,
        text.scale=c(1.5, 1.5, 2, 1.5, 2, 1.5),
        mainbar.y.label = "Number of reads classified", sets.x.label = "Hits per data base"
     ) 
grid.text(sampid, x = 0.65, y = 0.90,
              gp = gpar(fontsize = 20))
dev.off()
print("DONE")

 # List of names
print ("LISTING")
ItemsList <- venn(lt_byid, show.plot=FALSE)
int<-(attributes(ItemsList)$intersections)
#outfl=paste0(args[6], "_ids_sets.txt")
#cat("", file=outfl, append=FALSE, sep="\n");
for (i in 1:length(int)){
         outfl=paste0(args[6], "_", names(int[i]), "_", "ids.txt")
         #cat(paste0("#", names(int[i]), "\t", length(int[i][[1]])), file=outfl, append=TRUE, sep="\n");
         cat(paste0("#", names(int[i]), "\t", length(int[i][[1]])), file=outfl, append=FALSE, sep="\n");
         for (j in int[i]) {
                 cat(j,file=outfl, append=TRUE, sep="\n")
                 }
        }
print("DONE")


## INTERSECTION BY SAMPLE ## 
 # Plot
lt_bysp <- list( nr_euk=unique(as.vector(euk$V3)),  
                refseqs=unique(as.vector(refs$V3)),  
                rvdb=unique(as.vector(rvdbP$V3)), 
                viruses=unique(as.vector(vir$V3)), 
                captureSeqs=unique(as.vector(cap$V3))
                )

head (lt_bysp$nr_euk )
head( lt_bysp$captureSeqs )

png(filename=paste0(args[6], "_species_found_by_database.png"), , width = 1200, height = 700)
upset(fromList(lt_bysp), order.by = "freq",
 point.size = 4.5, line.size = 1, 
        text.scale=c(2, 1.5, 2, 1.5, 2, 2.5), 
        mainbar.y.label = "Number of species", sets.x.label = "Hits per data base"
     )
grid.text(sampid, x = 0.65, y = 0.90,
              gp = gpar(fontsize = 20))
dev.off()

 #List of names
ItemsList2 <- venn(lt_bysp, show.plot=FALSE)
int2<-(attributes(ItemsList2)$intersections)

for (i in 1:length(int)){
         outfl=paste0(args[6], "_", names(int[i]), "_", "spec.txt")
         #cat(paste0("#", names(int[i]), "\t", length(int[i][[1]])), file=outfl, append=TRUE, sep="\n");
         cat(paste0("#", names(int[i]), "\t", length(int[i][[1]])), file=outfl, append=FALSE, sep="\n");
         for (j in int[i]) {
                 cat(j,file=outfl, append=TRUE, sep="\n")
                 }
        }


#####   2. HEATMAP plot  ######

library(tidyverse)

colnames(euk)<-c("ID", "KG", "SP")
euk$DB<-rep("nr_euk", dim(euk)[1])

colnames(refs)<-c("ID", "KG", "SP")
refs$DB<-rep("refseqs", dim(refs)[1])
head(refs)

colnames(rvdbP)<-c("ID", "KG", "SP")
rvdbP$DB<-rep("rvdb", dim(rvdbP)[1])

colnames(vir)<-c("ID", "KG", "SP")
vir$DB<-rep("viruses", dim(vir)[1])

colnames(cap)<-c("ID", "KG", "SP")
cap$DB<-rep("CapSeq", dim(cap)[1])

head(vir)

DT <- rbind(euk, refs, rvdbP, vir, cap)  %>% 
        count(SP, DB)%>% 
        complete(SP, DB, fill=list(0)) %>%
        replace(is.na(.), 0)
print(dim(DT))
head(DT)
DT <- DT[DT$SP != "NAN",]
print(head(DT))

png(filename=paste0(args[6], "_counts_heatmap.png"), width = 1200, height = 700)
ggplot(DT, aes(x=DB, y=SP, fill=log(n))) + 
  geom_tile() +
  theme_bw() +
  scale_fill_gradient(low = "mintcream", high = "mediumseagreen", na.value="white") +
  scale_x_discrete(position="top") +
  #geom_text(aes(label=ifelse(as.numeric(n)>0,n,""))) +
  ggtitle(sampid) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(face="bold", size=13), 
        axis.text.y = element_blank(), 
        legend.position="bottom",
        plot.title = element_text(hjust=.5, face="bold", size=20)
   )
dev.off()


#fer un heatmap a partdir de DT
cnts<-DT %>%  spread(DB, n) %>% replace(is.na(.), 0)  
write.csv(cnts, paste0(args[6], "species_by_database.csv"), row.names=FALSE)
#Save CNTS as table

