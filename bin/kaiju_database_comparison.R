args = commandArgs(trailingOnly=TRUE)


#This script generates two static UpSet plots to visualize
#the interaction of kaiju outputs in 5 different databases
#in terms of species and number of reads respectively


library(gplots)
library(UpSetR)
## DEFINE INPUT FILES ##
euk<-read.table(args[1], sep=",") 
refs<-read.table(args[2], sep=",")
rvdbP<-read.table(args[3], sep=",")
vir<-read.table(args[4], sep=",")
cap<-read.table(args[5], sep=",")
# files args 1 to 5 are csv consisting con 3 columns:
        # column 1: read id
        # column 2: kingdom
        # column 3: specie
 #this csv files are obtained from kaiju output. (nee nextflow protocol)

## DEFINE FUNCTIONS ##

overlapGroups <- function (listInput, sort = TRUE) {
        #font: https://github.com/hms-dbmi/UpSetR/issues/85
  listInputmat    <- fromList(listInput) == 1
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
    myelements
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}

## INTERSECTION BY IDENTIFIER ##

lt_byid <- list( nr_euk=as.vector(euk$V1),  refseqs=as.vector(refs$V1),  rvdb=as.vector(rvdbP$V1), viruses=as.vector(vir$V1), captureSeqs=as.vector(cap$V1) )
  
  # Plot
print("PLOTTTING")
png(filename=paste0(args[6], "_reads_found_by_database.png"), width = 1200, height = 700)
upset(fromList(lt_byid), order.by = "freq",
        point.size = 4.5, line.size = 1, number.angles = 20,
        text.scale=c(1.5, 1.5, 2, 1.5, 2, 1.5), 
        mainbar.y.label = "Number of reads classified", sets.x.label = "Hits per data base"
     )
dev.off()
print("DONE")

 # List of names
print ("LISTING")
ItemsList <- venn(lt_byid, show.plot=FALSE)
int<-(attributes(ItemsList)$intersections)
#outfl=paste0(args[6], "_ids_sets.txt")
#cat("", file=outfl, append=FALSE, sep="\n");
for (i in 1:length(int)){
         print(i);
         outfl=paste0(args[6], "_", names(int[i]), "_", "ids.txt")
         #cat(paste0("#", names(int[i]), "\t", length(int[i][[1]])), file=outfl, append=TRUE, sep="\n");
         cat(paste0("#", names(int[i]), "\t", length(int[i][[1]])), file=outfl, append=FALSE, sep="\n");
         for (j in int[i]) {
                 cat(j,file=outfl, append=TRUE, sep="\n")
                 }
        }
print("DONE")


#q();

#li_id <- overlapGroups(lt_byid)
#head(li_id)
#save(li_id, file=paste0(args[6], "_ids_sets.txt"))

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
dev.off()

 #List of names
ItemsList2 <- venn(lt_bysp, show.plot=FALSE)
int2<-(attributes(ItemsList2)$intersections)

for (i in 1:length(int)){
         print(i);
         outfl=paste0(args[6], "_", names(int[i]), "_", "spec.txt")
         #cat(paste0("#", names(int[i]), "\t", length(int[i][[1]])), file=outfl, append=TRUE, sep="\n");
         cat(paste0("#", names(int[i]), "\t", length(int[i][[1]])), file=outfl, append=FALSE, sep="\n");
         for (j in int[i]) {
                 cat(j,file=outfl, append=TRUE, sep="\n")
                 }
        }
print("DONE")


