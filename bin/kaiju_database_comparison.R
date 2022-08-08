args = commandArgs(trailingOnly=TRUE)

library(UpSetR)
#library(ComplexHeatmap)
#library(remotes)
#library(ggplot2)
#euk<-read.table("taxonomy/kaiju/R01_C10_G02_all.nr_euk.ids", sep=",")
#refs<-read.table("taxonomy/kaiju/R01_C10_G02_all.refseq.ids", sep=",")
#rvdbP<-read.table("taxonomy/kaiju/R01_C10_G02_all.rvdb.ids", sep=",")
#vir<-read.table("taxonomy/kaiju/R01_C10_G02_all.viruses.ids", sep=",")

euk<-read.table(args[1], sep=",")
refs<-read.table(args[2], sep=",")
rvdbP<-read.table(args[3], sep=",")
vir<-read.table(args[4], sep=",")

lt_byid <- list( nr_euk=as.vector(euk$V1),  refseqs=as.vector(refs$V1),  rvdb=as.vector(rvdbP$V1), viruses=as.vector(vir$V1) )
#m1 = make_comb_mat(lt_byid)

png(filename=paste0(args[5], "_reads_found_by_database.png"), width = 1200, height = 700)
upset(fromList(lt_byid), order.by = "freq",
        point.size = 4.5, line.size = 1, 
        text.scale=c(2, 1.5, 2, 1.5, 2, 2.5), 
        mainbar.y.label = "Number of reads classified", sets.x.label = "Hits per data base"
     )
dev.off()


lt_bysp <- list( nr_euk=unique(as.vector(euk$V3)),  refseqs=unique(as.vector(refs$V3)),  rvdb=unique(as.vector(rvdbP$V3)), viruses=unique(as.vector(vir$V3)) )
#m2 = make_comb_mat(lt_bysp)
png(filename=paste0(args[5], "_species_found_by_database.png"), , width = 1200, height = 700)
upset(fromList(lt_bysp), order.by = "freq",
 point.size = 4.5, line.size = 1, 
        text.scale=c(2, 1.5, 2, 1.5, 2, 2.5), 
        mainbar.y.label = "Number of species", sets.x.label = "Hits per data base"
     )
dev.off()

