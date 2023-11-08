#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
gpsc.list<-c(10,2,22,26,5,8)
for (i in gpsc.list){
neutral.laneIDS<-read.csv(paste0("./GPSCs/GPSC",i,"/out_gpsc",i,"/gpsc",i,".",args[1],".laneIDS"),header=FALSE)
GPS.demes<-readRDS("./Metadata/GPS.demes.RData")
demes<-c("SOUTH_AFRICA","KENYA","MALAWI","THE_GAMBIA")

####84 neutral genes
##subset the lanes which include this gene
neutral.demes<-GPS.demes[which(GPS.demes$Lane%in%neutral.laneIDS$V1),]
###reorder to match the alignment
neutral.demes<-neutral.demes[match(neutral.laneIDS$V1,neutral.demes$Lane),]
###extract the indexes for each population
mal_neutral<-which(neutral.demes$Country=="MALAWI")
sa_neutral<-which(neutral.demes$Country=="SOUTH_AFRICA")
ken_neutral<-which(neutral.demes$Country=="KENYA")
gam_neutral<-which(neutral.demes$Country=="THE_GAMBIA")
write.table(sa_neutral,file=paste0("./GPSCs/GPSC",i,"/index/sa_",args[1],"_index.csv"),row.names = FALSE,col.names = FALSE)
write.table(mal_neutral,file=paste0("./GPSCs/GPSC",i,"/index/mal_",args[1],"_index.csv"),row.names = FALSE,col.names = FALSE)
write.table(ken_neutral,file=paste0("./GPSCs/GPSC",i,"/index/ken_",args[1],"_index.csv"),row.names = FALSE,col.names = FALSE)
write.table(gam_neutral,file=paste0("./GPSCs/GPSC",i,"/index/gam_",args[1],"_index.csv"),row.names = FALSE,col.names = FALSE)
}
