library(nlme)
library(mgcv)
library(metafor)
library(R.matlab)

InputPath <- '/Documents/R/Meta'
Data <- read.csv(file="Participant.csv", header=TRUE)
Data$GroupName=factor(Data$Group, levels = c(0,1),  labels = c( "HC","ASD"),ordered=TRUE)
contrasts(Data$GroupName)<-"contr.treatment"
contrasts(Data$GroupName)

MV_Zscore <-readMat(paste0(InputPath,'/MV_Zscore.mat'));
MV_Zscore <- data.frame(matrix(unlist(MV_Zscore),nrow = 939))
MV_Zscore_WithSiteIndex<- MV_Zscore 
MV_Zscore_WithSiteIndex$SiteID <- Data$SiteID

D_value<-matrix(NA,nrow=512,ncol=18)
SE_value<-matrix(NA,nrow=512,ncol=18)
for (Site in 1:18) {
  Data_Site <- subset(Data,SiteID == Site)
  MV_Zscore_Site <- subset(MV_Zscore_WithSiteIndex,SiteID == Site)
for (region in 1:512) {
  Data_Site$RegionMeasure<-MV_Zscore_Site[,region]
  mod_gam <- gam(RegionMeasure~GroupName+s(Age,k=4)+s(Age,by=GroupName,k=4)+MeanFD,data=Data_Site,method = "REML")
  ASD <- length(which(Data_Site$Group == 1))
  HC <- length(which(Data_Site$Group == 0))
  D <-summary(mod_gam)$p.t[2]*sqrt(1/ASD+1/HC)
  D_value[region,Site] <- D
  SE_value[region,Site]<-summary(mod_gam)[["se"]][["GroupNameASD"]]
  
 }
}
write.table( D_value, file = paste("D_value.txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
write.table( SE_value, file = paste("SE_value.txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )

MetaResults<-matrix(NA,nrow=512,ncol=3)
for (region in 1:512) {
  D<-D_value[region,]
  SE<-SE_value[region,]
  data_SiteN<-data.frame(D,SE)
res <- rma.uni(yi = D, sei = SE, method = "REML", data = data_SiteN)
MetaResults[region,]=c(res[["beta"]],res[["I2"]],res[["pval"]])
}

write.table(MetaResults, file = paste("MetaResults.txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )

