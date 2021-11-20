library(nlme)
library(mgcv)
library(R.matlab)

InputPath <- '/Documents/R/Mega'
Data <- read.csv(file="Participant.csv", header=TRUE)
Data$GroupName=factor(Data$Group, levels = c(0,1),  labels = c( "HC","ASD"),ordered=TRUE)
contrasts(Data$GroupName)<-"contr.treatment"
contrasts(Data$GroupName)

MV_Harmonized_Zscore <-readMat(paste0(InputPath,'/MV_Harmonized_Zscore.mat'))
MV_Harmonized_Zscore <-data.frame(matrix(unlist(MV_Harmonized_Zscore),nrow = 939))


p_value<-matrix(NA,nrow=512,ncol=3)
t_F_value<-matrix(NA,nrow=512,ncol=3)
D_value<-matrix(NA,nrow=512,ncol=1)
ASD<-440
HC<-499
for (region in 1:512) {
  Data$RegionMeasure<-MV_Harmonized_Zscore[,region]
  mod_gam <- gam(RegionMeasure~GroupName+s(Age,k=4)+s(Age,by=GroupName,k=4)+MeanFD,data=Data,method = "REML")
  summary_model<-summary(mod_gam)
  p_value[region,]<-c(summary_model$p.pv[2],summary_model$s.pv[1],summary_model$s.pv[2])
  t_F_value[region,]<-c(summary_model$p.t[2],summary_model[["s.table"]][5],summary_model[["s.table"]][6])
  D_value[region,1]<-summary(mod_gam)$p.t[2]*sqrt(1/ASD+1/HC)
  
}
write.table( p_value, file = paste( "MV_Harmonized_Zscore_p_value.txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
write.table( t_F_value, file = paste( "MV_Harmonized_Zscore_t_F_value.txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )
write.table(D_value, file = paste( "MV_Harmonized_Zscore_D_value.txt", sep="" ), row.names = FALSE, col.names = FALSE, quote = FALSE )




