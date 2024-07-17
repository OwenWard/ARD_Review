library(foreign)
library(car)
recode.approx <- function (a){
  recode (as.vector(a), "0=0; 1=1; 2=3.5; 3=8; 4=15")
}
rm(list=ls(all=TRUE))

setwd("C:/Documents and Settings/Tyler/Desktop/GSS Project/FullLikelihood_NoImp/acq-likelihood/Workspaces/")





##################with hispanics##################with hispanics##################with hispanics
load("acq_recall.RData")
#the variable 'y' is the data used in the simulation file; acq_run.R
#ind.n is the num of rows of y
#using the posterior mean
est.deg<-round(cbind(exp(monitor.sims.all[(1):(ind.n),c(1,5)])),3)

#now attach the id's
id.all<-read.table("C:/Documents and Settings/Tyler/Desktop/GSS Project/FullLikelihood_NoImp/id_numeric.txt",header=T)
rm(data)
rm(tempdata)
#remove all that did not respond to any questions in the section
data<-approx[,c(1,4:26)]
tempdata <- array(NA, c(nrow(data), ncol(data)))
for (j in 2:ncol(tempdata))
  tempdata[,j]<-as.numeric(as.vector(data[,j]))
tempdata<-cbind(id.all,tempdata[,-1])
#no values of either 7.5 or 98
#tempdata[,-1] <- ifelse ((tempdata[,-1]==98)|(tempdata[,-1]==7.5), NA, tempdata)
tempdata<-tempdata[raw$numknown==1,]
mis.data<-rowSums(is.na(tempdata[,-1]))
tempdata<-tempdata[mis.data<ncol(tempdata[,-1]),]
id.vec<-tempdata[,1]

#put all together
est.deg.out<-cbind(id.vec,est.deg)
colnames(est.deg.out)<-c("id","Posterior Mean","Posterior Median")
write.csv(est.deg.out,file="../Results/deg_acq_hisp.csv")


rm(list=ls(all=TRUE))


##################without hispanics##################without hispanics##################without hispanics
load("acq_recall_noHis.RData")
#the variable 'y' is the data used in the simulation file; acq_run.R
#ind.n is the num of rows of y
#using the posterior mean
est.deg<-round(cbind(exp(monitor.sims.all[(1):(ind.n),c(1,5)])),3)

#now attach the id's
id.all<-read.table("C:/Documents and Settings/Tyler/Desktop/GSS Project/FullLikelihood_NoImp/id_numeric.txt",header=T)
rm(data)
rm(tempdata)
#remove all that did not respond to any questions in the section
data<-approx[,c(1,4:11,14:26)]
tempdata <- array(NA, c(nrow(data), ncol(data)))
for (j in 2:ncol(tempdata))
  tempdata[,j]<-as.numeric(as.vector(data[,j]))
tempdata<-cbind(id.all,tempdata[,-1])
#no values of either 7.5 or 98
#tempdata[,-1] <- ifelse ((tempdata[,-1]==98)|(tempdata[,-1]==7.5), NA, tempdata)
tempdata<-tempdata[raw$numknown==1,]
mis.data<-rowSums(is.na(tempdata[,-1]))
tempdata<-tempdata[mis.data<ncol(tempdata[,-1]),]
id.vec<-tempdata[,1]

#put all together
est.deg.out<-cbind(id.vec,est.deg)
colnames(est.deg.out)<-c("id","Posterior Mean","Posterior Median")
write.csv(est.deg.out,file="../Results/deg_acq_noHisp.csv")
