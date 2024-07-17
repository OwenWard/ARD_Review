## set up log

#mccarty.log<-"c:/docume~1/tianzh~1/mydocu~1/Tian_research/Project_Michael/log/mccarty12.log"
#setwd("C:/Documents and Settings/Tyler/Desktop/ARD/mixing paper/Project_Name2")


#sink(mccarty.log)
print(date())
ct.trunc<-30

## the real thing

library(foreign)
mccarty <- read.dta("all.dta.dta")
population <- read.table ("population.txt", header=F, skip=1, sep=";")
popusize <- read.table ("popusize.csv", header=F, sep=",")

#mccarty <- read.dta(file.choose())
#population <- read.table (file.choose(), header=F, skip=1, sep=";")
#popusize <- read.table (file.choose(), header=F, sep=",")

## read in data
## can not find information about clergy or not

## check data consistence
### using only set 1 and set 2
mccarty.all<-mccarty
mccarty.12<-mccarty[mccarty[,'set']<3,]

mccarty<-mccarty.12

state.1<-table(c(mccarty[mccarty[,'set']==1, 68], 1:51), exclude=c(NA, 59))-rep(1, 51)
state.2<-table(c(mccarty[mccarty[,'set']==2, 68], 1:51), exclude=c(NA, 59))-rep(1, 51)
state.3<-table(c(mccarty[mccarty[,'set']==3, 68], 1:51), exclude=c(NA, 59))-rep(1, 51)

print(cor(cbind(state.1, state.2, state.3)))

age.1<-mccarty[mccarty[,'set']==1, 3]
age.2<-mccarty[mccarty[,'set']==2, 3]
#age.3<-mccarty[mccarty[,'set']==3, 3]

print(ks.test(age.1, age.2))
#print(ks.test(age.1, age.3))
#print(ks.test(age.2, age.3))

print(table(mccarty[, c(75, 60)]))
print(table(mccarty[, c(75, 61)]))

## get gender information

gender<-mccarty[,2]
age<-mccarty[,3]

## get rid of other variables

mccarty<-mccarty.12[, c(1, 5, 7:37)]

labels <- as.vector(population[,2])
category.short<-c('Michael', 'Christina', 'Christopher', 'Jacqueline', 'James', 'Jennifer', 'Anthony', 'Kimberly', 'Robert', 'Stephanie', 'David', 'Nicole', 'Amer. Indians', 'New Birth', 'Adoption', 'Widow(er)', 'Dialysis', 'Postal Worker', 'Pilot', 'Jaycees', 'Diabetic', 'Business', 'Twin', 'Gun Dealer', 'HIV positive', 'AIDS', 'Homeless', 'Rape', 'Prison', 'Homicide', 'Suicide', 'Auto Accident')
 
tempdata <- array(NA, c(nrow(mccarty), ncol(mccarty)-1))

for (j in 1:ncol(tempdata))
  tempdata[,j]<-as.numeric(as.vector(mccarty[,j+1]))

tempdata <- ifelse ((tempdata==98)|(tempdata==7.5), NA, tempdata)

## Start 

## making histograms of each question
par (mfrow=c(3,3))
for (j in 1:ncol(tempdata)){
  hist(tempdata[,j], main=labels[j], breaks=seq(-.5,max(tempdata[,j],na.rm=T)+.5),
        xlim=c(0, ct.trunc))
}

total <- apply (tempdata, 1, sum, na.rm=T)
variability<-apply(tempdata, 1, var, na.rm=T)
variability<-ifelse(is.na(variability), 0, variability)

tempdata[(variability==0),] <- rep(NA, ncol(tempdata))

# only those know people and know different number of people 
# from different subpopulations will be included since those 
# who answered a same number to all questions could having 
# been 'lying' or 'mis-recall'.

# calculate some data summaries

y <- tempdata
n <- nrow(y)
J <- ncol(y)

i.mean <- apply (y, 1, mean, na.rm=T)
j.mean <- apply (y, 2, mean, na.rm=T)
ij.mean <- mean(i.mean, na.rm=T)

# simple imputation

ok <- is.na(y)
y.imp <- ifelse (ok, outer (i.mean,j.mean)/ij.mean, y)

# simple estimates of averages by person and question

a.hat <- rep (NA, n)
for (i in 1:n){
  a.hat[i] <- weighted.mean (y.imp[i,], 1/sqrt(j.mean), na.rm=T)
}

b.hat <- rep (NA, J)
for (j in 1:J){
  b.hat[j] <- weighted.mean (y.imp[,j], 1/sqrt(i.mean), na.rm=T)
}

y.hat <- outer(a.hat,b.hat)
y.hat <- y.hat*mean(y[!ok])/mean(y.hat[!ok])

# calculate oversdispersion for each question

overdisp <- rep (NA, J)
for (j in 1:J){
  overdisp[j] <- mean ((y[,j]-y.hat[,j])^2/y.hat[,j], na.rm=T)
}

#

y <- ifelse (y>ct.trunc, ct.trunc, y)

i.mean <- apply (y, 1, mean, na.rm=T)
j.mean <- apply (y, 2, mean, na.rm=T)
ij.mean <- mean(i.mean, na.rm=T)

# simple imputation

ok <- is.na(y)
y.imp <- ifelse (ok, outer (i.mean,j.mean)/ij.mean, y)

# simple estimates of averages by person and question
a.hat <- rep (NA, n)
for (i in 1:n){
  a.hat[i] <- weighted.mean (y.imp[i,], 1/sqrt(j.mean), na.rm=T)
}
b.hat <- rep (NA, J)
for (j in 1:J){
  b.hat[j] <- weighted.mean (y.imp[,j], 1/sqrt(i.mean), na.rm=T)
}

y.hat <- outer(a.hat,b.hat)
y.hat <- y.hat*mean(y[!ok])/mean(y.hat[!ok])

# calculate oversdispersion for each question

overdisp2 <- rep (NA, J)
for (j in 1:J){
  overdisp2[j] <- mean ((y[,j]-y.hat[,j])^2/y.hat[,j], na.rm=T)
}

cat("\n category name \t overdisp \t overdisp (truncated)\n")

cbind (category.short, round(overdisp,1), round(overdisp2,1))

par (mfrow=c(1,1))
plot (j.mean, overdisp2, ylim=c(1,max(overdisp2)), type="n")
names <- rep(c(.7,.8),c(12,20))
text (j.mean, overdisp2, category.short, cex=names)

#sink()
