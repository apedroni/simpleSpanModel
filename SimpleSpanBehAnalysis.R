# Simple behavioral analysis of the test-retest pilot data 2017

rm(list=ls())
setwd("~/Dropbox/AA_Neurometric/ANALYSES/SimplespanAnalyses/Modeling/")
library(ggplot2)
require(R2jags) 
library(tidyr)
library(plyr)
library(reshape2)
#library(HDInterval)
## read and prepare data

data = as.matrix(read.table("SimpleSpanResponses.csv",sep=',',header = TRUE))

## all subjects
par(mfrow=c(2,4))
T1data = data[,1:6]==3
barplot(apply(T1data,2,mean), main ="Correct", ylim = c(0,1))

T2data = data[,7:12]==3
barplot(apply(T2data,2,mean), main ="Correct", ylim = c(0,1))

T1data = data[,1:6]==2
barplot(apply(T1data,2,mean), main ="WrongPos", ylim = c(0,1))

T2data = data[,7:12]==2
barplot(apply(T2data,2,mean), main ="WrongPos", ylim = c(0,1))

T1data = data[,1:6]==1
barplot(apply(T1data,2,mean), main ="Intrusion", ylim = c(0,1))

T2data = data[,7:12]==1
barplot(apply(T2data,2,mean), main ="Intrusion", ylim = c(0,1))

T1data = data[,1:6]==0
barplot(apply(T1data,2,mean), main ="Wrong", ylim = c(0,1))

T2data = data[,7:12]==0
barplot(apply(T2data,2,mean), main ="Wrong", ylim = c(0,1))


data = cbind(data, rep(c(1:23), each = 16))

## bad performers by median split
corr=array()
dd=array(NA,dim=c(23,6))
for (i in 1:23) {
  sdata = data[data[,13] == i,1:6]
  corr[i] = sum(sdata==3)
  dd[i,] = apply(sdata==3,2,mean)
}
spl1 = corr >= median(corr)

corr=array()
dd=array(NA,dim=c(23,6))
for (i in 1:23) {
  sdata = data[data[,13] == i,7:12]
  corr[i] = sum(sdata==3)
  dd[i,] = apply(sdata==3,2,mean)
}
spl2 = corr >= median(corr)

par(mfrow=c(2,3))
T1datag = dd[spl1,]
barplot(apply(T1datag,2,mean), main ="Correct good performers", ylim = c(0,1))
T1datab = dd[spl1==0,]
barplot(apply(T1datab,2,mean), main ="Correct bad performers", ylim = c(0,1))
barplot(apply(T1datag,2,mean)-apply(T1datab,2,mean), main ="Correct bad performers", ylim = c(0,1))

T2datag = dd[spl2,]
barplot(apply(T2datag,2,mean), main ="Correct good performers", ylim = c(0,1))
T2datab = dd[spl2==0,]
barplot(apply(T2datab,2,mean), main ="Correct bad performers", ylim = c(0,1))
barplot(apply(T2datag,2,mean)-apply(T2datab,2,mean), main ="Correct bad performers", ylim = c(0,1))

corr = apply(data[,1:6]==3,1,sum)

par(mfrow=c(1,3))
mspl1 = corr >= median(corr)
T1datag = data[mspl1,1:6]==3
barplot(apply(T1datag,2,mean), main ="Correct", ylim = c(0,1))

mspl2 = corr< median(corr)
T1datab = data[mspl2,1:6]==3
barplot(apply(T1datab,2,mean), main ="Correct", ylim = c(0,1))

mspl2 = corr< median(corr)
T1data = data[mspl2,1:6]==3
barplot(apply(T1datag,2,mean)-apply(T1datab,2,mean), main ="Correct", ylim = c(0,1))

