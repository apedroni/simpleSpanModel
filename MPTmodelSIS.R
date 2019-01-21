## Simple Span model for Neurometric data
# this script loads data and submits it to the model "SimpleSpan1noF.txt"
# The individual parameters are stored in JAGSModel$BUGSoutput$sims.list 


rm(list=ls())
setwd("~/Dropbox/AA_Neurometric/ANALYSES/SimplespanAnalyses/Modeling/")
#library(ggplot2)
require(R2jags) 
#library(tidyr)
#library(plyr)
#library(reshape2)

traceP <- T  # trace individual's parameters

## read and prepare data
#data = read.table("SimpleSpan4Modeling.csv",sep=',',header = TRUE)
R = read.table("SIS_2018.csv",sep=',',header = TRUE)

R <- R[,-1]

N = rowSums(R)
Nsub = length(R[,1])

## hierarchical model without F

  modelname = "SimpleSpan1noF.txt"
  
  # this is the base rate of categories
  ch <- c(1, 5, 6)
  
  dataList = list(R = R,Nsub = Nsub, N = N, ch = ch)
  
  subjectparameters = c("C" , "A"  )  # subject parameter distributions
  muparameters =      c( "muA" , "muC" ) # mean of group level posterior distirbution
  sgparameters =      c("sgA" ,  "sgC" )  # standard deviation of group level posterior distirbution
  parameters = c(muparameters,deltaparameters,subjectparameters)

  if (traceP == T) parameters <- c(parameters, "P") # Predicted value of the model
  
  JAGSModel <- jags.parallel(data=dataList, inits=NULL, parameters.to.save=parameters, model.file = modelname,
                             n.chains = 3, n.iter = 50000, n.burnin = 1000, n.thin = 5)
  
  P = JAGSModel$BUGSoutput$sims.list$P
  JAGSModel$BUGSoutput$DIC
  
  # Plot how well the predictions are
  Predicted = array(NA,dim=c(Nsub,3))  
  Predicted[,1] = apply(P[,,1],2,mean)
  Predicted[,2] = apply(P[,,2],2,mean)
  Predicted[,3] = apply(P[,,3],2,mean)

  plot(Predicted[,1],R[,1])
  plot(Predicted[,2],R[,2])
  plot(Predicted[,3],R[,3])
 
  cor(Predicted[,1],R[,1])
  cor(Predicted[,2],R[,2])
  cor(Predicted[,3],R[,3])