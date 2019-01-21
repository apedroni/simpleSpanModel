## Simple Span model with two timepoints 
# this script can be run with the 2017 pilot data for the test-retest study
# note that the different variations of the models are in archive.

rm(list=ls())
setwd("~/Dropbox/AA_Neurometric/ANALYSES/SimplespanAnalyses/Modeling/")
library(ggplot2)
require(R2jags) 
library(tidyr)
library(plyr)
library(reshape2)

hier =TRUE # do hierachical bayesian model
noF = TRUE # do include F as a intrusion parameter (not applicable in 2018 dataset)
traceP <- T  # trace individual's parameters

## read and prepare data
data = read.table("SimpleSpan4Modeling.csv",sep=',',header = TRUE)

colnames(data) = c("ID","correct_T1","wrong_pos_T1","previous_T1","not_remembered_T1","correct_T2","wrong_pos_T2","previous_T2","not_remembered_T2")

# put the dataframe to an array
R1 = array(NA, dim=c(length(data[,1]) ,4))
R1[,1] = data$correct_T1
R1[,2] = data$wrong_pos_T1
R1[,3] = data$previous_T1
R1[,4] = data$not_remembered_T1

R2 = array(NA, dim=c(length(data[,1]) ,4))
R2[,1] = data$correct_T2
R2[,2] = data$wrong_pos_T2
R2[,3] = data$previous_T2
R2[,4] = data$not_remembered_T2

R <- rbind(R1,R2)

N = rowSums(R)
Nsub = length(R[,1])

## hierarchical model without F
if (noF == TRUE) {

    # merge wrong position and previous to wrong position 
  R[,2] <- R[,2]+R[,3]
  R <- R[,-3]
  
  modelname = "SimpleSpan1noF.txt"
  
  # this is the base rate of categories
  ch <- c(1, 5, 6)
  
  dataList = list(R = R,Nsub = Nsub, N = N, ch = ch)
  
  subjectparameters = c("C" , "A"  )
  deltaparameters = c("DmuA", "DmuC")
  muparameters =      c( "muA" , "muC" )
  sgparameters =      c("sgA" ,  "sgC" )  
  parameters = c(muparameters,deltaparameters,subjectparameters)

  if (traceP == T) parameters <- c(parameters, "P")
  
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

  ## test retest DATA
  cor(R[1:23,1],R[24:46,1],method="spearman")
  cor(R[1:23,2],R[24:46,2],method="spearman")
  cor(R[1:23,3],R[24:46,3],method="spearman")
  plot(R[1:23,1],R[24:46,1])
  plot(R[1:23,2],R[24:46,2])
  plot(R[1:23,3],R[24:46,3])

  ## test retest Model Prediction
  cor(Predicted[1:23,1],Predicted[24:46,1])
  cor(Predicted[1:23,2],Predicted[24:46,2])
  cor(Predicted[1:23,3],Predicted[24:46,3])

  ## test retest Model Parameters
  A = JAGSModel$BUGSoutput$sims.list$A
  hist(A[,1])
  A = apply(A,2,median)
  cor(A[1:23],A[24:46],method="spearman")
  plot(A[1:23],A[24:46])
  
  C = JAGSModel$BUGSoutput$sims.list$C
  C = apply(C,2,median)
  cor(C[1:23],C[24:46],method="spearman")
  plot(C[1:23],C[24:46])
  
  plot(Predicted[1:23,1],Predicted[24:46,1])
  plot(Predicted[1:23,2],Predicted[24:46,2])
  plot(Predicted[1:23,3],Predicted[24:46,3])

}


## hierarchical model with F
if (noF == FALSE) {
## full model################################# RUN THE CHAINS  ##########################

# this is to calculate the base rate: 1 possible correct, 6 possible wrong, 6 possible wrong position
ch <- c(1, 5, 3, 3)

dataList = list(R = R,Nsub = Nsub, N = N, ch = ch)
traceP <- T  # trace individual's parameters

## with F and hierarchical
if (hier == 1) {
  subjectparameters = c("C" , "A" , "F" )
  deltaparameters = c("DmuA", "DmuC", "DmuF")
  muparameters =      c( "muA" , "muC" , "muF")
  sgparameters =      c("sgA" ,  "sgC" , "sgF")  
  parameters = c(muparameters,deltaparameters,subjectparameters)
  modelname = "SimpleSpan1.txt"
}    else {
  subjectparameters = c("C" , "A" , "F" )
  parameters = subjectparameters  
  modelname = "SimpleSpan1nh.txt"
}

if (traceP == T) parameters <- c(parameters, "P")

JAGSModel <- jags.parallel(data=dataList, inits=NULL, parameters.to.save=parameters, model.file = modelname,
                           n.chains = 3, n.iter = 50000, n.burnin = 1000, n.thin = 5)
JAGSModel$BUGSoutput$DIC

  checkConvergence = F
  if (traceP == T) checkConvergence = F
  if ( checkConvergence ) {
    print(JAGSModel)
    traceplot(JAGSModel)
    codaSamples <- as.mcmc(JAGSModel)
    xyplot(codaSamples)
    densityplot(codaSamples)
    mcmcChain = as.matrix( codaSamples )
    }
  codaData = as.mcmc(JAGSModel)
  codaSamples <- as.matrix(codaData)
  
  P = JAGSModel$BUGSoutput$sims.list$P
  
  
  setwd("~/Dropbox/AA_Neurometric/Paradigmen/SequenceMemory/Results/")
  # Plot how well the predictions are
  Predicted = array(NA,dim=c(Nsub,4))  
  Predicted[,1] = apply(P[,,1],2,mean)
  Predicted[,2] = apply(P[,,2],2,mean)
  Predicted[,3] = apply(P[,,3],2,mean)
  Predicted[,4] = apply(P[,,4],2,mean)
  
  plot(Predicted[,1],R[,1])
  plot(Predicted[,2],R[,2])
  plot(Predicted[,3],R[,3])
  plot(Predicted[,4],R[,4])
  
  ## test retest DATA
  cor(R[1:23,1],R[24:46,1])
  cor(R[1:23,2],R[24:46,2])
  cor(R[1:23,3],R[24:46,3])
  cor(R[1:23,4],R[24:46,4])  
  plot(R[1:23,1],R[24:46,1])
  plot(R[1:23,2],R[24:46,2])
  plot(R[1:23,3],R[24:46,3])
  plot(R[1:23,4],R[24:46,4])    
  
  ## test retest Model Prediction
  cor(Predicted[1:23,1],Predicted[24:46,1])
  cor(Predicted[1:23,2],Predicted[24:46,2])
  cor(Predicted[1:23,3],Predicted[24:46,3])
  cor(Predicted[1:23,4],Predicted[24:46,4])
  
  ## test retest Model Parameters
  A = JAGSModel$BUGSoutput$sims.list$A
  hist(A[,1])
  A = apply(A,2,median)
  cor(A[1:23],A[24:46],method="spearman")
  plot(A[1:23],A[24:46])
  
  C = JAGSModel$BUGSoutput$sims.list$C
  C = apply(C,2,median)
  cor(C[1:23],C[24:46],method="spearman")
  plot(C[1:23],C[24:46])
  
  Ff = JAGSModel$BUGSoutput$sims.list$F
  Ff = apply(Ff,2,median)
  cor(Ff[1:23],Ff[24:46],method="spearman")
  plot(Ff[1:23],Ff[24:46])
  
  plot(Predicted[1:23,1],Predicted[24:46,1])
  plot(Predicted[1:23,2],Predicted[24:46,2])
  plot(Predicted[1:23,3],Predicted[24:46,3])
  plot(Predicted[1:23,4],Predicted[24:46,4])  
}
