# Script to read the behavioral data out of the mat-files and aggregate them for the modeling
rm(list=ls())

library(rmatio)

setwd("~/Dropbox/AA_Neurometric/ANALYSES/SimplespanAnalyses/Modeling/")
dataDir = '/Volumes/methlab_data/Neurometric/2018/Test-Retest'

# get the files with behavioral results

folders = list.dirs(dataDir,recursive = FALSE)
R = array(NA, dim=c(length(folders) ,4))

q = 1
for (i in 1:length(folders)) {
  
  ID = substr(folders[i], nchar(folders[i])-2 , nchar(folders[i]))
  file = paste(folders[i], "/" ,ID ,"_SIS.mat",sep="")
  
  if (file.exists(file)) {
    
    data = read.mat(file)
    
    reference = rep(seq(from = 1, to = 6),16)
    responses = as.vector(t(data$FullResults))
    
    # check if correct, wrong pos, wrong
    correct = sum(responses  == reference)
    wrong = sum(responses > 6)
    wrongPosition = length(responses) - correct - wrong
    
    # save it as 
    R[q,] = c(ID,correct,wrongPosition,wrong)
    q = q + 1
  }
  else {
    print(folders[i])
  }
}

R = R[complete.cases(R),]

data = write.table(R,"/Users/apedroni/Dropbox/AA_Neurometric/ANALYSES/SimplespanAnalyses/Modeling/SIS_2018.csv",sep=",",row.names = F,quote = F  )
