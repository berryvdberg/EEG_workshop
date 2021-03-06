# clear workspace
rm(list = ls(all = TRUE))
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages) 
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}


# load required packages, install if missing
list.of.packages <- c("ggplot2", "lme4", 'car' , 'plyr' ,'dplyr','foreign', 'emmeans', 'multcomp',"lmerTest",'data.table','MuMIn')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)


# LOOP OVER SUBJECTS AND LOAD DATA 
modData = c()
for(i in c(2:4,6, 8:16)){
  dataIn <- read.csv(paste('../data/dataExport/',i,'_trialStruct.csv', sep = ""))
  dataIn <- data.table(dataIn)
  dataIn <- dataIn[(nrow(dataIn) - 642) :nrow(dataIn),]
  dataIn[,trialNumber := 1:nrow(dataIn)]
  modData <- rbind(modData,dataIn)
}

modData[,reactionTime := responseLatency - targetLatency]

### descriptives
modData[artThresh == 1 | artheog == 1,.(.N),.(cueType)]


### RTS ---- UNPOOLED -----
mUnpooled <- lmList(reactionTime ~ factor(cueType) | subjectID, modData[cueType != 6 & reactionTime >100 & reactionTime < 1200 &responseType %in% 3:4])
tmp <- coef(mUnpooled)
t.test(tmp$`factor(cueType)9`)

tmp <- modData[cueType != 6 & reactionTime >100 & reactionTime < 1200 &responseType %in% 3:4, 
        .(RT.mean = mean(reactionTime)),.(subjectID,cueType)]


t.test(tmp[cueType == 9,RT.mean],tmp[cueType == 7,RT.mean], paired = T)



### RTS ---- POOLED -----
mPooled <-lm(reactionTime ~ factor(cueType), modData[cueType != 6 & reactionTime >100 & reactionTime < 1200 &responseType %in% 3:4])
summary(mPooled)


### RTS ---- PARTIAL POOLED -----
mPartialPooled <-lmer(reactionTime ~ as.factor(cueType) + (1|subjectID), 
          data = modData[artThresh == 0 & cueType != 6 & reactionTime >100 & reactionTime < 1200 & responseType %in% 3:4]) 


summary(mPartialPooled)
ranef(mPartialPooled)

mPartialPooled <-lmer(reactionTime ~ as.factor(cueType) + (1|subjectID) + (0+cueType|subjectID), 
                      data = modData[cueType != 6 & reactionTime >100 & reactionTime < 1200 & responseType %in% 3:4 ]) 
summary(mPartialPooled)
ranef(mPartialPooled)

## lets construct polynomials... 
plot(poly(1:100,11)[,2],type = 'l')



