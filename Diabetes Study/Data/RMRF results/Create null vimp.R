#rm(list=ls(all=TRUE))

source('../../../R Functions/RF functions 20JAN19.R')

#Restrict to nights analyzed:
load(file="../nights/nights.RData")
nights <- nights[nights$psoPtPhase=='C)Clinical Trial Phase' & nights$valid6Hrs & nights$AlgoActive=="OFF",]

#Extract dates from SysStartDtTm.  Map times before noon to previous day
currentDt <- function(dttm){
  as.Date(dttm, tz="America/Los_Angeles") - (strftime(dttm, format="%H:%M:%S") < '12:00:00')
}
nights$nightNo <- unlist(sapply(unique(nights$DeidentID), function(id){currentDt(nights$SysStartDtTm[nights$DeidentID==id]) - currentDt(nights$SysStartDtTm[nights$DeidentID==id][1]) + 1}))

nights <- nights[with(nights, order(DeidentID, nightNo, SysStartDtTm)), ]
nights$nCGMHypo60_DaytimeReg <- nights$nCGMHypo60_Daytime
nights$nCGMHypo60_DaytimeReg[!nights$validDay] <- -999

#Convert character variables to factors:
nights$cgmHypo60_30min_Factor <- as.factor(ifelse(nights$cgmHypo60_30min,"Yes","No"))
nights$exer_Intensity <- as.factor(nights$exer_Intensity)
nights$cgmRateChgGrp <- as.factor(nights$cgmRateChgGrp)
nights$sysStartTmGrp <- as.factor(nights$sysStartTmGrp)
nights$daytimeHypo <- as.factor(nights$daytimeHypo)
nights$cgmRateChgGrp <- as.factor(nights$cgmRateChgGrp)
nights$sysStartTmGrp <- as.factor(nights$sysStartTmGrp)
dim(nights)


permute_vimp <- function(sim) {
  library(rpart)
  library(partykit)
  library(geeM)
  library(pROC)
  
  permute.response <- nights
  permute.response$cgmHypo60_30min_Factor <- permute.response$cgmHypo60_30min_Factor[sample.int(length(permute.response$cgmHypo60_30min_Factor))]
  
  #Manuscript also reports results when ignoring the dependence and using Gini gain.
  #To do this, simply remove the "| DeidentID" in the formula, and change sampleMethod="subsample"
  permute.ar <- growRF(ntrees=2000, formula=cgmHypo60_30min_Factor ~ ageAtEnroll + Gender + HbA1CTest + diabDuration + dailyBasal_Ins + dailyIns_U_kg_d +
                          bedtimeBG + Snack + iob_dailyIns + cgmRateChgGrp + sysStartTmGrp + exer_Intensity + daytimeHypo | DeidentID,
                        data=permute.response, search="ar", method="class", split="gini",
                        mtry=1, nsplit=1, minsplit=60, minbucket=20, maxdepth=10, sampleMethod='subsampleByID', corstr="exchangeable")
  permute.auc <- vimpRF(permute.ar, permute.response, checkCases=TRUE, vimpStat="auc")
  permute.misclass <- vimpRF(permute.ar, permute.response, checkCases=TRUE, vimpStat="misclass")
  return(list(permute.auc = permute.auc$scaledVI, permute.misclass = permute.misclass$scaledVI))
}

#nsim <- 100  #Manuscript uses 100 null variable importance simulations on supercomputers
nsim <- 2  #Note: this analysis takes a long time, even with parallel computing

num_cores <- detectCores()-1
cl <- makeCluster(num_cores)
clusterExport(cl,c("nights","expit","findCutpts","geem_mod","pvalue_check","pickCutpt","splitrule",
                   "obj.ttest","obj.MSE","obj.binary","ordinalize",
                   "partitionAR","partition","growTemp","growTree","growRF",
                   "predictTree","predictRF","vimpTree","vimpRF","permute_vimp"))
clusterSetRNGStream(cl=cl, iseed=435701)
nullVimp_dep <- parLapply(cl, 1:nsim, permute_vimp)
stopCluster(cl)

nullVimp_dep <- data.frame(nullVimp_dep)
names(nullVimp_dep) = paste0(c("null_auc", "null_misclass"), rep(1:nsim, each = 2))
nullVimp_dep <- cbind(factor = row.names(nullVimp_dep), nullVimp_dep)

#write.table(nullVimp_dep, file = "nullVimp_dep.csv", sep = ",", row.names=FALSE)

