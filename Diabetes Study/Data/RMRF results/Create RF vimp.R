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


#ntrees <- 2000  #Manuscript creates 2000 trees on supercomputers
ntrees <- 50  #Note: the RMRF analysis takes a while, even with parallel computing
seed <- 482041  #Manuscript uses 3 different seeds and combines the observed variable importance

arIndep <- growRF_Parallel(ntrees=ntrees, formula=cgmHypo60_30min_Factor ~ ageAtEnroll + Gender + HbA1CTest + diabDuration + dailyBasal_Ins + dailyIns_U_kg_d +
                             bedtimeBG + Snack + iob_dailyIns + cgmRateChgGrp + sysStartTmGrp + exer_Intensity + daytimeHypo,
                           data=nights, search="ar", method="class", split="gini",
                           mtry=1, nsplit=1, minsplit=60, minbucket=20, maxdepth=10, sampleMethod='subsample', iseed=seed)
vimp_arIndep_auc <- vimpRF(arIndep, nights, checkCases=TRUE, vimpStat="auc")
vimp_arIndep_misclass <- vimpRF(arIndep, nights, checkCases=TRUE, vimpStat="misclass")


arDep_ar <- growRF_Parallel(ntrees=ntrees, formula=cgmHypo60_30min_Factor ~ ageAtEnroll + Gender + HbA1CTest + diabDuration + dailyBasal_Ins + dailyIns_U_kg_d +
                              bedtimeBG + Snack + iob_dailyIns + cgmRateChgGrp + sysStartTmGrp + exer_Intensity + daytimeHypo | DeidentID,
                            data=nights, search="ar", method="class", split="gini",
                            mtry=1, nsplit=1, minsplit=60, minbucket=20, maxdepth=10, sampleMethod='subsampleByID', corstr="ar1", iseed=seed)
vimp_arDep_ar_auc <- vimpRF(arDep_ar, nights, checkCases=TRUE, vimpStat="auc")
vimp_arDep_ar_misclass <- vimpRF(arDep_ar, nights, checkCases=TRUE, vimpStat="misclass")


arDep_ex <- growRF_Parallel(ntrees=ntrees, formula=cgmHypo60_30min_Factor ~ ageAtEnroll + Gender + HbA1CTest + diabDuration + dailyBasal_Ins + dailyIns_U_kg_d +
                              bedtimeBG + Snack + iob_dailyIns + cgmRateChgGrp + sysStartTmGrp + exer_Intensity + daytimeHypo | DeidentID,
                            data=nights, search="ar", method="class", split="gini",
                            mtry=1, nsplit=1, minsplit=60, minbucket=20, maxdepth=10, sampleMethod='subsampleByID', corstr="exchangeable", iseed=seed)
vimp_arDep_ex_auc <- vimpRF(arDep_ex, nights, checkCases=TRUE, vimpStat="auc")
vimp_arDep_ex_misclass <- vimpRF(arDep_ex, nights, checkCases=TRUE, vimpStat="misclass")

diabetes_vimp <- as.data.frame(rbind(vimp_arIndep_auc$scaledVI, vimp_arIndep_misclass$scaledVI,
                                     vimp_arDep_ar_auc$scaledVI, vimp_arDep_ar_misclass$scaledVI,
                                     vimp_arDep_ex_auc$scaledVI, vimp_arDep_ex_misclass$scaledVI))

diabetes_vimp$vimp <- c("arIndep_auc","arIndep_misclass","arDep_ar_auc","arDep_ar_misclass","arDep_ex_auc","arDep_ex_misclass")

#save(arDep_ex, file="arDep_ex.RData")  #Save RF for risk profiles analysis
#write.table(diabetes_vimp, file = "diabetes_vimp_seed1.csv", sep = ",", row.names=FALSE)
