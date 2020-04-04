#rm(list=ls(all=TRUE))

source('../R Functions/RF functions 20JAN19.R')

#Restrict to nights analyzed:
load(file="Data/nights/nights.RData")
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

#Load RMRF
load(file="Data/RMRF results/arDep_ex.RData")


# Make a prediction for each night using the OOB sample
pred_ar_ex <- predictRF(arDep_ex, nights, prediction = "overall", checkCases = TRUE)

# Classify each night as either "Low Risk" (less than avg) or "High Risk" (greater than avg)
nights$risk <- ifelse(pred_ar_ex < mean(nights$cgmHypo60_30min), "1.Low", "2.High")

# In previous submission, but later removed and not shown.
library(plyr)
ddply(nights, ~ risk, summarise, nNights=length(risk), percHypoNights=round(mean(cgmHypo60_30min),3))
summary(geem((cgmHypo60_30min_Factor == "Yes") ~ risk, data = nights, id = DeidentID, corstr = "exchangeable"))

