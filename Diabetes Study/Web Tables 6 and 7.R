#rm(list=ls(all=TRUE))

library(plyr)

# Load nights dataset.  Restrict to control nights during clinical trial phase with at least 6 hrs of CGM data
load(file="Data/nights/nights.RData")
nights <- nights[nights$psoPtPhase=='C)Clinical Trial Phase' & nights$AlgoActive == "OFF" & nights$valid6Hrs, ]

# Function to categorize a continuous variable
createGrp <- function(var, cutpoints, forward=TRUE){
  temp <- numeric(length(var))
  if (forward){
    temp[var < cutpoints[1]] <- paste0("A:<",format(cutpoints[1], scientific=FALSE))
    temp[var >= cutpoints[length(cutpoints)]] <- paste0(rawToChar(as.raw(65+length(cutpoints))),":>=",format(cutpoints[length(cutpoints)],scientific=F))
    for(i in 2:length(cutpoints)){
      temp[var >= cutpoints[i-1] & var < cutpoints[i]] <- paste0(rawToChar(as.raw(64+i)),":[",format(cutpoints[i-1],scientific=F),",",format(cutpoints[i], scientific=FALSE),")")
    }
  } else {
    temp[var < cutpoints[1]] <- paste0(rawToChar(as.raw(65+length(cutpoints))),":<",format(cutpoints[1], scientific=FALSE))
    temp[var >= cutpoints[length(cutpoints)]] <- paste0("A:>=",format(cutpoints[length(cutpoints)],scientific=F))
    for(i in 2:length(cutpoints)){
      temp[var >= cutpoints[i-1] & var < cutpoints[i]] <- paste0(rawToChar(as.raw(66+length(cutpoints)-i)),":[",format(cutpoints[i-1],scientific=F),",",format(cutpoints[i], scientific=FALSE),")")
    }
  }
  return(temp)
}

nights$ageGrp <- createGrp(nights$ageAtEnroll,c(4, 13, 22, 46))
nights$HbA1cTestGrp <- createGrp(nights$HbA1CTest,c(6.5, 7.5, 8.0, 8.6), forward=FALSE)
nights$diabDurationGrp <- createGrp(nights$diabDuration,c(5, 10, 20))
nights$dailyBasal_InsGrp <- createGrp(nights$dailyBasal_Ins,c(0.4, 0.5, 0.6))
nights$dailyIns_U_kg_dGrp <- createGrp(nights$dailyIns_U_kg_d,c(0.4, 0.8))
nights$bedtimeBGGrp <- createGrp(nights$bedtimeBG,c(100, 130, 180), forward=FALSE)
nights$iob_dailyInsGrp <- createGrp(nights$iob_dailyIns,c(0, 2.5, 5, 7.5))
nights$Overall="Overall" #Create overall variable

### Web Table 5. Patient-Level Factors ###
patientVars <- c("Overall", "ageGrp","Gender","HbA1cTestGrp","diabDurationGrp","dailyBasal_InsGrp","dailyIns_U_kg_dGrp")
patientMetrics <- lapply(patientVars, function(x){ddply(.data=nights, c(x,"AlgoActive"), summarize,
                                                        nNights = length(cgmHypo60_30min), pctHypoNights = mean(cgmHypo60_30min))})
patientMetrics <- do.call("rbind.fill", patientMetrics)
patientMetrics$factor <- apply(patientMetrics, 1, function(x) { patientVars[!is.na(x[patientVars])]})
patientMetrics$level <- apply(patientMetrics, 1, function(x) { x[patientVars][!is.na(x[patientVars])]})
patientMetrics <- patientMetrics[c("factor","level","nNights","pctHypoNights")]
patientMetrics$pctHypoNights <- round(patientMetrics$pctHypoNights*100, digits=0)
patientMetrics

### Night-level factors ###
nightVars <- c("Overall","bedtimeBGGrp","Snack","iob_dailyInsGrp","cgmRateChgGrp","sysStartTmGrp","exer_Intensity","daytimeHypo")
nightMetrics <- lapply(nightVars, function(x){ddply(.data=nights, c(x,"AlgoActive"), summarize,
                                                    nNights = length(cgmHypo60_30min), pctHypoNights = mean(cgmHypo60_30min))})
nightMetrics <- do.call("rbind.fill", nightMetrics)
nightMetrics$factor <- apply(nightMetrics, 1, function(x) { nightVars[!is.na(x[nightVars])]})
nightMetrics$level <- apply(nightMetrics, 1, function(x) { x[nightVars][!is.na(x[nightVars])]})
nightMetrics <- nightMetrics[c("factor","level","nNights","pctHypoNights")]
nightMetrics$pctHypoNights <- round(nightMetrics$pctHypoNights*100, digits=0)
nightMetrics


