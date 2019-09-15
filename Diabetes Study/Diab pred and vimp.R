#rm(list=ls(all=TRUE))

library(lme4)

source('../R Functions/RF functions 20JAN19.R')
load(file="Data/nights/nights.RData")
nights <- nights[nights$psoPtPhase=='C)Clinical Trial Phase' & nights$valid6Hrs & nights$AlgoActive=="OFF", ]

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

scaleVar <- function(var){ (nightsScaled[[var]] - mean(nightsScaled[[var]]))/sd(nightsScaled[[var]]) }

nightsScaled <- nights
nightsScaled[["ageAtEnroll"]] <- scaleVar("ageAtEnroll")
nightsScaled[["HbA1CTest"]] <- scaleVar("HbA1CTest")
nightsScaled[["diabDuration"]] <- scaleVar("diabDuration")
nightsScaled[["dailyBasal_Ins"]] <- scaleVar("dailyBasal_Ins")
nightsScaled[["dailyIns_U_kg_d"]] <- scaleVar("dailyIns_U_kg_d")
nightsScaled[["bedtimeBG"]] <- scaleVar("bedtimeBG")
nightsScaled[["iob_dailyIns"]] <- scaleVar("iob_dailyIns")


# Function to perform 5-fold CV using either GLMM or RMRF #
kFoldCV <- function(currentVars, inData, method=c("glmm", "rmrf")) {
  
  method <- match.arg(method, c("glmm", "rmrf"))
  
  # Split data into 5 groups with equal number of IDs
  nIds <- length(unique(inData[["DeidentID"]]))
  randomIDs <- sample(unique(inData[["DeidentID"]]))
  
  # nIds isn't divisible by 5; make it as close to balanced as possible
  r  <- rep(1:5, each=floor(nIds/5))[1:nIds]
  addNum <- 1
  for (i in 1:4) {
    if(any(is.na(r))) {
      r[is.na(r)][1] <- addNum
      addNum <- addNum + 1
    }
  }
  
  foldData <- lapply(1:5, function(x) inData[inData[["DeidentID"]] %in% randomIDs[r == x], ])
  #length(unique(foldData[[2]]$DeidentID))
  
  auc.test <- rep(NA, 5)
  for (k in 1:5) {
    test <- foldData[[k]]
    training <- do.call("rbind", foldData[-k])
    
    if (method == "glmm") {
      modTrained <- glmer(currentVars, family=binomial, data=training, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e4)))
      pred.test <- predict(modTrained, newdata = test, type='response', re.form=NA)
      auc.test[k] <- auc(test[["cgmHypo60_30min"]], pred.test, direction="<")
    } else if (method == "rmrf") {
      ntrees <- 50  #Manuscript creates 2000 trees on supercomputers
      rmrfTrained <- growRF_Parallel(ntrees=ntrees, formula=currentVars,
                                     data=training, search="ar", method="class", split="gini",
                                     mtry=1, nsplit=1, minsplit=60, minbucket=20, maxdepth=10, sampleMethod='subsampleByID', corstr="exchangeable",
                                     iseed=seed)
      pred.test <- predictRF(rmrfTrained, test, prediction = "overall", checkCases = FALSE)
      auc.test[k] <- auc(test[["cgmHypo60_30min"]], pred.test, direction="<")
    }
  }
  
  return(mean(auc.test))
}


# Function to perform backwards selection.  Starts with the full model and sequentially removes #
# variables based on either the largest p-value or the smaller model with the greatest AUC #
backwardSelect <- function(step=1, currentVars, currentMod, maxAUC = NA, maxAUCVars = NA,
                           inData = nights, selection = c("p-value", "CV")) {
  
  selection <- match.arg(selection, c("p-value", "CV"))
  
  # Consider removing variable with the lowest F-value
  temp <- anova(currentMod)
  removeVar <- rownames(temp)[which.min(temp$F)]
  newVars <- as.formula(paste0("cgmHypo60_30min ~ ", 
                               paste(all.vars(currentVars[[3]])[all.vars(currentVars[[3]]) != removeVar & all.vars(currentVars[[3]]) != "DeidentID"], collapse = " + "),
                               " + (1 | DeidentID)"))
  
  if (selection == "p-value") {
    newMod <- glmer(newVars, family=binomial, data=inData, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e4)))
    pvalue <- anova(currentMod, newMod)$Pr[2]
    
    # If p-value is >0.05, then repeat function; otherwise, obtained final model
    if (pvalue >0.05) {
      print(paste0("STEP ", step, ": Removed ", removeVar, ".  P-value is ", round(pvalue, 3)))
      backwardSelect(step = step + 1, currentVars = newVars, currentMod = newMod, maxAUC = NA, maxAUCVars = NA,
                     inData = inData, selection = selection)
    } else {
      print(paste0("STEP ", step, ": All p-values >0.05"))
      return(currentVars)
    }
    
  } else if (selection == "CV") {
    newAUC <- kFoldCV(currentVars=newVars, inData=inData, method="glmm")
    
    # If new AUC has increased, update AUC
    if (maxAUC <= newAUC) {
      maxAUC <- newAUC
      maxAUCVars <- newVars
    }
    print(paste0("STEP ", step, ": Removed ", removeVar, ". New AUC is ", round(newAUC, 3), ".  Max AUC is ", round(maxAUC, 3)))
    
    # Once you get back to the intercept-only model, then stop
    if (length(all.vars(newVars)) == 3) { return(maxAUCVars) }
    
    # Need newMod to estimate the highest p-value
    newMod <- glmer(newVars, family=binomial, data=inData, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e4)))
    
    backwardSelect(step = step + 1, currentVars = newVars, currentMod = newMod, maxAUC = maxAUC, maxAUCVars = maxAUCVars,
                   inData = inData, selection = selection)
    
  }
}

# Function to perform step-wise selection.  Starts with the intercept-only model and sequentially adds #
# variables based on smallest p-value.  Performs approximated p-value calculation to show all variables #
# are still significant #
stepwiseSelect <- function(step=1, currentVars, remainingVars, currentMod, inData = nights) {
  
  # Consider adding each variable
  pvalues <- rep(NA, length(remainingVars))
  for (i in seq.int(remainingVars)) {
    newVars <- update(currentVars, paste0("~ . + ", remainingVars[i]))
    newMod <- glmer(newVars, family=binomial, data=inData, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e4)))
    pvalues[i] <- anova(currentMod, newMod)$Pr[2]
  }
  
  if (any(pvalues < 0.05)) {
    addedVar <- remainingVars[which.min(pvalues)]
    remainingVars <- remainingVars[-which.min(pvalues)] 
    newVars <- update(currentVars, paste0("~ . + ", addedVar))
    newMod <- glmer(newVars, family=binomial, data=inData, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e4)))
    print(paste0("STEP ", step, ": Added ", addedVar, ".  P-value is ", round(min(pvalues), 3)))
    
    # Checking if p-values are still significant is tricky with GLMM
    # Simplest is to give an approximate p-value using F table showing the variables are highly significant
    pvalues <- 1-pf(anova(newMod)$F, anova(newMod)$Df, nrow(inData) - sum(anova(newMod)$Df))
    print(paste0("STEP ", step+1, ": Highest p-value is approximately ", round(max(pvalues), 3)))
    
    stepwiseSelect(step = step + 2, currentVars = newVars, remainingVars = remainingVars, currentMod = newMod, inData = inData)
  } else {
    print(paste0("STEP ", step, ": No p-values <0.05"))
    return(currentVars)
  }
}

# Use either scaled data or original data
inData <- nightsScaled    # Manuscript considers nightsScaled; using nights yields similar results
# Determine final model either by p-value or AUC from CV
selection <- "p-value"  # Using CV instead of p-value slightly improves accuracy, but very similar results


########## GLMM ##########

##### VIMP #####

allVars <- c("ageAtEnroll", "Gender", "HbA1CTest", "diabDuration", "dailyBasal_Ins", "dailyIns_U_kg_d",
             "bedtimeBG", "Snack", "iob_dailyIns", "cgmRateChgGrp", "sysStartTmGrp", "exer_Intensity", "daytimeHypo")
vimpGLMM <- rep(NA, length(allVars))
names(vimpGLMM) <- allVars
for (i in seq.int(allVars)) {
  glmmMod <- glmer(as.formula(paste0("cgmHypo60_30min ~ ", allVars[i], " + (1 | DeidentID)")),
                   data=inData, family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e4)))
  vimpGLMM[i] <- summary(glmmMod)$coef[2,1]
}
vimpGLMM[order(abs(vimpGLMM), decreasing=TRUE)]


### Backwards Variable Selection ###

currentVars <- as.formula(cgmHypo60_30min ~ ageAtEnroll + Gender + HbA1CTest + diabDuration + dailyBasal_Ins + dailyIns_U_kg_d +
                            bedtimeBG + Snack + iob_dailyIns + cgmRateChgGrp + sysStartTmGrp + exer_Intensity + daytimeHypo + (1 | DeidentID))
currentMod <- glmer(currentVars, data=inData, family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e4)))
set.seed(470174)

if (selection == "p-value") { maxAUC <- NA
} else if (selection == "CV") { maxAUC <- kFoldCV(currentVars=currentVars, inData=inData, method="glmm") }

backwardGLMM <- backwardSelect(step=1, currentVars=currentVars, currentMod=currentMod, maxAUC=maxAUC, maxAUCVars=currentVars,
                               inData=inData, selection=selection)


### Step-wise Variable Selection ###

currentVars <- as.formula(cgmHypo60_30min ~ (1 | DeidentID))
remainingVars <- c("ageAtEnroll", "Gender", "HbA1CTest", "diabDuration", "dailyBasal_Ins", "dailyIns_U_kg_d",
                   "bedtimeBG", "Snack", "iob_dailyIns", "cgmRateChgGrp", "sysStartTmGrp", "exer_Intensity", "daytimeHypo")
currentMod <- glmer(cgmHypo60_30min ~ (1 | DeidentID), data=inData, family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e4)))
set.seed(470174)

stepwiseGLMM <- stepwiseSelect(step=1, currentVars=currentVars, remainingVars=remainingVars, currentMod=currentMod, inData=inData)


##### AUC for GLMM based on CV #####

# backwardGLMM and stepwiseGLMM when using classical p-value approach
#glmmMod <- as.formula(cgmHypo60_30min ~ bedtimeBG + exer_Intensity + daytimeHypo + (1 | DeidentID))

# backwardGLMM when using CV
#glmmMod <- as.formula(cgmHypo60_30min ~ bedtimeBG + iob_dailyIns + cgmRateChgGrp + exer_Intensity + daytimeHypo + (1 | DeidentID))

# Using the same model as RMRF
#glmmMod <- as.formula(cgmHypo60_30min ~ bedtimeBG + iob_dailyIns + sysStartTmGrp + exer_Intensity + daytimeHypo + (1 | DeidentID))

# Use three different seeds and take the average
set.seed(482041)
AUCseed1 <- kFoldCV(currentVars=glmmMod, inData=inData, method="glmm")
set.seed(482042)
AUCseed2 <- kFoldCV(currentVars=glmmMod, inData=inData, method="glmm")
set.seed(482043)
AUCseed3 <- kFoldCV(currentVars=glmmMod, inData=inData, method="glmm")

format(round(mean(c(AUCseed1, AUCseed2, AUCseed3)),3), nsmall=3)



########## RMRF ##########

##### VIMP already created in separate program #####

##### AUC for GLMM based on CV #####

# Takes a long time to do RMRF with CV #
currentVars <- cgmHypo60_30min_Factor ~ ageAtEnroll + Gender + HbA1CTest + diabDuration + dailyBasal_Ins + dailyIns_U_kg_d +
  bedtimeBG + Snack + iob_dailyIns + cgmRateChgGrp + sysStartTmGrp + exer_Intensity + daytimeHypo | DeidentID

seed <- 482041; set.seed(seed)
AUCseed1 <- kFoldCV(currentVars=currentVars, inData=inData, method="rmrf")
seed <- 482042; set.seed(seed)
AUCseed2 <- kFoldCV(currentVars=currentVars, inData=inData, method="rmrf")
seed <- 482043; set.seed(seed)
AUCseed3 <- kFoldCV(currentVars=currentVars, inData=inData, method="rmrf")

format(round(mean(c(AUCseed1, AUCseed2, AUCseed3)),3), nsmall=3)



##### AUC for GLMM based on OOB #####

# Using OOB is much faster since RMRF already calculated.  Also don't need to calculate 5 RMRF algorithms
load(file="Data/RMRF results/arDep_ex.RData")
# Make a prediction for each night using the OOB sample
pred_ar_ex <- predictRF(arDep_ex, nights, prediction = "overall", checkCases = TRUE)
auc(nights$cgmHypo60_30min, pred_ar_ex, direction="<")



