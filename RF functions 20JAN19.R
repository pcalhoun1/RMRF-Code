###############################################################################
###                                                                         ### 
### This program implements five random forests algorithms.                 ###
### The goal of this program is to provide users with the infrastructure to ###
### grow their own random forest without having to start from scratch.      ###
### Each function is documented to emphasize its use and which algorithms   ###
### use the function.  There are 5 algorithms implemented:                  ###
###                                                                         ###
### RF   - Random Forest.  Standard implementation by Breiman, but allows   ###
###        more control.                                                    ###
### SSS  - Smooth Sigmoid Surrogate Trees.  Replaces indicator function     ###
###        with logistic function.                                          ###
### ER   - Extremely Randomized Trees.  Implementation by Geurts et al.,    ###
###        but allows more control.                                         ###
### AR   - Acceptance-Rejection Trees.  Picks a random variable and then    ###
###        a random cutpoint, and performs quality control check            ###
###        (i.e., p-value < minpvalue)                                      ###
### RMRF - Repeated Measures Random Forest.  Extends AR to handle           ###
###        repeated measurements. Can also do repeated measures RF or ER    ###
###                                                                         ###
### Programmer: Peter Calhoun                                               ###
### Date: 10/2/18                                                           ###
###                                                                         ###
###############################################################################

### Issues ###
# 
# This code is significantly slower than randomForest and randomForestSRC
# R packages.  However, hopefully this code is easier to modify and
# understand for users wanting to build their own random forests.
# This code also implements other random forest algorithms not covered in
# any other R package.  Code cannot handle missing data or handle categorical
# outcomes with more than 2 levels.  Binary response or categorical predictors
# must be converted to a factor.  Tree converted to constparty object.  For
# repeated measures RF, must use predictTree or predictRF to handle
# correlation.

### Common input variables and description ###
# 
# ntrees - Number of trees grown
# formula - Formula to build RF.  Categorical predictors must be listed as
#           a factor.  Response with predictors and id variable (if applicable),
#           e.g., Y ~ x1 + x2 | id
# data - Data to grow RF/tree
# search - Method to search through candidate splits.
#          Options are "exhaustive","sss", "ar"
# method - Response type.
#          If continuous, method="anova"; if binary, method="class"
# split - Impurity measure splitting statistic
#         For continuous, split can be "MSE" or "ttest" (for "sss")
#         For binary, split can be "gini", "entropy" (for "sss"), or "information" 
# mtry - Number of variables randomly selected (for "exhaustive" or "sss")
# nsplit - Number of cutpoints selectioned (for "exhaustive")
# nsplit.random - Logical: indicates if process to select cutpts are random (for "exhaustive")
# minsplit - Number of observations required to continue growing tree
# minbucket - Number of observations required in each child node
# maxdepth - Maximum depth of tree
# a - Sigmoid approximation variable (for "sss")
# sampleMethod - Method to sample learning sample.
#                Options are "bootstrap", "subsample", "subsampleByID", "learning"
# scale.y - Logical: standardize y when creating splits
#           For "sss" to increase stability
# useSearch - Logical: indicates if optimization found by considering many values.
#             For "sss" to ensure minimum is found
# useOptim - Logical: indicates if optimization found by optim() function.
#            For "sss" to ensure minimum is found
# useRpart - Logical: uses rpart to find best split.  Rpart faster than this code.
#            For "exhaustive" with nsplit=NULL
# minpvalue - Minimum p-value to accept split.  For "ar"
# corstr - Covariance structure.  Only used if id variable specified in formula

### For each function, will specify if the function applies to RF, SSS, ER, AR, RMRF, or ALL ###

### Load Functions ###
#rm(list=ls(all=TRUE))
library(rpart)
library(partykit)
library(parallel)
library(geeM)
library(pROC)

### SSS: Logistic (expit) function.  Identical to plogis() function, but slightly faster ###
expit <- function(x) {(tanh(x/2)+1)/2}

### ALL: Finds candidate cutpoints ###
findCutpts <- function(x, minbucket) {
  nX <- length(x)
  #Due to ties, it's possible minbucket cannot be satisfied
  if (sort(x)[minbucket]==sort(x, decreasing=TRUE)[minbucket]) {cutpts=NULL
  } else {
    #Cutpoints considered must satisfy minbucket
    cutpts <- unique(sort(x)[minbucket:(nX-minbucket+1)])
    #Take average of distinct points to determine cutoff (like rpart)
    if(length(cutpts)==1){stop(paste0("Only 1 cutpt??? ", cutpts, x))}
    cutpts <- (cutpts[1:(length(cutpts)-1)]+cutpts[2:length(cutpts)])/2
  }
  return(cutpts)
}

### RMRF: Extract pvalue from geem model ###
geem_mod <- function(yTemp, x, ids, cutpt, corstr, splitStat=c("pvalue","teststat")){
  splitStat <- match.arg(splitStat, c("pvalue","teststat"))
  splitVar <- (x <= cutpt)
  options(warn=2) #Warning message can give strange estimates
  # Try fitting marginal model
  mod1 <- try(summary(geem(yTemp ~ splitVar, family = binomial("logit"),
                           corstr=corstr, id=ids)), silent=TRUE)
  # Try fitting logistic regression model
  mod2 <- try(anova(glm(yTemp ~ splitVar, family = binomial("logit")), test="Chisq"), silent=TRUE)
  # geeM can give erroneously low p-values.  If independence p-value > 100 * robust p-value, use independence stats
  if (inherits(mod2, "try-error")){return(NA)  # If glm fails (unlikely), don't trust geem model
  } else if (inherits(mod1, "try-error") || is.na(mod1$p[2]) || (mod2[["Pr(>Chi)"]][2] > (100*mod1$p[2]))){
    #While p-value comes from a chisquare dist'n, convert to wald statistic to make test statistics comparable
    out <- ifelse(splitStat=="teststat", qnorm(mod2[["Pr(>Chi)"]][2]/2), mod2[["Pr(>Chi)"]][2])
  } else {out <- ifelse(splitStat=="teststat", -abs(mod1$wald[2]), 2*(1-pnorm(abs(mod1$wald[2]))))}
  options(warn=0)
  return(out)
}

### AR and RMRF: Check if cutpoint satisfies p-value threshold ###
pvalue_check <- function(y, x, ids, cutpt, method, minpvalue, corstr){
  method <- match.arg(method, c("anova","class"))
  if (method=="anova") {
    
    #Unfortunately, t-tests splits out errors if constant values or not enough values.  Create new t-test
    new.t.test <- function(x, y){
      #Try regular t-test
      obj <- try(t.test(x, y, alternative="two.sided"), silent=TRUE)
      if (inherits(obj, "try-error")) {
        #If regular t-test fails, try adding noise to see if all values are near constant
        if (length(x) > length(y)) {obj2 <- try(t.test(jitter(x),y,alternative="two.sided"), silent=TRUE)
        } else {obj2 <- try(t.test(x,jitter(y),alternative="two.sided"), silent=TRUE)}
        if (inherits(obj2, "try-error")) {return(NA)
        } else {return(obj2$p.value)}
      } else {return(obj$p.value)}
    }
    
    pvalue <- new.t.test(y[x <= cutpt], y[x > cutpt])
  } else {
    if (!is.null(ids)) {pvalue <- geem_mod(yTemp=(y==levels(y)[2]), x=x, ids=ids, cutpt=cutpt, corstr=corstr, splitStat="pvalue")
    } else {
      obj <- try(suppressWarnings(chisq.test((x <= cutpt),y)), silent=TRUE)
      if (inherits(obj, "try-error")) {pvalue <- NA
      } else {pvalue <- obj$p.value}
    }
  }
  #Check if split has p-value <= minpvalue
  if (is.na(pvalue) || pvalue > minpvalue) {return(NULL)
  } else {return(cutpt)}
}

### AR and RMRF: Pick a random variable and then a random cutpoint.  Check if cutpoint satisfies p-value threshold ###
pickCutpt <- function(allVars, y, ids, method, minsplit, minbucket, nGuesses, minpvalue, corstr){

  # Note: it is possible no candidate splits satisfy p-value threshold.  Thus, code will try nGuesses (default 50) before stopping
  for (guess in 1:nGuesses) {
    v <- sample.int(NCOL(allVars), 1)
    #If variable selected is a factor, pick a random subset of categories.
    #Do not sort by mean.  Otherwise, will be very likely to select categorical variable
    xTemp <- ordinalize(x=allVars[,v], y, ids, sortCat=FALSE, corstr=corstr)
    x <- xTemp$x
    #If all x values the same, do not check optimal split
    if (abs(max(x) - min(x)) > 1e-8) {
      cutpts <- findCutpts(x, minbucket)
      if (!is.null(cutpts)) {
        randCutpt <- ifelse(length(cutpts) == 1, cutpts, sample(cutpts,1))
        #Check whether random cutpoint satisfies quality control
        cutpt <- pvalue_check(y=y, x=x, ids=ids, cutpt=randCutpt, method=method, minpvalue=minpvalue, corstr=corstr)
        if (!is.null(cutpt)) {
          if (is.factor(allVars[,v])) {return(list(varid=v, cutpt=cutpt, x=x, cutToLvl=xTemp$cutToLvl))
          } else {return(list(varid=v, cutpt=cutpt, x=x, cutToLvl=xTemp$cutToLvl))}
        }
      }
    }
  }
}

### RF, ER, AR, and RMRF: calculate splitting statistic ###
splitrule <- function(y, x, ids, cutpts, method, split, corstr){
  
  if (method=="anova") {
    #MSE: Q_m(T)=1/n_L*sum((y-ybar)^2)
    stat <- vapply(cutpts, function(cutpt){
      splitVar <- (x <= cutpt)
      sum((y[splitVar]-mean(y[splitVar]))^2)+sum((y[!splitVar]-mean(y[!splitVar]))^2)
    }, numeric(1))
  } else {
    split <- match.arg(split, c("gini", "information"))
    class1 <- levels(y)[2]
    
    if (!is.null(ids)) {
      #chisq.test((x <= cutpts[1]),(y==class1))
      stat <- vapply(cutpts, geem_mod, yTemp=(y==class1), x=x, ids=ids, corstr=corstr, splitStat="pvalue", numeric(1))
    } else if (split=="gini") {
      #Gini index: Q_m(T)=p*(1-p)
      stat <- vapply(cutpts, function(cutpt){
        splitVar <- (x <= cutpt)
        pL <- mean(y[splitVar]==class1); pR<-mean(y[!splitVar]==class1)
        #Weight each daughter node by the number of observations (i.e. n_L or n_R)
        sumL <- ifelse(pL %in% c(0,1),0,sum(splitVar)*(pL*(1-pL)))
        sumR <- ifelse(pR %in% c(0,1),0,sum(!splitVar)*(pR*(1-pR)))
        sumL+sumR
      }, numeric(1))
    } else {
      #Information index: Q_m(T)=-p*log(p)-(1-p)*log(1-p)
      stat <- vapply(cutpts,function(cutpt){
        splitVar <- (x <= cutpt)
        pL <- mean(y[splitVar]==class1); pR<-mean(y[!splitVar]==class1)
        #Weight each daughter node by the number of observations (i.e. n_L or n_R)
        sumL <- ifelse(pL %in% c(0,1),0,sum(splitVar)*(-pL*log(pL)-(1-pL)*log(1-pL)))
        sumR <- ifelse(pR %in% c(0,1),0,sum(!splitVar)*(-pR*log(pR)-(1-pR)*log(1-pR)))
        sumL+sumR
      }, numeric(1))
    }
  }
  if (all(is.na(stat))) {return(list(cutoff=NA, stat=NA))
  } else {return(list(cutoff = cutpts[which.min(stat)], stat = min(stat,na.rm=TRUE)))}
}

### SSS: Optimization using t-test as splitting rule ###
obj.ttest <- function(cutoffpts, y, x, a=50, scale.y=TRUE, ...){
  if(is.null(cutoffpts)){return(NA)}
  if(scale.y){y <- scale(y, center = TRUE, scale = TRUE)}  # Standardize Y for stability
  score <- NA; n <- length(y)
  grp <- expit(a*outer(as.vector(x),cutoffpts,"-"))
  n.L <- colSums(grp); n.R <- n-n.L
  y1Sum <- as.vector(t(y)%*%grp)
  ybar1 <- y1Sum/n.L; ybar0 <- (colSums(y)-y1Sum)/n.R
  sp2 <- (colSums(y^2)-n.L*ybar1^2-n.R*ybar0^2)/(n-2)  # Compute pooled S2
  #Due to rounding errors, it's possible sp2 is negative or 0 (i.e. sp2=-1x10^-16).  Suppress warning
  tStat <- suppressWarnings((ybar1-ybar0)/sqrt(sp2*(1/n.L+1/n.R)))
  tStat[is.infinite(tStat)] <- NA
  score <- tStat^2
  return(-score)
}

### SSS: Optimization using MSE as splitting rule ###
obj.MSE <- function(cutoffpts, y, x, a=50,...){
  if(is.null(cutoffpts)){return(NA)}
  SSE <- rep(NA,length(cutoffpts))
  grp <- expit(a*outer(as.vector(x),cutoffpts,"-"))
  n <- length(y)
  n.L <- colSums(grp); n.R<-n-n.L
  ySum.L <- as.vector(t(y)%*%grp); ySum.R <- (sum(y)-ySum.L)
  for(i in seq_along(cutoffpts)){
    SSE[i] <- (y-ySum.L[i]/n.L[i])^2%*%grp[,i]+(y-ySum.R[i]/n.R[i])^2%*%(1-grp[,i])
  }
  return(SSE)
}

### SSS: Optimization using entropy or Gini as splitting rule ###
obj.binary <- function(cutoffpts, y, x, a=50, split=c("entropy", "gini"),...){
  if(is.null(cutoffpts)){return(NA)}
  split <- match.arg(split,c("entropy", "gini"))
  Delta.i <- rep(NA,length(cutoffpts))
  n <- length(y); n1 <- sum(y==1); n0 <- n-n1
  grp <- expit(a*outer(as.vector(x),cutoffpts,"-"))
  n.L <- colSums(grp); n.L1 <- colSums(y*grp)
  n.R <- n-n.L; n.R1 <- n1-n.L1
  calcIn <- rep(TRUE, length(n.L))
  #calcIn <- (round(n.L,8) != round(n.L1,8) & round(n.R,8) != round(n.R1,8))
  if(split=="entropy"){
    Delta.i[calcIn] <- n.L1[calcIn]*log(n.L1[calcIn])+(n.L[calcIn]-n.L1[calcIn])*log(n.L[calcIn]-n.L1[calcIn])+
      n.R1[calcIn]*log(n.R1[calcIn])+(n.R[calcIn]-n.R1[calcIn])*log(n.R[calcIn]-n.R1[calcIn])-
      n.L[calcIn]*log(n.L[calcIn])-n.R[calcIn]*log(n.R[calcIn])
  } else if(split=="gini"){
    Delta.i[calcIn] <- -n.L1[calcIn]*(1-n.L1[calcIn]/n.L[calcIn])-n.R1[calcIn]*(1-n.R1[calcIn]/n.R[calcIn])
  }
  return(-Delta.i)
}

### ALL: Convert factors to numerical value. ###
### If performing exhaustive search (for RF, SSS), convert using ranks ###
### If picking random subset of levels (for ER, AR, RMRF), convert to random values ###
ordinalize <- function(x, y, ids, sortCat=TRUE, corstr){
  if (is.factor(x)) {
    x <- factor(x) #Remove factors not listed
    #One can randomly assign a category a distinct numerical value
    if (!sortCat) {
      cutToLvl <- t(sample.int(length(levels(x))))
      colnames(cutToLvl)=levels(x)
    } else {
      #For binary, sort data by proportion in class 1.  For continuous, sort by means
      if (is.factor(y)) {
        if (!is.null(ids)) {
          options(warn=2) #Warning message can give strange estimates
          #Beta coef for glm should be very similar to geem.  However, geem can give incorrect beta, so use glm
          mod <- try(coef(glm((y==levels(y)[2]) ~ x, family = binomial("logit"))), silent=TRUE)
          options(warn=0)
          if (inherits(mod, "try-error")) {cutToLvl <- prop.table(table(y,x),2)[1,,drop=FALSE]
          } else {
            cutToLvl <- t(c(mod[1], mod[1] + mod[2:length(mod)]))
            colnames(cutToLvl) <- levels(x)
          }
        } else {cutToLvl <- prop.table(table(y,x),2)[1,,drop=FALSE]}
      } else {cutToLvl <- t(vapply(levels(x), function(z){mean(y[x==z])}, numeric(1)))}
    }
    #Convert lvls to numerical value. Slow method. Make this faster later.
    xTemp <- rep(NA,length(x))
    for (lvls in levels(x)) {
      xTemp[x==lvls] <- cutToLvl[colnames(cutToLvl)==lvls]
    }
  } else {
    xTemp <- x
    cutToLvl <- NULL
  }
  return(list(x=xTemp, cutToLvl=cutToLvl))
}

### AR and RMRF: Partition data when randomly selecting one variable then one cutpoint method ###
partitionAR <- function(allVars, y, ids, subset, search, method, split, minsplit, minbucket, minpvalue, corstr){
  if (search != "ar") {stop("search not AR?")}
  if (sum(subset) < minsplit) {return(NULL)}
  allVars <- allVars[subset,,drop=FALSE]
  y <- y[subset]
  ids <- ids[subset]
  
  if (NROW(allVars) < 2*minbucket) {stop("Can't split tree to satisfy minbucket")}
  #If all y values are the same, rpart will give an error message.  Do not split tree
  if (length(unique(y)) == 1) {return(NULL)}
  
  randPt <- pickCutpt(allVars, y, ids, method, minsplit, minbucket, nGuesses=50, minpvalue=minpvalue, corstr=corstr)
  #It is possible (unlikely) no cutpoint selected or none can satisfy threshold criterion
  stat <- NA; cutoff <- NA
  if (!is.null(randPt)) {
    mod <- splitrule(y=y, x=randPt$x, ids=ids, cutpts=randPt$cutpt, method=method, split=split, corstr=corstr)
    stat <- mod$stat
    if(is.factor(allVars[,randPt$varid])){
      breakLeft <- rep(NA, length(levels(allVars[,randPt$varid])))
      breakLeft[levels(allVars[,randPt$varid]) %in% colnames(randPt$cutToLvl)[randPt$cutToLvl <= mod$cutoff]]=1L
      breakLeft[levels(allVars[,randPt$varid]) %in% colnames(randPt$cutToLvl)[randPt$cutToLvl > mod$cutoff]]=2L
      if(all(is.na(breakLeft)) & length(unique(breakLeft))<=1){stop("Did not find correct cutpoints")}
    } else {cutoff <- mod$cutoff}
  }
  
  #If each candidate variable cannot be split (e.g. cannot satisfy minbucket), return null
  if (is.na(stat)) {return(NULL)}
  if (is.na(cutoff)) {
    #Index is used for categorical variable splits
    return(partysplit(varid=randPt$varid,
                      index=breakLeft,
                      info=list(stats=stat)))
  } else {
    #Breaks is used for continuous variable splits
    return(partysplit(varid=randPt$varid,
                      breaks=cutoff,
                      info=list(stats=stat)))
  }
}


### RF, ER, and SSS: Partition data when randomly selecting variables then considering all or nsplit candidate splits ###
partition <- function(vars, y, ids, subset, search, method, split, nsplit, nsplit.random,
                      minsplit, minbucket, a, scale.y, useSearch, useOptim, useRpart, allVars, corstr){
  
  if (sum(subset) < minsplit) {return(NULL)}
  
  vars <- vars[subset,,drop=FALSE]
  y <- y[subset]
  ids <- ids[subset]
  
  if (search=="sss") {
    #For binary response, convert response to 0-1
    if(split %in% c("entropy", "gini")){
      stopifnot(length(unique(y)) <= 2)
      yTemp <- (y==levels(y)[2])
    } else {yTemp=y}
  }
  
  if (NROW(vars) < 2*minbucket) {stop("Can't split tree to satisfy minbucket")}
  #If all y values are the same, rpart will give an error message.  Do not split tree
  if (length(unique(y))==1) {return(NULL)}
  
  nVars <- NCOL(vars)
  stats <- rep(NA, nVars)
  cutoff <- rep(NA, nVars)
  breakLeft <- vector("list",nVars)
  
  for (v in 1:nVars) {
    #If randomly picking a subset of categories, do not sort by mean.  Would be more likely to select variables when sorted
    if (search=="exhaustive" && !is.null(nsplit) && nsplit.random) {xTemp <- ordinalize(x=vars[,v], y, ids, sortCat=FALSE, corstr=corstr)
    } else {xTemp <- ordinalize(x=vars[,v], y, ids, sortCat=TRUE, corstr=corstr)}
    x <- xTemp$x
    #If all x values the same, do not check optimal split
    if (abs(max(x) - min(x)) > 1e-8) {
      #The SSS partition deals with problems when there is a very small number of observations
      #Use exhaustive search in this case (or set minsplit >= 5)
      if (search=="sss") {
        obj<-switch(split, ttest=obj.ttest, MSE=obj.MSE, obj.binary)
        
        # Standarize x to apply a constant 'a'
        if (!is.factor(vars[,v])) {sigma <- sd(x); mu <- mean(x); x <- scale(x)}
        cutpts <- findCutpts(x, minbucket)
        if (!is.null(cutpts)) {
          #Go from first possible cutpoint to last possible cutpoint
          pts <- seq(cutpts[1], cutpts[length(cutpts)], length.out=100)
          
          #Decide whether to find global minimum by plugging in several possible values of x, doing optimize, or both
          #plot(pts, opt)
          opt <- Inf; optRev <- data.frame(objective=Inf)
          if (!useSearch && !useOptim) {stop("Need method to search for min")}
          if (useSearch || length(unique(pts)) == 1) {opt <- obj(cutoffpts=pts, y=yTemp, x=x, a=a, scale.y=scale.y, split=split)}
          #Refine minimum by using optimize function.  If only one cutpoint, do not optimize function
          if (useOptim && length(unique(pts)) > 1) {optRev <- suppressWarnings(optimize(obj, lower=min(pts), upper=max(pts), maximum=FALSE,
                                                                                      y=yTemp, x=x, a=a, scale.y=scale.y, split=split))}
          #Change NA to Inf (maybe change later)
          if (all(is.na(opt))) {opt=Inf}; if(is.na(optRev$objective)){optRev$objective=Inf}
        
          #It's possible minbucket cannot be satisfied
          if (!is.infinite(min( min(opt, na.rm=TRUE), optRev$objective, na.rm=TRUE))) {
            stats[v] <- min( min(opt, na.rm=TRUE), optRev$objective, na.rm=TRUE)
            minMethod <- which.min( c( min(opt, na.rm=TRUE), optRev$objective))
          
            if (is.factor(vars[,v])) {
              breakLeft[[v]] <- rep(NA, length(levels(vars[,v])))
              if (minMethod==1) {
                breakLeft[[v]][levels(vars[,v]) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl <= pts[which.min(opt)]]]=1L
                breakLeft[[v]][levels(vars[,v]) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl > pts[which.min(opt)]]]=2L
              } else if (minMethod==2) {
                breakLeft[[v]][levels(vars[,v]) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl <= optRev$minimum]]=1L
                breakLeft[[v]][levels(vars[,v]) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl > optRev$minimum]]=2L
              } else {stop("minMethod unknown")}
              if (all(is.na(breakLeft[[v]])) & length(unique(breakLeft[[v]]))<=1) {stop("Did not find correct cutpoints")}
            } else {
              if (minMethod==1) {cutoff[v] <- pts[which.min(opt)]
              } else if (minMethod==2) {cutoff[v] <- optRev$minimum
              } else {stop("minMethod unknown")}
            }
            if (!is.factor(vars[,v])) {cutoff[v] <- cutoff[v]*sigma + mu}	# Transform back
          }
        }
      } else if (search=="exhaustive") {
        #Note: Rpart uses impurity gain
        if (!useRpart || !is.null(nsplit)) {
          cutpts <- findCutpts(x, minbucket)
          #Take nsplit cutpoints (if applicable)
          if (!is.null(nsplit) && !is.null(cutpts) && length(cutpts) > 1) {
            #If nsplit.random is TRUE, take nsplit cutpts randomly.  Otherwise, take nsplit cutpts equally spread out across cutpts
            if (!nsplit.random & length(cutpts) > nsplit) {
              cutpts <- unique(cutpts[seq(1, length(cutpts), length.out=nsplit)])
            } else {
              cutpts <- sort(sample(cutpts, min(c(nsplit, length(cutpts))), replace=FALSE))
            }
          }
        
          #It is possible (unlikely) no cutpoint can satisfy minbucket
          if (!is.null(cutpts)) {
            mod <- splitrule(y=y, x=x, ids=ids, cutpts=cutpts, method=method, split=split, corstr=corstr)
            if (!is.na(mod$stat)) {
              stats[v] <- mod$stat
              if (is.factor(vars[,v])) {
                breakLeft[[v]] <- rep(NA, length(levels(vars[,v])))
                breakLeft[[v]][levels(vars[,v]) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl <= mod$cutoff]]=1L
                breakLeft[[v]][levels(vars[,v]) %in% colnames(xTemp$cutToLvl)[xTemp$cutToLvl > mod$cutoff]]=2L
                if (all(is.na(breakLeft[[v]])) & length(unique(breakLeft[[v]]))<=1) {stop("Did not find correct cutpoints")}
              } else {cutoff[v] <- mod$cutoff}
            }
          }
        } else { #Can use rpart to do exhaustive search (can't do nsplit)
          x <- vars[,v]
          #Note: use cp=-1 to enforce exactly one split (otherwise, may only have root node)
          #Use minsplit=2 (already checked minsplit) and maxcompete=1, maxsurrogate=0 to speed up convergence
          mod <- rpart(y~x, method=method, parms=list(split=split),
                       control=rpart.control(minsplit=2, minbucket=minbucket, cp=-1, xval=0, maxdepth=1,
                                             maxcompete=1, maxsurrogate=0))
          #It is possible (very unlikely) that no single split can satisfy minbucket constraint
          if(!is.null(mod$splits)){
            stats[v] <- -mod$splits[3]
            if(is.factor(vars[,v])){
              #csplit - 1 left, 2 not present, 3 right.  Change to 1, NA, 2 
              breakLeft[[v]] <- as.vector(mod$csplit)
              breakLeft[[v]][mod$csplit==2]=NA
              breakLeft[[v]][mod$csplit==3]=2L
            } else {cutoff[v] <- mod$splits[4]}
          }
        }
      } else {stop("Unexpected search")}
    }
  }
  #If each candidate variable cannot be split (e.g. cannot satisfy minbucket), return null
  if (all(is.na(stats))) {return(NULL)}
  if (is.na(cutoff[which.min(stats)])) {
    #Index is used for categorical variable splits
    return(partysplit(varid=as.integer(colnames(vars)[which.min(stats)]),
                      index=breakLeft[[which.min(stats)]],
                      info=list(stats=stats)))
  } else {
    #Breaks is used for continuous variable splits
    return(partysplit(varid=as.integer(colnames(vars)[which.min(stats)]),
                      breaks=cutoff[which.min(stats)],
                      info=list(stats=stats)))
  }
}

### ALL: Grow tree by using partition() function several times in a recursive loop ###
growTemp <- function(id=1L, depth=1L, data, response, idVar, subset, search, method, split,
                     mtry, nsplit, nsplit.random, minsplit, minbucket, maxdepth,
                     a, scale.y, useSearch, useOptim, useRpart, minpvalue, corstr){
  #print(c(depth,id))
  if (depth > maxdepth) {return(partynode(id=id))}
  
  y <- data[[response]]
  ids <- data[[idVar]]

  if (search=="ar") {
    sp <- partitionAR(allVars=data[,!(names(data) %in% c(response,idVar)),drop=FALSE], y=y, ids=ids, subset=subset, search=search,
                        method=method, split=split, minsplit=minsplit, minbucket=minbucket, minpvalue=minpvalue, corstr=corstr)
  } else {
    #Select candidate variables
    varSelected <- sort(sample.int(ncol(data)-1-(idVar != ""), mtry))
    vars <- data[varSelected]
    colnames(vars) <- varSelected #Have columns represent varid
    
    sp <- partition(vars=vars, y=y, ids=ids, subset=subset,
                    search=search, method=method, split=split, nsplit=nsplit, nsplit.random=nsplit.random,
                    minsplit=minsplit, minbucket=minbucket, a=a, scale.y=scale.y,
                    useSearch=useSearch, useOptim=useOptim, useRpart=useRpart, corstr=corstr)
  }
  
  if (is.null(sp)) {return(partynode(id=id))}
  
  # Split the data
  kidids <- kidids_split(sp, data=data)
  depth <- depth + 1
  
  kids <- vector(mode="list", length=max(kidids, na.rm=TRUE))
  for (kidid in seq_along(kids)) {
    s <- subset
    s[kidids != kidid] <- FALSE
    # Node ID
    if (kidid > 1) {myid <- max(nodeids(kids[[kidid-1]]))
    } else {myid <- id}
    # Start recursion on this daugther node
    kids[[kidid]] <- growTemp(id=as.integer(myid+1), depth=depth, data=data, response=response, idVar=idVar,
                              subset=s, search=search, method=method, split=split, mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random,
                              minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth,
                              a=a, scale.y=scale.y, useSearch=useSearch, useOptim=useOptim, useRpart=useRpart,
                              minpvalue=minpvalue, corstr=corstr)
  }
  return(partynode(id=as.integer(id), split=sp, kids=kids,
                   info=list(stats=min(info_split(sp)$stats, na.rm=TRUE))))
}

### ALL: Grows a tree with formula interface. Converts to constparty object ###
growTree <- function(formula, data, subset=NULL, search=c("exhaustive","sss","ar"),
                     method=c("anova","class"),
                     split=c("MSE", "ttest", "gini", "entropy", "information"),
                     mtry=NULL, nsplit=NULL, nsplit.random=TRUE, minsplit=20, minbucket=round(minsplit/3), maxdepth=30,
                     a=50, scale.y=TRUE, useSearch=TRUE, useOptim=TRUE, useRpart=FALSE, minpvalue=0.05, corstr="ar1"){
  search <- match.arg(search,c("exhaustive","sss","ar"))
  method <- match.arg(method,c("anova","class"))
  split <- match.arg(split,c("MSE", "ttest", "gini", "entropy", "information"))
  if((method=="anova" && split %in% c("gini", "entropy", "information")) || 
     (method=="class" && split %in% c("MSE", "ttest"))){stop("Split not compatible with method")}
  stopifnot(is.logical(nsplit.random), is.logical(scale.y), is.logical(useSearch), is.logical(useOptim), is.logical(useRpart))
  if (is.numeric(nsplit) && !nsplit.random && nsplit < 5) {"Selecting <5 ordered splits may yield unexpected results"}
  
  response <- all.vars(formula)[1]
  
  if(grepl("\\|", as.character(formula)[3])){
    idVar <- trimws(strsplit(as.character(formula)[3], "\\|")[[1]][2], which="both")
    data <- data[order(data[[idVar]]), ]
  } else {idVar <- ""}
  
  data <- data[c(all.vars(formula)[-1], response)] #Rearrange data so that response comes last
  if (!all(complete.cases(data)) & !is.null(subset)) { paste0("Specifying subset with missing data can yield unexpected results") }
  data <- data[complete.cases(data),] #Complete cases only
  
  if (is.null(mtry)){mtry <- length(all.vars(formula[[3]]))-(idVar != "")}
  
  #if(is.factor(data[[response]])){data[[response]]=as.numeric(data[[response]]==levels(data[[response]])[1])}
  
  if (is.null(subset)){subset <- rep(TRUE, nrow(data))}
  
  # Grow tree
  nodes <- growTemp(id=1L, depth=1L, data=data, response=response, idVar=idVar, subset=subset, search=search, method=method, split=split,
                    mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket,
                    maxdepth=maxdepth, a=a, scale.y=scale.y, useSearch=useSearch, useOptim=useOptim, useRpart=useRpart,
                    minpvalue=minpvalue, corstr=corstr)
  
  # Compute terminal node number for each observation
  fitted <- fitted_node(nodes, data=data)
  # Return rich constparty object
  ret <- party(nodes, data = data,
               fitted = data.frame("(fitted)" = fitted,
                                   "(response)" = data[[response]],
                                   "(weights)" = subset+0,
                                   check.names = FALSE),
               terms = terms(formula))
  as.constparty(ret)
}

### ALL: Grows a random forest ###
growRF <- function(ntrees, formula, data, subset=NULL, search=c("exhaustive","sss","ar"),
                   method=c("anova","class"),
                   split=c("MSE", "ttest", "gini", "entropy", "information"),
                   mtry=NULL, nsplit=NULL, nsplit.random=TRUE, minsplit=20, minbucket=round(minsplit/3), maxdepth=30,
                   a=50, sampleMethod=c('bootstrap','subsample','subsampleByID','learning'),
                   scale.y=TRUE, useSearch=TRUE, useOptim=TRUE, useRpart=FALSE, minpvalue=0.05, corstr="ar1"){
  sampleMethod <- match.arg(sampleMethod, c('bootstrap','subsample','subsampleByID','learning'))

  if(grepl("\\|", as.character(formula)[3])){
    idVar <- trimws(strsplit(as.character(formula)[3], "\\|")[[1]][2], which="both")
  } else {idVar <- ""}
  
  #Construct random forest
  randFor <- lapply(1:ntrees,function(b){
    if(b%%100==0){print(paste0("Tree Number: ",b))}
    obs.b <- switch(sampleMethod,
                    bootstrap = sample.int(nrow(data), size=nrow(data), replace=TRUE),
                    subsample = sample.int(nrow(data), size=round(nrow(data)*0.632), replace=FALSE),
                    subsampleByID = {nIds <- length(unique(data[[idVar]]))
                    unlist(lapply(sample(unique(data[[idVar]]), size=round(nIds*0.632), replace=FALSE),
                                  function(x){which(data[[idVar]] == x)}))},
                    learning = 1:nrow(data))
    sample.b <- data[obs.b,]
    tree.b <- growTree(formula=formula, data=sample.b, subset=subset, search=search, method=method, split=split,
                       mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth, a=a,
                       scale.y=scale.y, useSearch=useSearch, useOptim=useOptim, useRpart=useRpart, minpvalue=minpvalue, corstr=corstr)
    list(tree=tree.b,cases=sort(unique(obs.b)))
  })
  return(randFor)
}

### ALL: Return predicted values from a tree. ###
### Note: users may not want to use predict.party function.  Predict.party assumes indepence (invalid with RMRF). ###
predictTree <- function(tree, newdata=tree$data, type=c("response","prob","node"), corstr="ar1"){
  
  type <- match.arg(type, c("response","prob","node"))
  
  #Note: predict.party works for returning terminal node
  if (type=="node") {return(predict(tree, newdata, type="node"))}
  
  #Determine if tree is classification or regression
  isclass <- is.factor(fitted(tree)[,2])
  
  #Extract formula from tree
  formulaTree <- formula(tree$terms)
  
  if (grepl("\\|", as.character(formulaTree)[3])) {
    idVar <- trimws(strsplit(as.character(formulaTree)[3], "\\|")[[1]][2], which="both")
  } else {idVar <- ""}

  if (!isclass) {
    if (type=="prob") {stop("Cannot return probability for regression")}
    if (type=="response") {return(predict(tree, newdata=newdata, type="response"))}
  } else {
    #If should handle dependence
    if (idVar != "") {
      y <- fitted(tree)[["(response)"]]
      options(warn=2) #Warning message can give strange estimates
      #Beta coef for glm should be very similar to geem.  However, geem can give incorrect beta, so use glm
      mod_glm <- try(coef(glm((y==levels(y)[2]) ~ as.factor(fitted(tree)[["(fitted)"]]), family = binomial("logit"))), silent=TRUE)
      options(warn=0)
      if (inherits(mod_glm, "try-error")) {predictions <- predict(tree, newdata=newdata, type="prob")[,2]
      } else {
        nodesNewdata <- predict(tree, newdata=newdata, type="node")
        nodesToPred <- 1/(1+exp(-c(mod_glm[1],mod_glm[1]+mod_glm[2:length(mod_glm)])))
        names(nodesToPred) <- sort(unique(as.factor(fitted(tree)[["(fitted)"]])))
        predictions <- nodesToPred[as.character(nodesNewdata)]
      }
    } else {predictions <- predict(tree, newdata=newdata, type="prob")[,2]}
    
    if (type=="prob") {return(predictions)}
    if (type=="response") {
      lvls <- levels(fitted(tree)[,2])
      return(factor(lvls[(predictions >= 0.5) + 1], levels=lvls))}
  }
}

### ALL: Return predicted values from a random forest ###
### Note: For binary response, return probability of new obs in 2nd class using mean prediction (not majority vote) ###
predictRF <- function(rf, newdata, prediction=c("overall","by iter"), checkCases=FALSE, corstr="ar1"){
  
  #Return prediction using all trees ("overall") or using first i trees ("by iter")
  prediction <- match.arg(prediction, c("overall","by iter"))
  
  #Cumulative mean with NAs
  cumMeanNA <- function(x){
    xTemp<-x;
    xTemp[is.na(xTemp)] <- 0
    cumsum(xTemp)/cumsum(!is.na(x))
  }
  
  #Determine if RF is classification or regression
  isclass <- is.factor(fitted(rf[[1]]$tree)[,2])
  
  #Each tree predicts prob of new obs in 2nd class and then avg prob across trees
  predictMat <- matrix(NA, ncol=length(rf), nrow=NROW(newdata))
  for (i in seq_along(rf)) {
    #If binary, take probability of class 2
    if (isclass) {
      if (!checkCases) {predictMat[,i] <- predictTree(rf[[i]]$tree, newdata=newdata, type="prob", corstr=corstr)
      } else {predictMat[-rf[[i]]$cases,i] <- predictTree(rf[[i]]$tree, newdata=newdata[-rf[[i]]$cases,], type="prob", corstr=corstr)}
    } else {
      if (!checkCases) {predictMat[,i] <- predictTree(rf[[i]]$tree,newdata=newdata,type="response", corstr=corstr)
      } else {predictMat[-rf[[i]]$cases,i] <- predictTree(rf[[i]]$tree,newdata=newdata[-rf[[i]]$cases,],type="response", corstr=corstr)}
    }
  }
  
  if (prediction=="overall") {return(rowMeans(predictMat,na.rm=TRUE))  #Take mean of predictions
    # return(apply(predictMat,1,median,na.rm=TRUE))  #Take median of predictions (similar to majority vote)
    #Let each col represent prediction based on number of trees used (i.e. col=5 means 5 trees used)
  } else {return(t(apply(predictMat, 1, cumMeanNA)))}
}

### ALL: Minimal depth ###
### WARNING - takes too long for large trees ###
rfMinDepth <- function(rf){
  dataTemp <- names(data_party(rf[[1]]$tree))
  vars <- dataTemp[1:(length(dataTemp)-4)]
  
  mindepth <- rep(0, length(vars))
  for (t in seq_along(rf)) {
    intNodes <- nodeids(rf[[t]]$tree)[-nodeids(rf[[t]]$tree, terminal = TRUE)]
    varsInTree <- vars[unique(unlist(nodeapply(rf[[t]]$tree, ids=intNodes, FUN=function(n){split_node(n)$varid})))]
    varsAtNode <- unlist(nodeapply(rf[[t]]$tree, ids=intNodes, FUN=function(n){split_node(n)$varid}))
    #Root node should be counted as 0
    depthAtNode <- table(unlist(lapply(intNodes, function(x) intersect(intNodes, nodeids(rf[[t]]$tree, from=x)))))-1
    treeDepth <- depth(rf[[t]]$tree)
    
    for (j in seq_along(vars)) {
      if (is.element(vars[j], varsInTree)) { #If variable is in tree
        mindepth[j]=mindepth[j]+min(depthAtNode[varsAtNode==j])
      } else {  #If variable not in tree, set mindepth to maximum depth+1
        mindepth[j]=mindepth[j]+treeDepth+1
      }
    }
  }
  mindepth <- mindepth/length(rf)
  names(mindepth) <- vars
  return(mindepth)
}

### ALL: Calculate variable importance (VIMP) for a given tree and test dataset ###
vimpTree <- function(tree, test, corstr, vimpStat){

  isclass <- is.factor(fitted(tree)[,2])
  
  #Extract formula from tree
  formulaTree <- formula(tree$terms)
  response <- all.vars(formulaTree)[1]
  
  if (grepl("\\|", as.character(formulaTree)[3])) {
    idVar <- trimws(strsplit(as.character(formulaTree)[3], "\\|")[[1]][2], which="both")
  } else {idVar <- ""}
  
  if (isclass) {
    if (vimpStat=="misclass") {
      pred.test <- predictTree(tree, newdata=test, type="response", corstr=corstr)
      err.test <- mean(test[[response]] != pred.test)
    } else if (vimpStat=="auc") {
      pred.test <- predictTree(tree, newdata=test, type="prob", corstr=corstr)
      #Note: higher auc is better, so take area over the curve (or negative sign, will take difference)
      err.test <- 1-auc(test[[response]], pred.test, direction="<")}
  } else {
    pred.test <- predictTree(tree, newdata=test, type="response", corstr=corstr)
    nobs<-length(test[[response]])
    err.test <- 1/nobs*sum((test[[response]]-pred.test)^2)
  }
  
  vars <- names(data_party(tree))[1:(ncol(data_party(tree)) - 4 - (idVar != ""))]
  intNodes <- nodeids(tree)[-nodeids(tree, terminal = TRUE)]
  varsInTree <- vars[unique(unlist(nodeapply(tree, ids=intNodes, FUN=function(n){split_node(n)$varid}), use.names=FALSE))]
  
  permVI <- rep(NA, length(vars))
  for (j in seq_along(vars)) {
    if (is.element(vars[j], varsInTree)) { #if variable is in tree
      permute.test <- test
      permute.test[[vars[j]]] <- permute.test[[vars[j]]][sample.int(nrow(test), nrow(test), replace=FALSE)]
      if (isclass) {
        if (vimpStat=="misclass") {
          pred.permute <- predictTree(tree, newdata=permute.test, type="response", corstr=corstr)
          err.permute <- mean(permute.test[[response]] != pred.permute)
        } else if (vimpStat=="auc") {
          pred.permute <- predictTree(tree, newdata=permute.test, type="prob", corstr=corstr)
          err.permute <- 1-auc(permute.test[[response]], pred.permute, direction="<")}
      } else {
        pred.permute <- predictTree(tree, newdata=permute.test, type="response", corstr=corstr)
        err.permute <- 1/nobs*sum((permute.test[[response]] - pred.permute)^2)}
    } else {err.permute <- err.test}
    permVI[j] <- err.permute - err.test #Take difference in prediction error
  }
  names(permVI) <- vars
  return(list(permVI=permVI))
}

### ALL: Calculate VIMP for a given RF and test dataset ###
vimpRF <- function(rf, test, checkCases=FALSE, corstr="ar1", vimpStat=c("misclass", "auc")){
  vimpStat <- match.arg(vimpStat, c("misclass", "auc"))
  for (t in seq_along(rf)) {
    #if(t%%50==0){print(paste0("Tree Number: ",t))}
    if (checkCases == FALSE) {treeVI <- vimpTree(rf[[t]]$tree, test=test, corstr=corstr, vimpStat=vimpStat)    #VI using test sample
    } else {treeVI <- vimpTree(tree=rf[[t]]$tree, test=test[-rf[[t]]$cases,], corstr=corstr, vimpStat=vimpStat)}    #VI using OOB sample
    if (any(is.na(treeVI$permVI))) {stop(paste0("Variable has NA VI, check tree: ",t))}
    #Instead of adding rows to matrix, preset dimensions of matrix
    if (t == 1) {
      permVI <- matrix(NA, nrow=length(rf), ncol=length(treeVI$permVI))
      colnames(permVI) <- names(treeVI$permVI)
    }
    permVI[t,] <- treeVI$permVI
  }
  #unscaledvimp: VI = xbar; scaledvimp: Z = xbar/(sigma/sqrt(n))
  return(list(unscaledVI = colMeans(permVI),
              scaledVI = colMeans(permVI)/(apply(permVI, 2, sd)/sqrt(length(rf)))))
}

### ALL: Plots VIMP ###
vimpPlot <- function(vimp, title="") {
  vimpSeq<-seq_along(vimp)
  par(mar=c(8, 4, 1, 2) + 0.1)
  plot(range(vimpSeq)+c(-.15,.15), range(c(0,vimp)),  type="n", yaxt="s", xaxt="n", xlab="", ylab="Importance")
  mtext("Predictor", 1, line=7)
  title(title)
  for (j in vimpSeq) {
    polygon(rep(c(vimpSeq[j]-0.15, vimpSeq[j]+0.15), rep(2,2)), c(0, vimp[j], vimp[j], 0), col="lightgray")
  }
  axis(1, at=vimpSeq, labels = names(vimp), tick=FALSE, las=2)
}


### ALL: Grows tree setup for parallel computing ###
growRF_Parallel_temp <- function(b, formula, data, search, method,
                                 split, mtry, nsplit, nsplit.random,
                                 minsplit, minbucket, maxdepth,
                                 a, sampleMethod,
                                 scale.y, useSearch, useOptim, useRpart, minpvalue, corstr) {
  library(rpart)
  library(partykit)
  library(geeM)
  
  if (b %% 100 == 0) {print(paste0("Tree Number: ",b))}
  
  if (grepl("\\|", as.character(formula)[3])) {
    idVar <- trimws(strsplit(as.character(formula)[3], "\\|")[[1]][2], which="both")
  } else {idVar <- ""}
  
  obs.b <- switch(sampleMethod,
                  bootstrap = sample.int(nrow(data), size=nrow(data), replace=TRUE),
                  subsample = sample.int(nrow(data), size=round(nrow(data)*0.632), replace=FALSE),
                  subsampleByID = {nIds <- length(unique(data[[idVar]]))
                                   unlist(lapply(sample(unique(data[[idVar]]), size=round(nIds*0.632), replace=FALSE),
                                                 function(x){which(data[[idVar]] == x)}))},
                  learning = 1:nrow(data))
  sample.b <- data[obs.b,]
  tree.b <- growTree(formula=formula, data=sample.b, search=search, method=method, split=split,
                     mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random, minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth, a=a,
                     scale.y=scale.y, useSearch=useSearch, useOptim=useOptim, useRpart=useRpart, minpvalue=minpvalue, corstr=corstr)
  list(tree=tree.b,cases=sort(unique(obs.b)))
}


### ALL: Grows a random forest using parallel computing ###
### Note: not as many checks performed.  Recommend using growRF() function first and get code to work, ###
### then switch to growRF_Parallel() function ###
growRF_Parallel <- function(ntrees, formula, data, search, method,
                             split=c("MSE", "ttest", "gini", "entropy", "information"),
                             mtry, nsplit=NULL, nsplit.random=TRUE, minsplit, minbucket, maxdepth,
                             a=50, sampleMethod=c('bootstrap','subsample','subsampleByID','learning'),
                             scale.y=TRUE, useSearch=TRUE, useOptim=TRUE, useRpart=FALSE, minpvalue=0.05, corstr="ar1",
                             iseed=111111){
  sampleMethod <- match.arg(sampleMethod, c('bootstrap','subsample','subsampleByID','learning'))
  
  #Construct random forest
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl,c("expit","findCutpts","geem_mod","pvalue_check","pickCutpt","splitrule",
                     "obj.ttest","obj.MSE","obj.binary","ordinalize",
                     "partitionAR","partition","growTemp","growTree","growRF_Parallel_temp"))
  clusterSetRNGStream(cl=cl, iseed=iseed)
  randFor <- parLapply(cl,1:ntrees, growRF_Parallel_temp,
                       formula=formula, data=data, search=search, method=method,
                       split=split, mtry=mtry, nsplit=nsplit, nsplit.random=nsplit.random,
                       minsplit=minsplit, minbucket=minbucket, maxdepth=maxdepth,
                       a=a, sampleMethod=sampleMethod,
                       scale.y=scale.y, useSearch=useSearch, useOptim=useOptim, useRpart=useRpart, minpvalue=minpvalue, corstr=corstr)
  stopCluster(cl)
  return(randFor)
}

