#rm(list=ls(all=TRUE))

library(parallel)

singleSim_RFvsGLM <- function(x, model=LETTERS[1:4], intercept, beta, sigma, rho){
  
  library(SimCorMultRes)
  source('C:/Users/pcalh/Desktop/Research/RF code/RF functions 02OCT18.R')
  #source('RF functions 02OCT18.R')
  
  model <- match.arg(model, LETTERS[1:4])
  
  N <- 800   #Number of patients
  clsize <- 4   #Number of nights
  x1 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  x2 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  x3 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  x4 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  x5 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  x1_sine <- sin(2*pi*x1)
  x4_sine <- sin(2*pi*x4)

  w <- rep(rnorm(N, mean=0, sd=sigma), each=clsize)
  if (model=="A") {prob <- 1/(1+exp(-(intercept + beta*x1 + beta*x4 + w)))
  } else if (model=="B") {prob <- 1/(1+exp(-(intercept + beta*x1_sine + beta*x4_sine + w)))
  } else if (model=="C") {prob <- 1/(1+exp(-(intercept + beta*x1*x4 + w)))
  } else if (model=="D") {prob <- 1/(1+exp(-(intercept + beta*x1_sine*x4_sine + w)))
  } else { stop(paste0("Unexpected model",model)) }
  dat <- data.frame(id=rep(seq.int(N), each=clsize), time=rep(seq.int(clsize), N), x1, x2, x3, x4, x5,
                    Class=as.factor(rbinom(clsize*N, size=1, prob=prob)))

  print(paste0("Simulation: ", x))
  
  # Standardize variables to obtain variable importance
  for (x in 1:5) {
    dat[[paste0("x",x,"_standard")]] <- (dat[[paste0("x",x)]] - mean(dat[[paste0("x",x)]])) / sd(dat[[paste0("x",x)]])
  }
  
  # Naive Logistic Regression
  form <- (dat$Class==1) ~ x1_standard + x2_standard + x3_standard + x4_standard + x5_standard 
  glmIndep <- coef(glm(form, data=dat, family = binomial("logit")))
  glmIndep_vimp <- rank(-abs(glmIndep)[-1])
  
  # Robust Logistic Regression
  form <- (dat$Class==1) ~ x1_standard + x2_standard + x3_standard + x4_standard + x5_standard 
  glmDep <- coef(geem(form, data=dat, family = binomial("logit"), corstr='exchangeable', id=dat$id))
  glmDep_vimp <- rank(-abs(glmDep)[-1])
  
  # Build RFs
  rfIndep <- growRF(ntrees=500, formula=Class ~ x1 + x2 + x3 + x4 + x5,
                    data=dat, search="exhaustive", method="class", split="gini",
                    mtry=2, nsplit=NULL, minsplit=3, minbucket=1, maxdepth=2, sampleMethod="subsample")
  
  rfDep <- growRF(ntrees=500, formula=Class ~ x1 + x2 + x3 + x4 + x5 | id,
                  data=dat, search="exhaustive", method="class", split="gini",
                  mtry=2, nsplit=NULL, minsplit=3, minbucket=1, maxdepth=2, sampleMethod="subsampleByID", corstr="exchangeable")
  
  erIndep <- growRF(ntrees=500, formula=Class ~ x1 + x2 + x3 + x4 + x5,
                    data=dat, search="exhaustive", method="class", split="gini",
                    mtry=2, nsplit=1, minsplit=3, minbucket=1, maxdepth=2, sampleMethod="subsample")
  
  erDep <- growRF(ntrees=500, formula=Class ~ x1 + x2 + x3 + x4 + x5 | id,
                  data=dat, search="exhaustive", method="class", split="gini",
                  mtry=2, nsplit=1, minsplit=3, minbucket=1, maxdepth=2, sampleMethod="subsampleByID", corstr="exchangeable")
  
  rmrfIndep <- growRF(ntrees=500, formula=Class ~ x1 + x2 + x3 + x4 + x5,
                      data=dat, search="ar", method="class", split="gini",
                      mtry=2, nsplit=1, minsplit=3, minbucket=1, maxdepth=2, sampleMethod="subsample")
  
  rmrfDep <- growRF(ntrees=500, formula=Class ~ x1 + x2 + x3 + x4 + x5 | id,
                    data=dat, search="ar", method="class", split="gini",
                    mtry=2, nsplit=1, minsplit=3, minbucket=1, maxdepth=2, sampleMethod="subsampleByID", corstr="exchangeable")
  
  rfIndep_vimp <- rank(-vimpRF(rfIndep, dat, checkCases=TRUE, vimpStat="auc")$scaledVI)
  rfDep_vimp <- rank(-vimpRF(rfDep, dat, checkCases=TRUE, vimpStat="auc")$scaledVI)
  erIndep_vimp <- rank(-vimpRF(erIndep, dat, checkCases=TRUE, vimpStat="auc")$scaledVI)
  erDep_vimp <- rank(-vimpRF(erDep, dat, checkCases=TRUE, vimpStat="auc")$scaledVI)
  rmrfIndep_vimp <- rank(-vimpRF(rmrfIndep, dat, checkCases=TRUE, vimpStat="auc")$scaledVI)
  rmrfDep_vimp <- rank(-vimpRF(rmrfDep, dat, checkCases=TRUE, vimpStat="auc")$scaledVI)

  names(rfIndep_vimp) <- names(rfDep_vimp) <- names(erIndep_vimp) <- names(erDep_vimp) <- names(rmrfIndep_vimp) <- names(rmrfDep_vimp) <- paste0("x",1:5)
  
  return(list(glmIndep_vimp=glmIndep_vimp, glmDep_vimp=glmDep_vimp, rfIndep_vimp=rfIndep_vimp, rfDep_vimp=rfDep_vimp, erIndep_vimp=erIndep_vimp, erDep_vimp=erDep_vimp, rmrfIndep_vimp=rmrfIndep_vimp, rmrfDep_vimp=rmrfDep_vimp))
}

#singleSim_RFvsGLM(1, model="C", intercept=0, beta=log(2), rho=NA, sigma=3)


simPar_RFvsGLM <- function(nsim, model, intercept, beta, rho, sigma){
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl,"singleSim_RFvsGLM")
  clusterSetRNGStream(cl=cl, iseed=431897)
  simuls <- parLapply(cl, 1:nsim, singleSim_RFvsGLM, model=model,
                      intercept=intercept, beta=beta, rho=rho, sigma=sigma)
  stopCluster(cl)
  
  glmIndep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  glmDep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  rfIndep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  rfDep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  erIndep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  erDep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  rmrfIndep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  rmrfDep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  for(sim in 1:nsim){
    glmIndep_vimp[sim,] <- simuls[[sim]]$glmIndep_vimp
    glmDep_vimp[sim,] <- simuls[[sim]]$glmDep_vimp
    rfIndep_vimp[sim,] <- simuls[[sim]]$rfIndep_vimp
    rfDep_vimp[sim,] <- simuls[[sim]]$rfDep_vimp
    erIndep_vimp[sim,] <- simuls[[sim]]$erIndep_vimp
    erDep_vimp[sim,] <- simuls[[sim]]$erDep_vimp
    rmrfIndep_vimp[sim,] <- simuls[[sim]]$rmrfIndep_vimp
    rmrfDep_vimp[sim,] <- simuls[[sim]]$rmrfDep_vimp
  }
  return(list(glmIndep_vimp=glmIndep_vimp, glmDep_vimp=glmDep_vimp,
              rfIndep_vimp=rfIndep_vimp, rfDep_vimp=rfDep_vimp,
              erIndep_vimp=erIndep_vimp, erDep_vimp=erDep_vimp, 
              rmrfIndep_vimp=rmrfIndep_vimp, rmrfDep_vimp=rmrfDep_vimp))
}

#nsim <- 100  # Manuscript uses 100 simulations
nsim <- 2
modLetter <- "A"
beta <- log(2); sigma=3; outputText <- paste0("mod",modLetter,"_betaLog2sigma3")

sim_mod <- simPar_RFvsGLM(nsim=nsim, model=modLetter, intercept=0, beta=beta, rho=NA, sigma=sigma)

# Function takes results and returns rank of X1 and X4
vi_max <- function(simData){
  data.frame(lapply(names(simData), function(x){setNames(data.frame(simData[[x]][,c(1,4)]), paste0(x, c("1", "4")))}))
}

sim_RFvsGLM <- data.frame(vi_max(sim_mod), model=modLetter, beta=beta, sigma=sigma)

write.table(sim_RFvsGLM, file = paste0("glmVsRF_", outputText, ".csv"), sep = ",", row.names=FALSE)


