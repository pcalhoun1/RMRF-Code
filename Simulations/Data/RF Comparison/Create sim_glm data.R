#rm(list=ls(all=TRUE))

library(parallel)

singleSim_GLM <- function(x, method=c("random effects", "marginal"),
                         intercept, var, beta, sigma, rho){
  
  library(SimCorMultRes)
  source('../../../R Functions/RF functions 20JAN19.R')
  
  method <- match.arg(method, c("random effects", "marginal"))
  
  N <- 100   #Number of patients
  clsize <- 4   #Number of nights
  #Check if patient-level or night-level
  if (var == "night and patient") {
    x1 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x2 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x3 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x4 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
    x5 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  } else if (var == "only night") {
    x1 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x2 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x3 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x4 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x5 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  } else if (var == "only patient") {
    x1 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
    x2 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
    x3 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
    x4 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
    x5 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  } else {print(paste0("Unexpected var? ", var))}
  x1_Ind <- (x1 <= 0.5)
  
  if (method=="random effects") {
    w <- rep(rnorm(N, mean=0, sd=sigma), each=clsize)
    prob <- 1/(1+exp(-(intercept+beta*x1_Ind+w)))
    dat <- data.frame(id=rep(seq.int(N), each=clsize), time=rep(seq.int(clsize), N), x1, x2, x3, x4, x5,
                      Class=as.factor(rbinom(clsize*N, size=1, prob=prob)))
  } else if (method=="marginal") {
    cor.matrix <- toeplitz(c(1, rep(rho,clsize-1)))  #Correlation matrix
    dat <- rbin(clsize = clsize, intercepts = intercept, beta = beta, xformula = ~ x1_Ind,
                cor.matrix = cor.matrix, link = "logit")$simdata
    dat$x1 <- x1; dat$x2 <- x2; dat$x3 <- x3; dat$x4 <- x4; dat$x5 <- x5
    dat$Class <- as.factor(dat$y)
    dat$y <- dat$x1_Ind <- NULL
  }
  
  # Standardize variables to obtain variable importance
  for (x in 1:5) {
    dat[[paste0("x",x,"_standard")]] <- (dat[[paste0("x",x)]] - mean(dat[[paste0("x",x)]])) / sd(dat[[paste0("x",x)]])
  }
  
  # Naive Logistic Regression using All Variables
  form <- (dat$Class==1) ~ x1_standard + x2_standard + x3_standard + x4_standard + x5_standard 
  glmIndepAllVars <- coef(glm(form, data=dat, family = binomial("logit")))
  glmIndepAllVars_vimp <- rank(-abs(glmIndepAllVars)[-1])
  
  # Robust Logistic Regression using All Variables
  form <- (dat$Class==1) ~ x1_standard + x2_standard + x3_standard + x4_standard + x5_standard 
  glmDepAllVars <- coef(geem(form, data=dat, family = binomial("logit"), corstr='exchangeable', id=dat$id))
  glmDepAllVars_vimp <- rank(-abs(glmDepAllVars)[-1])
  
  # Naive Logistic Regression using one variable at a time
  glmIndepSingleVars <- rep(NA, 5)
  for (x in 1:5) {
    form <- (dat$Class==1) ~ dat[[paste0("x",x,"_standard")]]
    glmIndepSingleVarsTemp <- coef(glm(form, data=dat, family = binomial("logit")))
    glmIndepSingleVars[x] <- glmIndepSingleVarsTemp[2]
  }
  names(glmIndepSingleVars) <- paste0("x",1:5,"_standard")
  glmIndepSingleVars_vimp <- rank(-abs(glmIndepSingleVars))
  
  # Robust Logistic Regression using one variable at a time
  glmDepSingleVars <- rep(NA, 5)
  for (x in 1:5) {
    form <- (dat$Class==1) ~ dat[[paste0("x",x,"_standard")]]
    glmDepSingleVarsTemp <- coef(geem(form, data=dat, family = binomial("logit"), corstr='exchangeable', id=dat$id))
    glmDepSingleVars[x] <- glmDepSingleVarsTemp[2]
  }
  names(glmDepSingleVars) <- paste0("x",1:5,"_standard")
  glmDepSingleVars_vimp <- rank(-abs(glmDepSingleVars))

  return(list(glmIndepAllVars_vimp=glmIndepAllVars_vimp, glmIndepSingleVars_vimp=glmIndepSingleVars_vimp,
              glmDepAllVars_vimp=glmDepAllVars_vimp, glmDepSingleVars_vimp=glmDepSingleVars_vimp))
}

#singleSim_GLM(1, method="marginal", intercept=0, var="night and patient", beta=log(2), rho=0, sigma=NA)


simPar_GLM <- function(nsim, method, intercept, var, beta, rho, sigma){
  
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl,"singleSim_GLM")
  clusterSetRNGStream(cl=cl, iseed=431897)
  simuls <- parLapply(cl, 1:nsim, singleSim_GLM, method=method,
                      intercept=intercept, var=var, beta=beta, rho=rho, sigma=sigma)
  stopCluster(cl)
  
  glmIndepAllVars_vimp <- matrix(NA, ncol=5, nrow=nsim)
  glmIndepSingleVars_vimp <- matrix(NA, ncol=5, nrow=nsim)
  glmDepAllVars_vimp <- matrix(NA, ncol=5, nrow=nsim)
  glmDepSingleVars_vimp <- matrix(NA, ncol=5, nrow=nsim)
  for(sim in 1:nsim){
    glmIndepAllVars_vimp[sim,] <- simuls[[sim]]$glmIndepAllVars_vimp
    glmIndepSingleVars_vimp[sim,] <- simuls[[sim]]$glmIndepSingleVars_vimp
    glmDepAllVars_vimp[sim,] <- simuls[[sim]]$glmDepAllVars_vimp
    glmDepSingleVars_vimp[sim,] <- simuls[[sim]]$glmDepSingleVars_vimp
  }
  return(list(glmIndepAllVars_vimp=glmIndepAllVars_vimp, glmIndepSingleVars_vimp=glmIndepSingleVars_vimp,
              glmDepAllVars_vimp=glmDepAllVars_vimp, glmDepSingleVars_vimp=glmDepSingleVars_vimp))
}


#nsim <- 500  #Manuscript uses 500 simulations.  Takes around 10min
nsim <- 2
#var <- "only night"; outputText <- "night_betaLog2"  # Web Table 4 has only night-level factors
var <- "night and patient"; outputText <- "nightAndPt_betaLog2"  # Table 3 has both night and patient level factors
beta <- log(2)

sim_modA <- simPar_GLM(nsim=nsim, method="marginal", intercept=0, var=var, beta=beta, rho=0, sigma=NA)
sim_modB <- simPar_GLM(nsim=nsim, method="random effects", intercept=0, var=var, beta=beta, rho=NA, sigma=1)
sim_modC <- simPar_GLM(nsim=nsim, method="random effects", intercept=0, var=var, beta=beta, rho=NA, sigma=3)
sim_modD <- simPar_GLM(nsim=nsim, method="marginal", intercept=0, var=var, beta=beta, rho=0.4, sigma=NA)
sim_modE <- simPar_GLM(nsim=nsim, method="marginal", intercept=0, var=var, beta=beta, rho=0.9, sigma=NA)

# Function takes glm results and returns rank of X1 for each simulation
vi_max1 <- function(simData){
  data.frame(lapply(names(simData), function(x){setNames(data.frame(simData[[x]][,1]), x)}))
}

sim_GLM <- data.frame(rbind(vi_max1(sim_modA), vi_max1(sim_modB), vi_max1(sim_modC), vi_max1(sim_modD), vi_max1(sim_modE)),
                     model=rep(LETTERS[1:5], each=nsim), var=var, beta=beta)
#write.table(sim_GLM, file = paste0("sim_GLM_", outputText, ".csv"), sep = ",", row.names=FALSE)


mean(vi_max1(sim_modA)$glmIndepSingleVars_vimp == 1)
mean(vi_max1(sim_modA)$glmIndepSingleVars_vimp == 1)
mean(vi_max1(sim_modA)$glmDepAllVars_vimp == 1)
mean(vi_max1(sim_modA)$glmDepSingleVars_vimp == 1)


