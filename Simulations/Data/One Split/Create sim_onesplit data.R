#rm(list=ls(all=TRUE))

library(parallel)

### Look at patient vs night biases
singleSim <- function(x, method=c("random effects", "marginal"),
                      var, cutpoint,
                      intercept, beta, sigma, rho){
  
  library(SimCorMultRes)
  source('../../../R Functions/RF functions 20JAN19.R')
  
  method <- match.arg(method, c("random effects", "marginal"))
  
  N <- 100   #Number of patients
  clsize <- 4   #Number of nights
  #Specify whether x1 and x2 are patient or night-level (x3 and x4 be converse)
  if (var=="patient-level") {
    x1 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
    x2 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
    x3 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x4 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  } else if (var=="night-level") {
    x1 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x2 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
    x3 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
    x4 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  } else {stop("Unexpected var")}
  x1_Ind <- (x1 <= cutpoint)
    
  if (method == "random effects") {
    w <- rep(rnorm(N, mean=0, sd=sigma), each=clsize)
    prob <- 1/(1+exp(-(intercept+beta*x1_Ind+w)))
    dat <- data.frame(id=rep(seq.int(N), each=clsize), time=rep(seq.int(clsize), N), x1, x2, x3, x4,
                      Class=as.factor(rbinom(clsize*N, size=1, prob=prob)))
  }
  if (method == "marginal") {
    cor.matrix <- toeplitz(c(1, rep(rho,clsize-1)))  #Correlation matrix
    dat <- rbin(clsize = clsize, intercepts = intercept, beta = beta, xformula = ~ x1_Ind,
                cor.matrix = cor.matrix, link = "logit")$simdata
    dat$x1 <- x1; dat$x2 <- x2; dat$x3 <- x3; dat$x4 <- x4; dat$Class <- as.factor(dat$y)
    dat$y <- dat$x1_Ind <- NULL
  }
  
  treeIndep <- growTree(formula=Class ~ x1 + x2 + x3 + x4,
                        data=dat, search="exhaustive", method="class", split="gini",
                        minsplit=3, minbucket=1, maxdepth=1)
  
  treeDep <- growTree(formula=Class ~ x1 + x2 + x3 + x4 | id,
                      data=dat, search="exhaustive", method="class", split="gini", corstr="exchangeable",
                      minsplit=3, minbucket=1, maxdepth=1)
  
  treeIndep_vimp <- prop.table(table(factor(unlist(nodeapply(treeIndep, ids=1, FUN=function(n){split_node(n)$varid})), levels=1:4)))
  treeDep_vimp <- prop.table(table(factor(unlist(nodeapply(treeDep, ids=1, FUN=function(n){split_node(n)$varid})), levels=1:4)))
  names(treeIndep_vimp) <- names(treeDep_vimp) <- paste0("x",1:4)
  
  return(list(treeIndep_vimp=treeIndep_vimp, treeDep_vimp=treeDep_vimp))
}

#singleSim(x=1, method="marginal", var="night-level", cutpoint=0.5, intercept=0, beta=log(2), rho=0, sigma=NA)


simPar<-function(nsim, method, var, cutpoint, intercept, beta, rho, sigma){
  
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl,"singleSim")
  clusterSetRNGStream(cl=cl, iseed=431897)
  simuls <- parLapply(cl, 1:nsim, singleSim, method=method, var=var, cutpoint=cutpoint,
                      intercept=intercept, beta=beta, rho=rho, sigma=sigma)
  stopCluster(cl)

  treeIndep_vimp <- matrix(NA, ncol=4, nrow=nsim)
  treeDep_vimp <- matrix(NA, ncol=4, nrow=nsim)
  for(sim in 1:nsim){
    treeIndep_vimp[sim,] <- simuls[[sim]]$treeIndep_vimp
    treeDep_vimp[sim,] <- simuls[[sim]]$treeDep_vimp
  }
  return(list(treeIndep_vimp=treeIndep_vimp, treeDep_vimp=treeDep_vimp))
}


#nsim <- 5000  #Manuscript uses 5,000 simulations on supercomputers
nsim <- 50
var <- "night-level"  #Manuscript reports results for informative night-level factor
cutpoint <- 0.3
#beta <- log(4); outputText <- "night_cutpt3_betaLog2"  #Manuscript considers beta equal to log(2) and log(4)
beta <- log(2); outputText <- "night_cutpt3_betaLog2" 


sim_modA <- simPar(nsim=nsim, method="marginal", var=var, cutpoint=cutpoint,
                   intercept=0, beta=beta, rho=0, sigma=NA)
sim_modB <- simPar(nsim=nsim, method="random effects", var=var, cutpoint=cutpoint,
                   intercept=0, beta=beta, rho=NA, sigma=1)
sim_modC <- simPar(nsim=nsim, method="random effects", var=var, cutpoint=cutpoint,
                   intercept=0, beta=beta, rho=NA, sigma=3)
sim_modD <- simPar(nsim=nsim, method="marginal", var=var, cutpoint=cutpoint,
                   intercept=0, beta=beta, rho=0.4, sigma=NA)
sim_modE <- simPar(nsim=nsim, method="marginal", var=var, cutpoint=cutpoint,
                   intercept=0, beta=beta, rho=0.9, sigma=NA)

sim_onesplit <- data.frame(treeIndep=rbind(sim_modA$treeIndep_vimp, sim_modB$treeIndep_vimp, sim_modC$treeIndep_vimp, sim_modD$treeIndep_vimp, sim_modE$treeIndep_vimp),
                           treeDep=rbind(sim_modA$treeDep_vimp, sim_modB$treeDep_vimp, sim_modC$treeDep_vimp, sim_modD$treeDep_vimp, sim_modE$treeDep_vimp),
                           model=rep(LETTERS[1:5], each=nsim), var=var, cutpoint=cutpoint, beta=beta)

#write.table(sim_onesplit, file = paste0("sim_onesplit_", outputText, ".csv"), sep = ",", row.names=FALSE)



