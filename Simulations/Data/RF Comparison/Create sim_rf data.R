#rm(list=ls(all=TRUE))

library(parallel)

singleSim_RF <- function(x, method=c("random effects", "marginal"),
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
  
  # BUILD RF
  rfIndep <- growRF(ntrees=1000, formula=Class ~ x1 + x2 + x3 + x4 + x5,
                    data=dat, search="exhaustive", method="class", split="gini",
                    mtry=2, nsplit=NULL, minsplit=3, minbucket=1, maxdepth=1, sampleMethod="subsample")
  
  rfDep <- growRF(ntrees=1000, formula=Class ~ x1 + x2 + x3 + x4 + x5 | id,
                  data=dat, search="exhaustive", method="class", split="gini",
                  mtry=2, nsplit=NULL, minsplit=3, minbucket=1, maxdepth=1, sampleMethod="subsampleByID", corstr="exchangeable")
  
  erIndep <- growRF(ntrees=1000, formula=Class ~ x1 + x2 + x3 + x4 + x5,
                    data=dat, search="exhaustive", method="class", split="gini",
                    mtry=2, nsplit=1, minsplit=3, minbucket=1, maxdepth=1, sampleMethod="subsample")
  
  erDep <- growRF(ntrees=1000, formula=Class ~ x1 + x2 + x3 + x4 + x5 | id,
                  data=dat, search="exhaustive", method="class", split="gini",
                  mtry=2, nsplit=1, minsplit=3, minbucket=1, maxdepth=1, sampleMethod="subsampleByID", corstr="exchangeable")
  
  rmrfIndep <- growRF(ntrees=1000, formula=Class ~ x1 + x2 + x3 + x4 + x5,
                      data=dat, search="ar", method="class", split="gini",
                      mtry=1, nsplit=1, minsplit=3, minbucket=1, maxdepth=1, sampleMethod="subsample")
  
  rmrfDep <- growRF(ntrees=1000, formula=Class ~ x1 + x2 + x3 + x4 + x5 | id,
                    data=dat, search="ar", method="class", split="gini",
                    mtry=1, nsplit=1, minsplit=3, minbucket=1, maxdepth=1, sampleMethod="subsampleByID", corstr="exchangeable")

  rfIndep_vimp <- prop.table(table(factor(unlist(lapply(rfIndep, function(x){nodeapply(x$tree, ids=1, FUN=function(n){split_node(n)$varid})})), levels=1:5)))
  rfDep_vimp <- prop.table(table(factor(unlist(lapply(rfDep, function(x){nodeapply(x$tree, ids=1, FUN=function(n){split_node(n)$varid})})), levels=1:5)))
  erIndep_vimp <- prop.table(table(factor(unlist(lapply(erIndep, function(x){nodeapply(x$tree, ids=1, FUN=function(n){split_node(n)$varid})})), levels=1:5)))
  erDep_vimp <- prop.table(table(factor(unlist(lapply(erDep, function(x){nodeapply(x$tree, ids=1, FUN=function(n){split_node(n)$varid})})), levels=1:5)))
  rmrfIndep_vimp <- prop.table(table(factor(unlist(lapply(rmrfIndep, function(x){nodeapply(x$tree, ids=1, FUN=function(n){split_node(n)$varid})})), levels=1:5)))
  rmrfDep_vimp <- prop.table(table(factor(unlist(lapply(rmrfDep, function(x){nodeapply(x$tree, ids=1, FUN=function(n){split_node(n)$varid})})), levels=1:5)))
  
  names(rfIndep_vimp) <- names(rfDep_vimp) <- names(erIndep_vimp) <- names(erDep_vimp) <- names(rmrfIndep_vimp) <- names(rmrfDep_vimp) <- paste0("x",1:5)
  
  return(list(rfIndep_vimp=rfIndep_vimp, rfDep_vimp=rfDep_vimp, erIndep_vimp=erIndep_vimp, erDep_vimp=erDep_vimp, rmrfIndep_vimp=rmrfIndep_vimp, rmrfDep_vimp=rmrfDep_vimp))
}

#singleSim_RF(1, method="marginal", intercept=0, var="night and patient", beta=log(2), rho=0, sigma=NA)




simPar_RF<-function(nsim, method, intercept, var, beta, rho, sigma){
  
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl,"singleSim_RF")
  clusterSetRNGStream(cl=cl, iseed=431897)
  simuls <- parLapply(cl, 1:nsim, singleSim_RF, method=method,
                      intercept=intercept, var=var, beta=beta, rho=rho, sigma=sigma)
  stopCluster(cl)
  
  rfIndep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  rfDep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  erIndep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  erDep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  rmrfIndep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  rmrfDep_vimp <- matrix(NA, ncol=5, nrow=nsim)
  for(sim in 1:nsim){
    rfIndep_vimp[sim,] <- simuls[[sim]]$rfIndep_vimp
    rfDep_vimp[sim,] <- simuls[[sim]]$rfDep_vimp
    erIndep_vimp[sim,] <- simuls[[sim]]$erIndep_vimp
    erDep_vimp[sim,] <- simuls[[sim]]$erDep_vimp
    rmrfIndep_vimp[sim,] <- simuls[[sim]]$rmrfIndep_vimp
    rmrfDep_vimp[sim,] <- simuls[[sim]]$rmrfDep_vimp
  }
  return(list(rfIndep_vimp=rfIndep_vimp, rfDep_vimp=rfDep_vimp,
              erIndep_vimp=erIndep_vimp, erDep_vimp=erDep_vimp, 
              rmrfIndep_vimp=rmrfIndep_vimp, rmrfDep_vimp=rmrfDep_vimp))
}

#nsim <- 500  #Manuscript uses 500 simulations on supercomputers
nsim <- 2
#var <- "only night"; outputText <- "night_betaLog2"  # Web Table 4 has only night-level factors
var <- "night and patient"; outputText <- "nightAndPt_betaLog2"  # Table 3 has both night and patient level factors
beta <- log(2)

sim_modA <- simPar_RF(nsim=nsim, method="marginal", intercept=0, var=var, beta=beta, rho=0, sigma=NA)
sim_modB <- simPar_RF(nsim=nsim, method="random effects", intercept=0, var=var, beta=beta, rho=NA, sigma=1)
sim_modC <- simPar_RF(nsim=nsim, method="random effects", intercept=0, var=var, beta=beta, rho=NA, sigma=3)
sim_modD <- simPar_RF(nsim=nsim, method="marginal", intercept=0, var=var, beta=beta, rho=0.4, sigma=NA)
sim_modE <- simPar_RF(nsim=nsim, method="marginal", intercept=0, var=var, beta=beta, rho=0.9, sigma=NA)

#Function takes random forest and returns whether X1 has highest VI
vi_max1 <- function(simData){
  data.frame(lapply(names(simData), function(x){setNames(data.frame(apply(-simData[[x]], 1, rank)[1,]), x)}))
}

sim_rf <- data.frame(rbind(vi_max1(sim_modA), vi_max1(sim_modB), vi_max1(sim_modC), vi_max1(sim_modD), vi_max1(sim_modE)),
                     model=rep(LETTERS[1:5], each=nsim), var=var, beta=beta)
#write.table(sim_rf, file = paste0("sim_rf_", outputText, ".csv"), sep = ",", row.names=FALSE)



