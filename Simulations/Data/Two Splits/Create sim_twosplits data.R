#rm(list=ls(all=TRUE))

library(parallel)

singleSim <- function(x, method=c("random effects", "marginal"),
                      intercept, betas, sigma, rho){

  library(SimCorMultRes)
  source('../../../R Functions/RF functions 20JAN19.R')
  
  method <- match.arg(method, c("random effects", "marginal"))

  N <- 100   #Number of patients
  clsize <- 4   #Number of nights
  #Note: x1, x3 are patient-level factors.  x2 and x4 and night-level factors
  x1 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  x3 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  x2 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  x4 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  x1_Ind <- (x1 <= 0.3)
  x2_Ind <- (x2 <= 0.3)
  
  if (method == "random effects") {
    w <- rep(rnorm(N, mean=0, sd=sigma), each=clsize)
    prob <- 1/(1+exp(-(intercept+betas[1]*x1_Ind+betas[2]*x2_Ind+w)))
    dat <- data.frame(id=rep(seq.int(N), each=clsize), time=rep(seq.int(clsize), N), x1, x2, x3, x4,
                      Class=as.factor(rbinom(clsize*N, size=1, prob=prob)))
  }
  if (method == "marginal") {
    cor.matrix <- toeplitz(c(1, rep(rho,clsize-1)))  #Correlation matrix
    dat <- rbin(clsize = clsize, intercepts = intercept, beta = betas, xformula = ~ x1_Ind + x2_Ind,
                cor.matrix = cor.matrix, link = "logit")$simdata
    dat$x1 <- x1; dat$x2 <- x2; dat$x3 <- x3; dat$x4 <- x4
    dat$Class <- as.factor(dat$y)
    dat$y <- dat$x1_Ind <- dat$x2_Ind <- NULL
  }
  
  treeIndep <- growTree(formula=Class ~ x1 + x2 + x3 + x4,
                        data=dat, search="exhaustive", method="class", split="gini",
                        minsplit=3, minbucket=1, maxdepth=2)
  
  treeDep <- growTree(formula=Class ~ x1 + x2 + x3 + x4 | id,
                      data=dat, search="exhaustive", method="class", split="gini", corstr="exchangeable",
                      minsplit=3, minbucket=1, maxdepth=2)
  
  #Need to convert id value into factor to observe 0 counts
  treeIndep_vimp <- table(factor(unlist(nodeapply(treeIndep, ids=nodeids(treeIndep)[-nodeids(treeIndep, terminal = TRUE)],
                                                           FUN=function(n){split_node(n)$varid})), levels=1:4))
  treeDep_vimp <- table(factor(unlist(nodeapply(treeDep, ids=nodeids(treeDep)[-nodeids(treeDep, terminal = TRUE)],
                                              FUN=function(n){split_node(n)$varid})), levels=1:4))
  
  if (all(treeIndep_vimp==0) || all(treeDep_vimp==0)) {stop("No splits selected? ")}

  
  # Create test dataset for accuracy #
  N <- 1000
  clsize <- 4
  x1 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  x3 <- rep(sample((1:10)/10, size=N, replace=TRUE), each=clsize)
  x2 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  x4 <- sample((1:10)/10, size=clsize*N, replace=TRUE)
  x1_Ind <- (x1 <= 0.3)
  x2_Ind <- (x2 <= 0.3)
  if (method == "random effects") {
    w <- rep(rnorm(N, mean=0, sd=sigma), each=clsize)
    prob <- 1/(1+exp(-(intercept+betas[1]*x1_Ind+betas[2]*x2_Ind+w)))
    test <- data.frame(id=rep(seq.int(N), each=clsize), time=rep(seq.int(clsize), N), x1, x2, x3, x4,
                       Class=as.factor(rbinom(clsize*N, size=1, prob=prob)))
  }
  if (method == "marginal") {
    cor.matrix <- toeplitz(c(1, rep(rho,clsize-1)))  #Correlation matrix
    test <- rbin(clsize = clsize, intercepts = intercept, beta = betas, xformula = ~ x1_Ind + x2_Ind,
                 cor.matrix = cor.matrix, link = "logit")$simdata
    test$x1 <- x1; test$x2 <- x2; test$x3 <- x3; test$x4 <- x4
    test$Class <- as.factor(test$y)
    test$y <- test$x1_Ind <- test$x2_Ind <- NULL
  }
  pred_Indep <- predictTree(treeIndep, test, type="prob")
  pred_Dep <- predictTree(treeDep, test, type="prob", corstr="exchangeable")
  lvls <- levels(test$Class)
  treeIndep_acc <- mean(lvls[1+(pred_Indep > 0.5)] == test$Class)
  treeDep_acc <- mean(lvls[1+(pred_Dep > 0.5)] == test$Class)
  treeIndep_auc <- auc(test$Class, pred_Indep, direction="<")
  treeDep_auc <- auc(test$Class, pred_Dep, direction="<")
  
  return(data.frame(treeIndep_inc=all(treeIndep_vimp[1:2] > 0), treeIndep_exc=all(treeIndep_vimp[1:2] > 0 & treeIndep_vimp[3:4] == 0),
              treeDep_inc=all(treeDep_vimp[1:2] > 0), treeDep_exc=all(treeDep_vimp[1:2] > 0 & treeDep_vimp[3:4] == 0),
              treeIndep_acc=treeIndep_acc, treeDep_acc=treeDep_acc,
              treeIndep_auc=treeIndep_auc, treeDep_auc=treeDep_auc))
}

#singleSim(x=1, method="marginal", intercept=0, betas=c(log(2), log(2)), rho=0, sigma=NA)

simPar<-function(nsim, method, intercept, betas, rho, sigma){
  
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl,"singleSim")
  clusterSetRNGStream(cl=cl, iseed=431897)
  simuls <- parLapply(cl, 1:nsim, singleSim, method=method,
                      intercept=intercept, betas=betas, rho=rho, sigma=sigma)
  stopCluster(cl)
  
  return(do.call("rbind", simuls))
}


#nsim <- 5000  #Manuscript uses 5,000 simulations on supercomputers
nsim <- 50

##Manuscript considers betas equal to log(2) or log(4)
#betas <- c(log(2), log(4)); outputText <- "beta1Log2_beta2Log4"
#betas <- c(log(4), log(2)); outputText <- "beta1Log4_beta2Log2"
betas <- c(log(2),log(2)); outputText <- "beta1Log2_beta2Log2"

sim_modA <- simPar(nsim=nsim, method="marginal", intercept=0, betas=betas, rho=0, sigma=NA)
sim_modB <- simPar(nsim=nsim, method="random effects", intercept=0, betas=betas, rho=NA, sigma=1)
sim_modC <- simPar(nsim=nsim, method="random effects", intercept=0, betas=betas, rho=NA, sigma=3)
sim_modD <- simPar(nsim=nsim, method="marginal", intercept=0, betas=betas, rho=0.4, sigma=NA)
sim_modE <- simPar(nsim=nsim, method="marginal", intercept=0, betas=betas, rho=0.9, sigma=NA)

sim_twosplits <- data.frame(rbind(sim_modA, sim_modB, sim_modC, sim_modD, sim_modE),
                            model=rep(LETTERS[1:5], each=nsim), beta1=betas[1], beta2=betas[2])

#write.table(sim_twosplits, file = paste0("sim_twosplits_", outputText, ".csv"), sep = ",", row.names=FALSE)


