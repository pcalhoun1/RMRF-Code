#rm(list=ls(all=TRUE))

library(parallel)

singleSim <- function(x, method=c("random effects", "marginal"),
                      var, cutpoint,
                      intercept, beta, sigma, rho, N, clsize){
  
  library(SimCorMultRes)
  source('../../../R Functions/RF functions 20JAN19.R')
  
  method <- match.arg(method, c("random effects", "marginal"))
  
  if(var=="patient-level"){x1 <- rep(sample((1:50)/50, size=N, replace=TRUE), each=clsize)
  } else if (var=="night-level"){x1 <- sample((1:50)/50, size=clsize*N, replace=TRUE)
  } else {stop("Unexpected var")}
  x1_Ind <- (x1 <= cutpoint)
  
  if(method=="random effects"){
    w <- rep(rnorm(N, mean=0, sd=sigma), each=clsize)
    prob <- 1/(1+exp(-(intercept+beta*x1_Ind+w)))
    dat <- data.frame(id=rep(seq.int(N), each=clsize), time=rep(seq.int(clsize), N), x1,
                      Class=as.factor(rbinom(clsize*N, size=1, prob=prob)))
  }
  if(method=="marginal"){
    cor.matrix <- toeplitz(c(1, rep(rho,clsize-1)))  #Correlation matrix
    dat <- rbin(clsize = clsize, intercepts = intercept, beta = beta, xformula = ~ x1_Ind,
                cor.matrix = cor.matrix, link = "logit")$simdata
    dat$x1 <- x1; dat$Class <- as.factor(dat$y)
    dat$y <- dat$x1_Ind <- NULL
  }
  
  treeIndep <- growTree(formula=Class ~ x1,
                    data=dat, search="exhaustive", method="class", split="gini",
                    minsplit=3, minbucket=1, maxdepth=1)
  
  treeDep <- growTree(formula=Class ~ x1 | id,
                  data=dat, search="exhaustive", method="class", split="gini", corstr="exchangeable",
                  minsplit=3, minbucket=1, maxdepth=1)
  
  indepCut <- as.numeric(unlist(nodeapply(treeIndep, ids=1, FUN=function(n){split_node(n)$breaks})))
  depCut <- as.numeric(unlist(nodeapply(treeDep, ids=1, FUN=function(n){split_node(n)$breaks})))
  
  return(data.frame(indepCut=indepCut, depCut=depCut))
}

#singleSim(x=1, method="marginal", var="night-level", cutpoint=0.5, intercept=0, beta=log(2), rho=0, sigma=NA, N=100, clsize=4)


simPar<-function(nsim, method, var, cutpoint, intercept, beta, rho, sigma, N, clsize){
  
  num_cores <- detectCores()-1
  cl <- makeCluster(num_cores)
  clusterExport(cl,"singleSim")
  clusterSetRNGStream(cl=cl, iseed=431897)
  simuls <- parLapply(cl, 1:nsim, singleSim, method=method, var=var, cutpoint=cutpoint,
                      intercept=intercept, beta=beta, rho=rho, sigma=sigma, N=N, clsize=clsize)
  stopCluster(cl)
  
  return(do.call("rbind", simuls))
}


#nsim <- 5000  #Manuscript uses 5,000 simulations on supercomputers
nsim <- 50
N <- 100; clsize <- 4

#Manuscript considers cutpoints at 0.5 and 0.3, and beta = log(2) or log(4)
#cutpoint <- 0.3; beta <- log(2); outputText = "cutpt3_betaLog2"
#cutpoint <- 0.3; beta <- log(4); outputText = "cutpt3_betaLog4"
#cutpoint <- 0.5; beta <- log(4); outputText = "cutpt5_betaLog4"
cutpoint <- 0.5; beta <- log(2); outputText = "cutpt5_betaLog2"

sim_modA_night <- simPar(nsim=nsim, method="marginal", var="night-level", cutpoint=cutpoint,
                         intercept=0, beta=beta, rho=0, sigma=NA, N=N, clsize=clsize)
sim_modB_night <- simPar(nsim=nsim, method="random effects", var="night-level", cutpoint=cutpoint,
                         intercept=0, beta=beta, rho=NA, sigma=1, N=N, clsize=clsize)
sim_modC_night <- simPar(nsim=nsim, method="random effects", var="night-level", cutpoint=cutpoint,
                         intercept=0, beta=beta, rho=NA, sigma=3, N=N, clsize=clsize)
sim_modD_night <- simPar(nsim=nsim, method="marginal", var="night-level", cutpoint=cutpoint,
                         intercept=0, beta=beta, rho=0.4, sigma=NA, N=N, clsize=clsize)
sim_modE_night <- simPar(nsim=nsim, method="marginal", var="night-level", cutpoint=cutpoint,
                         intercept=0, beta=beta, rho=0.9, sigma=NA, N=N, clsize=clsize)

sim_cutpt <- data.frame(rbind(sim_modA_night, sim_modB_night, sim_modC_night, sim_modD_night, sim_modE_night),
                        model=rep(LETTERS[1:5], each=nsim), cutpoint=cutpoint, beta=beta, N=N, clsize=clsize)
#write.table(sim_cutpt, file = paste0("sim_cutpoint_",outputText,".csv"), sep = ",", row.names=FALSE)

