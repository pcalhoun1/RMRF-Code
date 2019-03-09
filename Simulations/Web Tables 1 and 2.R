#rm(list=ls(all=TRUE))

library(plyr)

bias <- function(value, cutpt) mean(value-cutpt)
sdCutpt <- function(value, cutpt) sd(value)
rmse <- function(value, cutpt) sqrt(mean((value-cutpt)^2))

# Web Table 1 sets beta = log(2)
sim_cutpt5 <- read.csv(paste0("Data/Cutpoint/sim_cutpoint_cutpt5_betaLog2.csv"), header=TRUE)
sim_cutpt3 <- read.csv(paste0("Data/Cutpoint/sim_cutpoint_cutpt3_betaLog2.csv"), header=TRUE)

# Web Table 2 sets beta = log(4).  Can load following data instead to show those results
#sim_cutpt5 <- read.csv(paste0("Data/Cutpoint/sim_cutpoint_cutpt5_betaLog4.csv"), header=TRUE)
#sim_cutpt3 <- read.csv(paste0("Data/Cutpoint/sim_cutpoint_cutpt3_betaLog4.csv"), header=TRUE)

sim_cutpt5Vert <- data.frame(cutpt=c(sim_cutpt5$depCut, sim_cutpt5$indepCut),
                            method=rep(c("1.Robust","2.Naive"), each=nrow(sim_cutpt5)),
                            model=rep(sim_cutpt5$model,2), cutpoint=rep(sim_cutpt5$cutpoint,2),
                            beta=rep(sim_cutpt5$beta,2), N=rep(sim_cutpt5$N,2), clsize=rep(sim_cutpt5$clsize,2))

sim_cutpt3Vert <- data.frame(cutpt=c(sim_cutpt3$depCut, sim_cutpt3$indepCut),
                            method=rep(c("1.Robust","2.Naive"), each=nrow(sim_cutpt3)),
                            model=rep(sim_cutpt3$model,2), cutpoint=rep(sim_cutpt3$cutpoint,2),
                            beta=rep(sim_cutpt3$beta,2), N=rep(sim_cutpt3$N,2), clsize=rep(sim_cutpt3$clsize,2))


#Produce Web Table:
ddply(sim_cutpt5Vert, ~ model + cutpoint + method, summarise, Bias=round(mean(cutpt - cutpoint),3), SD=round(sd(cutpt),3), RMSE=round(sqrt(mean((cutpt-cutpoint)^2)),3))
ddply(sim_cutpt3Vert, ~ model + cutpoint + method, summarise, Bias=round(mean(cutpt - cutpoint),3), SD=round(sd(cutpt),3), RMSE=round(sqrt(mean((cutpt-cutpoint)^2)),3))


