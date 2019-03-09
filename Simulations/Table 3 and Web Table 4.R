#rm(list=ls(all=TRUE))

library(dplyr)

# Table 3 has both night and patient-level factors
sim_rf <- read.csv(paste0("Data/RF Comparison/sim_rf_nightAndPt_betaLog2.csv"), header=TRUE)
sim_glm <- read.csv(paste0("Data/RF Comparison/sim_glm_nightAndPt_betaLog2.csv"), header=TRUE)

# Web Table 4 only considers night-level factors.  Can load following data instead to reproduce Web Table 4
#sim_rf <- read.csv(paste0("Data/RF Comparison/sim_rf_night_betaLog2.csv"), header=TRUE)
#sim_glm <- read.csv(paste0("Data/RF Comparison/sim_glm_night_betaLog2.csv"), header=TRUE)

sim <- cbind(glmIndep_vimp=sim_glm$glmIndepAllVars_vimp, glmDep_vimp=sim_glm$glmDepAllVars_vimp, sim_rf)

metrics <- paste0(c("arDep", "arIndep", "erDep", "erIndep", "rfDep", "rfIndep", "glmDep", "glmIndep"),"_vimp")

#Function records the rank of X1 (the informative variable) and transforms dataset to vertical)
vi_max1 <- function(var){setNames(aggregate(list(sim[[var]]), list(sim$model), FUN=function(x){mean(floor(x)==1)}), c("Model", var))}
datVert <- data.frame(lapply(metrics, vi_max1))
datIndep <- datVert[c("arIndep_vimp", "erIndep_vimp", "rfIndep_vimp", "glmIndep_vimp")]
colnames(datIndep) <- c("RMRF", "ER", "RF", "GLM")
datDep <- data.frame(datVert[c("arDep_vimp", "erDep_vimp", "rfDep_vimp", "glmDep_vimp")])
colnames(datDep) <- c("RMRF", "ER", "RF", "GLM")
table <- data.frame(model = rep(LETTERS[1:5], 2), method = rep(c("Robust","Naive"), each=5))
table <- cbind(table, rbind(datDep, datIndep))
arrange(table, model)




