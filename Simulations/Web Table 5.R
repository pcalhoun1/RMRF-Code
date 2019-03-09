#rm(list=ls(all=TRUE))

library(dplyr)

simA <- read.csv(paste0("Data/GLM vs RF/sim_glmVsRF_modA_betaLog2sigma3.csv"))
simB <- read.csv(paste0("Data/GLM vs RF/sim_glmVsRF_modB_betaLog2sigma3.csv"))
simC <- read.csv(paste0("Data/GLM vs RF/sim_glmVsRF_modC_betaLog6sigma3.csv"))
simD <- read.csv(paste0("Data/GLM vs RF/sim_glmVsRF_modD_betaLog6sigma3.csv"))

sim <- rbind(simA, simB, simC, simD)



metrics <- c("rmrfDep", "rmrfIndep", "erDep", "erIndep", "rfDep", "rfIndep", "glmDep", "glmIndep")

#Record whether both X1 and X2 are both in the top 2 ranks
vi_top2 <- function(var, model) { setNames(mean(apply(apply(sim[sim$model == model, gsub("_vimp[12345]", "", names(sim)) %in% var], 2, function(x){x %in% c(1,2)}), 1, all)), var)}
datVert <- NULL
for (mod in LETTERS[1:4]) {
  datVert <- rbind(datVert, unlist(lapply(metrics, vi_top2, model=mod)))
}
datVert <- as.data.frame(datVert*100)

datIndep <- datVert[c("rmrfIndep", "erIndep", "rfIndep", "glmIndep")]
colnames(datIndep) <- c("RMRF", "ER", "RF", "GLM")
datDep <- data.frame(datVert[c("rmrfDep", "erDep", "rfDep", "glmDep")])
colnames(datDep) <- c("RMRF", "ER", "RF", "GLM")
table <- data.frame(model = rep(c("Linear", "Quadratic", "Linear Interaction", "Quadratic Interaction"), 2),
                    method = rep(c("Robust","Naive"), each=4))
table <- cbind(table, rbind(datDep, datIndep))
arrange(table, model)


# Can manually check calculations too:
#Dependent
mean(simA$rmrfDep_vimp1 %in% c(1,2) & simA$rmrfDep_vimp2 %in% c(1,2))
mean(simA$erDep_vimp1 %in% c(1,2) & simA$erDep_vimp2 %in% c(1,2))
mean(simA$rfDep_vimp1 %in% c(1,2) & simA$rfDep_vimp2 %in% c(1,2))
mean(simA$glmDep_vimp1 %in% c(1,2) & simA$glmDep_vimp2 %in% c(1,2))


#Independent
mean(simA$rmrfIndep_vimp1 %in% c(1,2) & simA$rmrfIndep_vimp2 %in% c(1,2))
mean(simA$erIndep_vimp1 %in% c(1,2) & simA$erIndep_vimp2 %in% c(1,2))
mean(simA$rfIndep_vimp1 %in% c(1,2) & simA$rfIndep_vimp2 %in% c(1,2))
mean(simA$glmIndep_vimp1 %in% c(1,2) & simA$glmIndep_vimp2 %in% c(1,2))


