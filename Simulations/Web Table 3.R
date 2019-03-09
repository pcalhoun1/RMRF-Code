#rm(list=ls(all=TRUE))

library(dplyr)

# Web Table 3 sets beta1 = log(2) and beta2 = log(2)
sim_twosplits <- read.csv(paste0("Data/Two Splits/sim_twosplits_beta1Log2_beta2Log2.csv"), header=TRUE)

# Manuscript text also reports results when beta1 = log(2) and beta2 = log(4), and when beta1 = log(4) and beta2 = log(2).
# Can load following data instead to show those results
#sim_twosplits <- read.csv(paste0("Data/Two Splits/sim_twosplits_beta1Log4_beta2Log2.csv"), header=TRUE)
#sim_twosplits <- read.csv(paste0("Data/Two Splits/sim_twosplits_beta1Log2_beta2Log4.csv"), header=TRUE)

metrics <- c("treeDep_inc", "treeDep_exc", "treeIndep_inc", "treeIndep_exc", "treeDep_acc", "treeIndep_acc", "treeDep_auc", "treeIndep_auc")
vi_mean <- function(var){setNames(aggregate(list(sim_twosplits[[var]]), list(sim_twosplits$model), FUN=function(x){round(mean(x),3)}), c("Model", var))}
datVert <- data.frame(lapply(metrics, vi_mean))
datIndep <- datVert[c("treeIndep_inc","treeIndep_exc","treeIndep_acc","treeIndep_auc")]
colnames(datIndep) <- c("Inclusive", "Exclusive", "Average Accuracy", "Average AUC")
datDep <- data.frame(datVert[c("treeDep_inc","treeDep_exc","treeDep_acc","treeDep_auc")])
colnames(datDep) <- c("Inclusive", "Exclusive", "Average Accuracy", "Average AUC")
table <- data.frame(model = rep(LETTERS[1:5], 2), method = rep(c("Robust","Naive"), each=5))
table <- cbind(table, rbind(datDep, datIndep))
arrange(table, model)


