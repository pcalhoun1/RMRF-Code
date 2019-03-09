#rm(list=ls(all=TRUE))

### Diabetes VIMP ###

diabetes_vimp1 <- read.csv("Data/RMRF results/diabetes_vimp_seed1.csv", header=TRUE)
diabetes_vimp2 <- read.csv("Data/RMRF results/diabetes_vimp_seed2.csv", header=TRUE)
diabetes_vimp3 <- read.csv("Data/RMRF results/diabetes_vimp_seed3.csv", header=TRUE)

#Average 3 VI metrics
diabetes_vimp <- cbind((diabetes_vimp1[1:13] + diabetes_vimp2[1:13] + diabetes_vimp3[1:13])/3, diabetes_vimp1[14])

#Null variable importance when specifying dependence structure
nullDep <- read.csv("Data/RMRF results/nullVimp_dep.csv")
#Calculate p-values for each variable using dependence structure
vimp_dep_pvalue <- function(factor) {
  nullDep_auc <- c(unlist(nullDep[nullDep$factor == factor, seq(2, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", factor])
  nullDep_misclass <- c(unlist(nullDep[nullDep$factor == factor, seq(3, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_misclass", factor])
  list(pvalue_auc = mean(nullDep_auc >= diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", factor]),
       pvalue_misclass = mean(nullDep_misclass >= diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_misclass", factor]))
}
pvaluesDep <- sapply(as.character(nullDep$factor), vimp_dep_pvalue)
pvaluesDepAUC <- pvaluesDep[1, ]
pvaluesDepMisclass <- pvaluesDep[2, ]


#Null variable importance when ignoring dependence structure
null_indep <- read.csv("Data/RMRF results/nullVimp_indep.csv")
#Calculate p-values for each variable when ignoring dependence structure
vimp_indep_pvalue <- function(factor) {
  null_auc <- c(unlist(null_indep[null_indep$factor == factor, seq(2, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arIndep_auc", factor])
  null_misclass <- c(unlist(null_indep[null_indep$factor == factor, seq(3, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arIndep_misclass", factor])
  list(pvalue_auc = mean(null_auc >= diabetes_vimp[diabetes_vimp$vimp == "arIndep_auc", factor]),
       pvalue_misclass = mean(null_misclass >= diabetes_vimp[diabetes_vimp$vimp == "arIndep_misclass", factor]))
}
pvaluesIndep <- sapply(as.character(nullDep$factor), vimp_indep_pvalue)
pvaluesIndepAUC <- pvaluesIndep[1, ]
pvaluesIndepMisclass <- pvaluesIndep[2, ]


##### Results in manuscript text #####

names(pvaluesDepAUC)[pvaluesDepAUC < 0.05]  #Variables that are significant (p-value < 0.05) accounting for dependence and using AUC
names(pvaluesDepMisclass)[pvaluesDepMisclass < 0.05]  #Variables that are significant (p-value < 0.05) accounting for dependence and using Misclassification
names(pvaluesIndepAUC)[pvaluesIndepAUC < 0.05]  #Variables that are significant (p-value < 0.05) ignorning for dependence and using AUC


##### Figure 1 #####

vimp <- diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", names(diabetes_vimp) != 'vimp']
vimpSeq <- seq_along(vimp)


#tiff("Figures/Fig1.tiff", width = 5, height = 4, units = "in", pointsize = 12, res=600)
par(mar=c(5.7, 2.5, 0, 0) + 0.1)
plot(NA, xlim=c(1,13), ylim=c(0,35),  type="n", axes=FALSE, xlab="", ylab="")
mtext("Importance", 2, font=2, line=1.6, cex=1.2)
for (j in vimpSeq){
  if(j < 7){polygon(rep(c(vimpSeq[j]-0.25, vimpSeq[j]+0.25), rep(2,2)), c(0, vimp[j], vimp[j], 0), col="gray45", xpd=TRUE)
  } else {polygon(rep(c(vimpSeq[j]-0.25, vimpSeq[j]+0.25), rep(2,2)), c(0, vimp[j], vimp[j], 0), col="gray80", xpd=TRUE)}
}
axis(1, at=vimpSeq, labels = c("Age", "Gender", "HbA1c", "", "", "Total Daily Dose",
                               "Bedtime BG", "", "IOB", "CGM ROC", "Time Activated", "Exercise", "Daytime Hypo"), tick=FALSE, las=2, cex.axis=0.9, line=-1.35)
axis(1, at=vimpSeq, labels = c(rep("",7),"Bedtime Snack", rep("",5)), tick=FALSE, las=2, cex.axis=0.9, line=-0.3)
axis(1, at=vimpSeq, labels = c(rep("",4),"Basal Insulin",rep("",8)), tick=FALSE, las=2, cex.axis=0.9, line=0.2)
axis(1, at=vimpSeq, labels = c(rep("",3),"Diab Duration", rep("",9)), tick=FALSE, las=2, cex.axis=0.9, line=-1.1)
axis(2, at=seq(0,35,5), las=1, font.axis=2, mgp=c(3, .6, 0), xpd=TRUE)

legend(0.7, 36, c("Patient-Level Factor","Night-Level Factor"),pch=22, pt.bg=c("gray45","gray80"), bty='n', pt.cex=1.8, cex=1, adj=c(0,0.5), horiz=TRUE, xpd=TRUE,
       text.width=c(0,5.7))

#Add asterisk for significant variables
text((1:13)[pvaluesDepAUC < 0.05], vimp[pvaluesDepAUC < 0.05]+2, "*", cex=2)

dev.off()



##### Web Figure 2 #####


#tiff("Figures/WebFig2.tiff", width = 6, height = 8, units = "in", pointsize = 12, res=600)
par(mfrow=c(3, 2))

# HbA1c #
factor <- as.character(nullDep$factor)[3]
nullDep_auc <- c(unlist(nullDep[nullDep$factor == factor, seq(2, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", factor])
x <- hist(unlist(nullDep_auc), breaks=20, plot=FALSE)

par(mar=c(2.2, 3, 0, 0.1) + 0.1)
plot(x, xlim=c(-15, 30), ylim=c(0,21), axes=FALSE, xlab="", ylab="", main="")
axis(1, at=seq(-15, 30, 5), tick=TRUE, cex.axis=1.3, line=-0.5, mgp=c(3, .6, 0))
axis(1, at=c(-10, 15, 25), tick=TRUE, cex.axis=1.3, line=-0.5, mgp=c(3, .6, 0))
axis(2, at=seq(0,20,5), las=1, font.axis=2, mgp=c(3, .6, 0), xpd=TRUE)
mtext("Frequency", 2, font=2, line=1.6, cex=1.2)
mtext("Importance for HbA1c", 1, font=2, line=1.2, cex=1.1)

polygon(rep(c(26,28), each=2), c(0,1,1,0), col='darkgray')
legend(21, 15, paste0("p-value = ", round(pvaluesDepAUC[[factor]],2)), pch=NA, xjust=0.5, yjust=0.5, x.intersp=0, cex=1.3, bg='white')


# Bedtime BG #
factor <- as.character(nullDep$factor)[7]
nullDep_auc <- c(unlist(nullDep[nullDep$factor == factor, seq(2, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", factor])
x <- hist(unlist(nullDep_auc), breaks=20, plot=FALSE)

par(mar=c(2.2, 3, 0, 0.1) + 0.1)
plot(x, xlim=c(-10, 25), ylim=c(0,21), axes=FALSE, xlab="", ylab="", main="")
axis(1, at=seq(-10, 25, 5), tick=TRUE, cex.axis=1.3, line=-0.5, mgp=c(3, .6, 0))
axis(2, at=seq(0,20,5), las=1, font.axis=2, mgp=c(3, .6, 0), xpd=TRUE)
mtext("Frequency", 2, font=2, line=1.6, cex=1.2)
mtext("Importance for Bedtime BG", 1, font=2, line=1.2, cex=1.1)

polygon(rep(c(22,24), each=2), c(0,1,1,0), col='darkgray')
legend(18, 15, paste0("p-value = ", round(pvaluesDepAUC[[factor]],2)), pch=NA, xjust=0.5, yjust=0.5, x.intersp=0, cex=1.3, bg='white')


# IOB #
factor <- as.character(nullDep$factor)[9]
nullDep_auc <- c(unlist(nullDep[nullDep$factor == factor, seq(2, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", factor])
x <- hist(unlist(nullDep_auc), breaks=20, plot=FALSE)

par(mar=c(2.2, 3, 0, 0.1) + 0.1)
plot(x, xlim=c(-15, 15), ylim=c(0, 21), axes=FALSE, xlab="", ylab="", main="")
axis(1, at=seq(-15, 15, 5), tick=TRUE, cex.axis=1.3, line=-0.5, mgp=c(3, .6, 0))
axis(2, at=seq(0,20,5), las=1, font.axis=2, mgp=c(3, .6, 0), xpd=TRUE)
mtext("Frequency", 2, font=2, line=1.6, cex=1.2)
mtext("Importance for IOB", 1, font=2, line=1.2, cex=1.1)

polygon(rep(c(12,13), each=2), c(0,1,1,0), col='darkgray')
legend(9, 15, paste0("p-value = ", round(pvaluesDepAUC[[factor]],2)), pch=NA, xjust=0.5, yjust=0.5, x.intersp=0, cex=1.3, bg='white')


# Time Activated #
factor <- as.character(nullDep$factor)[11]
nullDep_auc <- c(unlist(nullDep[nullDep$factor == factor, seq(2, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", factor])
x <- hist(unlist(nullDep_auc), breaks=20, plot=FALSE)

par(mar=c(2.2, 3, 0, 0.1) + 0.1)
plot(x, xlim=c(-15, 15), ylim=c(0,21), axes=FALSE, xlab="", ylab="", main="")
axis(1, at=seq(-15, 15, 5), tick=TRUE, cex.axis=1.3, line=-0.5, mgp=c(3, .6, 0))
axis(2, at=seq(0,20,5), las=1, font.axis=2, mgp=c(3, .6, 0), xpd=TRUE)
mtext("Frequency", 2, font=2, line=1.6, cex=1.2)
mtext("Importance for Time Activated", 1, font=2, line=1.2, cex=1.1)

polygon(rep(c(12,13), each=2), c(0,1,1,0), col='darkgray')
legend(9, 15, paste0("p-value = ", round(pvaluesDepAUC[[factor]],2)), pch=NA, xjust=0.5, yjust=0.5, x.intersp=0, cex=1.3, bg='white')


# Exercise #
factor <- as.character(nullDep$factor)[12]
nullDep_auc <- c(unlist(nullDep[nullDep$factor == factor, seq(2, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", factor])
x <- hist(unlist(nullDep_auc), breaks=20, plot=FALSE)

par(mar=c(2.2, 3, 0, 0.1) + 0.1)
plot(x, xlim=c(-15, 20), ylim=c(0,21), axes=FALSE, xlab="", ylab="", main="")
axis(1, at=seq(-15, 20, 5), tick=TRUE, cex.axis=1.3, line=-0.5, mgp=c(3, .6, 0))
axis(2, at=seq(0,20,5), las=1, font.axis=2, mgp=c(3, .6, 0), xpd=TRUE)
mtext("Frequency", 2, font=2, line=1.6, cex=1.2)
mtext("Importance for Exercise", 1, font=2, line=1.2, cex=1.1)

polygon(rep(c(15,16), each=2), c(0,1,1,0), col='darkgray')
legend(13, 15, paste0("p-value = ", round(pvaluesDepAUC[[factor]],2)), pch=NA, xjust=0.5, yjust=0.5, x.intersp=0, cex=1.3, bg='white')


# Daytime Hypo #
factor <- as.character(nullDep$factor)[13]
nullDep_auc <- c(unlist(nullDep[nullDep$factor == factor, seq(2, 201, 2)]), diabetes_vimp[diabetes_vimp$vimp == "arDep_ex_auc", factor])
x <- hist(unlist(nullDep_auc), breaks=20, plot=FALSE)

par(mar=c(2.2, 3, 0, 0.1) + 0.1)
plot(x, xlim=c(-15, 35), ylim=c(0,21), axes=FALSE, xlab="", ylab="", main="")
axis(1, at=seq(-15, 35, 5), tick=TRUE, cex.axis=1.3, line=-0.5, mgp=c(3, .6, 0))
axis(1, at=c(-10, 10, 20, 30), tick=TRUE, cex.axis=1.3, line=-0.5, mgp=c(3, .6, 0))
axis(2, at=seq(0,20,5), las=1, font.axis=2, mgp=c(3, .6, 0), xpd=TRUE)
mtext("Frequency", 2, font=2, line=1.6, cex=1.2)
mtext("Importance for Daytime Hypo", 1, font=2, line=1.2, cex=1.1)

polygon(rep(c(32,34), each=2), c(0,1,1,0), col='darkgray')
legend(24, 15, paste0("p-value = ", round(pvaluesDepAUC[[factor]],2)), pch=NA, xjust=0.5, yjust=0.5, x.intersp=0, cex=1.3, bg='white')

dev.off()
