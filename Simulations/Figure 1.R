#rm(list=ls(all=TRUE))

# Figure 1 sets beta = log(2)
sim_onesplit <- read.csv(paste0("Data/One Split/sim_onesplit_night_cutpt3_betaLog2.csv"), header=TRUE)

# Manuscript text also reports results when beta = log(4).  Can load following data instead to show those results
#sim_onesplit <- read.csv(paste0("Data/One Split/sim_onesplit_night_cutpt3_betaLog4.csv"), header=TRUE)

metrics <- c(paste0("treeDep.",1:4), paste0("treeIndep.",1:4))
sim_pctSims <- data.frame(lapply(metrics, function(x){setNames(aggregate(list(sim_onesplit[[x]]), list(sim_onesplit$model), FUN=mean), c("Model", x))}))[c(1,seq(2,16,2))]
sim_pctSims[ , c(1, 2, 6)] #Percentage of simulations with informative variable split


# Create Figure 1:

#tiff("Figures/Fig1.tiff", width = 4, height = 4, units = "in", pointsize = 12, res=600)
par(mar=c(1.6, 3.9, 0, 0) + 0.1)
plot(NA, xlim=c(1.8,10.7), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
mtext("Model", 1, line=0.8, font=2, cex=1.3)
mtext("% Simulations with", 2, line=3, font=2, cex=1.2)
mtext("Informative Variable Split", 2, line=2, font=2, cex=1.2)
for (j in seq.int(5)){
  polygon(rep(c(2*j-0.2, 2*j+0.2), each=2), c(0, sim_pctSims$treeDep.1[j], sim_pctSims$treeDep.1[j], 0), col="gray45")
  polygon(rep(c(2*j+0.3, 2*j+0.7), each=2), c(0, sim_pctSims$treeIndep.1[j], sim_pctSims$treeIndep.1[j], 0), col="gray80")
}
axis(1, at=seq(2.25,10.25,2), labels = LETTERS[1:5], tick=FALSE, line=-1.3)
axis(2, at=seq(0,1,.25), labels = paste0(seq(0,1,.25)*100,'%'), tick=TRUE, las=2, line=-.4, mgp=c(3, .6, 0))
legend(3.5, 1.05, c("Robust","Naïve"),pch=22, pt.bg=c("gray45","gray80"), bty='n', pt.cex=1.8, cex=1, adj=c(0,0.5))
dev.off()


