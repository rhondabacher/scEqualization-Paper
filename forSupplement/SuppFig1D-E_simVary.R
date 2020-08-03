# Make figures of simulations ran on hipergator:
setwd("EqualizationManuscript")

load("RDATA/simulateFromMatching_calcMAD_varyDepth_trimmed.RData")

pcntRange <- c(.25,.5,1,2,3)

library(ggsci)
library("scales")

useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,trackMAD)))
exprStat <- do.call(c,trackMAD)
eqVal <- rep(1:5, each=length(trackMAD[[1]]))

pdf("PLOTS/MAD_simSet_varyDepth_Statistic_trimmed.pdf", height=2, width=3)
par(mar=c(2,2,2,1))
plot(eqVal, exprStat, pch=16, col=useCols.Dots, xaxt='n', yaxt='n',
     xlab="", ylab = "Median Absolute Deviation", cex=.5, cex.lab=1, ylim=c(.2, 1), xlim=c(.5,5.5))
axis(2, cex.axis=1)
axis(1, at=c(1:5), label=as.character(pcntRange), cex.axis=1)
useCols.Main <- pal_npg("nrc", alpha = .8)(1)
points(unique(eqVal), sapply(trackMAD, mean), pch="-", cex=3, col=useCols.Main)
abline(h=seq(0,max(exprStat)+.2, by=.1), lwd=1, lty=3, col="gray80")
dev.off()

############################################################################################################
rm(list=ls())

load("RDATA/simulateFromMatching_calcMAD_varyEQ_trimmed_try.RData")

library(ggsci)
library("scales")

trackMAD <- trackMAD[-4]
pcntRange <- c(0,.125,.25,.50,-1)
pcntRange <- pcntRange[-4]

useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,trackMAD)))
exprStat <- do.call(c,trackMAD)
eqVal <- rep(1:4, each=length(trackMAD[[1]]))

pdf("PLOTS/MAD_simSet_varyEQ_Statistic_trimmed.pdf", height=2, width=3)
par(mar=c(2,2,2,1))
plot(eqVal, exprStat, pch=16, col=useCols.Dots, xaxt='n', yaxt='n',
     xlab="", ylab = "Median Absolute Deviation", cex=.5, cex.lab=1, ylim=c(0, 1), xlim=c(.5,4.5))
axis(2, cex.axis=1)
axis(1, at=c(1:4), label=as.character(pcntRange), cex.axis=1)
useCols.Main <- pal_npg("nrc", alpha = .8)(1)
points(unique(eqVal), sapply(trackMAD, mean), pch="-", cex=3, col=useCols.Main)
abline(h=seq(0,max(exprStat)+.2, by=.1), lwd=1, lty=3, col="gray80")
dev.off()
