setwd("EqualizationManuscript")

load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(ggsci)
library("scales")
library(SingleCellExperiment)

orig.all <- counts(filtered)[,which(filtered$Experiment == "Unequalized")]

ec.cells.qc = grep("EC", colnames(orig.all), value=T)

## ORIG SPLIT

useG <- rownames(orig.all)
length(ec.cells.qc)
set.seed(1212)
split1 <- sample(ec.cells.qc, 20)
split2 <- sample(setdiff(ec.cells.qc, split1), 20)

sum(colSums(orig.all[,split1]) - colSums(orig.all[,split2]) < 0)
sum(colSums(orig.all[,split1]) - colSums(orig.all[,split2]) > 0)

MedExpr <- log(apply(orig.all[useG,ec.cells.qc], 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr)
grps <- length(splitby) / 4
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

countsNEW = orig.all[useG,split1]
countsOLD = orig.all[useG,split2]

zeroNEW <- apply(countsNEW, 1, function(x) mean(x!=0))
zeroOLD <- apply(countsOLD, 1, function(x) mean(x!=0))

meanDiff <- sort(zeroNEW - zeroOLD)

diffList <- list()
for (i in length(sreg_orig1):1) {
	topgenes <- names(sreg_orig1[[i]])
  diffList[[i]] <- meanDiff[topgenes]
}

pdf("PLOTS/changeInZeros_proportion_Sub637_EC_SPLIT-unEQ.pdf", height=4, width=12, useDingbats = F)
MIN = min(meanDiff) - .02 
MAX = max(meanDiff)  + .02
par(oma=c(1,6,1,1), las=1, mgp=c(5,2,0))

par(mfrow=c(1,4), mar=c(4,.4,.1,.4))

NEG <- sum(unique(diffList[[1]]) < 0)
POS <- sum(unique(diffList[[1]]) >= 0)
plot.ecdf(diffList[[1]], xlim=c(MIN, MAX), main="", ylab="", 
cex.axis=3, cex=1.5, ylim=c(-.01,1.01), col=c(rep("brown1", NEG), rep("dodgerblue3", POS)));abline(v=0)
abline(h=mean(diffList[[1]] < 0), lty=2, lwd=2)
abline(h=1 - mean(diffList[[1]] > 0), lty=2, lwd=2)

NEG <- sum(unique(diffList[[2]]) < 0)
POS <- sum(unique(diffList[[2]]) >= 0)
plot.ecdf(diffList[[2]], xlim=c(MIN, MAX), main="", ylab="", yaxt='n',
cex.axis=3, cex=1.5, ylim=c(-.01,1.01), col=c(rep("brown1", NEG), rep("dodgerblue3", POS)));abline(v=0)
abline(h=mean(diffList[[2]] < 0), lty=2, lwd=2)
abline(h=1 - mean(diffList[[2]] > 0), lty=2, lwd=2)

NEG <- sum(unique(diffList[[3]]) < 0)
POS <- sum(unique(diffList[[3]]) >= 0)
plot.ecdf(diffList[[3]], xlim=c(MIN, MAX), main="", ylab="", yaxt='n',
cex.axis=3, cex=1.5, ylim=c(-.01,1.01), col=c(rep("brown1", NEG), rep("dodgerblue3", POS)));abline(v=0)
abline(h=mean(diffList[[3]] < 0), lty=2, lwd=2)
abline(h=1 - mean(diffList[[3]] > 0), lty=2, lwd=2)

NEG <- sum(unique(diffList[[4]]) < 0)
POS <- sum(unique(diffList[[4]]) >= 0)
plot.ecdf(diffList[[4]], xlim=c(MIN, MAX), main="", ylab="", yaxt='n',
cex.axis=3, cex=1.5, ylim=c(-.01,1.01), col=c(rep("brown1", NEG), rep("dodgerblue3", POS)));abline(v=0)
abline(h=mean(diffList[[4]] < 0), lty=2, lwd=2)
abline(h=1 - mean(diffList[[4]] > 0), lty=2, lwd=2)

dev.off()

detectG <- rowMeans(orig.all[useG,ec.cells.qc]==0)

length(intersect(names(which(detectG == 0)), names(diffList[[1]]))) / length(diffList[[1]])
length(intersect(names(which(detectG == 0)), names(diffList[[2]]))) / length(diffList[[2]])
length(intersect(names(which(detectG == 0)), names(diffList[[3]]))) / length(diffList[[3]])
length(intersect(names(which(detectG == 0)), names(diffList[[4]]))) / length(diffList[[4]])



mean(diffList[[1]] == 0)
mean(diffList[[1]] < 0)
mean(diffList[[1]] > 0)

mean(diffList[[2]] == 0)
mean(diffList[[2]] < 0)
mean(diffList[[2]] > 0)

mean(diffList[[3]] == 0)
mean(diffList[[3]] < 0)
mean(diffList[[3]] > 0)

mean(diffList[[4]] == 0)
mean(diffList[[4]] < 0)
mean(diffList[[4]] > 0)

XX = c(mean(diffList[[1]] < 0) /mean(diffList[[1]] > 0), 
mean(diffList[[2]] < 0) /mean(diffList[[2]] > 0),
mean(diffList[[3]] < 0) /mean(diffList[[3]] > 0),
mean(diffList[[4]] < 0) /mean(diffList[[4]] > 0))

sd(XX) / mean(XX)


########################################################################################
########################################################################################
tb.cells.qc = grep("TB", colnames(orig.all), value=T)

## Sub637 versus Orig, TB
MedExpr <- log(apply(orig.all[useG,tb.cells.qc], 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr)
grps <- length(splitby) / 4
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

length(tb.cells.qc)
set.seed(821)
split1 <- sample(tb.cells.qc, 25)
split2 <- sample(setdiff(tb.cells.qc, split1), 25)

sum(colSums(orig.all[,split1]) - colSums(orig.all[,split2]) < 0)
sum(colSums(orig.all[,split1]) - colSums(orig.all[,split2]) > 0)

countsNEW = orig.all[useG,split1]
countsOLD = orig.all[useG,split2]

countsNEW <- t((t(countsNEW) / colSums(countsNEW)) * 1e6)
countsOLD <- t((t(countsOLD) / colSums(countsOLD)) * 1e6)

zeroNEW <- apply(countsNEW, 1, function(x) mean(x!=0))
zeroOLD <- apply(countsOLD, 1, function(x) mean(x!=0))

meanDiff <- sort(zeroNEW - zeroOLD)

## The signs here are flipped more zeros in 637 means the meanDiff is negative.
sub637.more.zeros.tb <- meanDiff[which(meanDiff < 0)]
sub637.less.zeros.tb <- meanDiff[which(meanDiff > 0)]



diffList <- list()
for (i in length(sreg_orig1):1) {
	topgenes <- names(sreg_orig1[[i]])
  diffList[[i]] <- meanDiff[topgenes]
}

pdf("PLOTS/changeInZeros_proportion_Sub637_TB_SPLIT-unEQ.pdf", height=4, width=12, useDingbats = F)
MIN = min(meanDiff) - .02 
MAX = max(meanDiff)  + .02
par(oma=c(1,6,1,1), las=1, mgp=c(5,2,0))

par(mfrow=c(1,4), mar=c(4,.4,.1,.4))

NEG <- sum(unique(diffList[[1]]) < 0)
POS <- sum(unique(diffList[[1]]) >= 0)
plot.ecdf(diffList[[1]], xlim=c(MIN, MAX), main="", ylab="", 
cex.axis=3, cex=1.5, ylim=c(-.01,1.01), col=c(rep("brown1", NEG), rep("dodgerblue3", POS)));abline(v=0)
abline(h=mean(diffList[[1]] < 0), lty=2, lwd=2)
abline(h=1 - mean(diffList[[1]] > 0), lty=2, lwd=2)

NEG <- sum(unique(diffList[[2]]) < 0)
POS <- sum(unique(diffList[[2]]) >= 0)
plot.ecdf(diffList[[2]], xlim=c(MIN, MAX), main="", ylab="", yaxt='n',
cex.axis=3, cex=1.5, ylim=c(-.01,1.01), col=c(rep("brown1", NEG), rep("dodgerblue3", POS)));abline(v=0)
abline(h=mean(diffList[[2]] < 0), lty=2, lwd=2)
abline(h=1 - mean(diffList[[2]] > 0), lty=2, lwd=2)

NEG <- sum(unique(diffList[[3]]) < 0)
POS <- sum(unique(diffList[[3]]) >= 0)
plot.ecdf(diffList[[3]], xlim=c(MIN, MAX), main="", ylab="", yaxt='n',
cex.axis=3, cex=1.5, ylim=c(-.01,1.01), col=c(rep("brown1", NEG), rep("dodgerblue3", POS)));abline(v=0)
abline(h=mean(diffList[[3]] < 0), lty=2, lwd=2)
abline(h=1 - mean(diffList[[3]] > 0), lty=2, lwd=2)

NEG <- sum(unique(diffList[[4]]) < 0)
POS <- sum(unique(diffList[[4]]) >= 0)
plot.ecdf(diffList[[4]], xlim=c(MIN, MAX), main="", ylab="", yaxt='n',
cex.axis=3, cex=1.5, ylim=c(-.01,1.01), col=c(rep("brown1", NEG), rep("dodgerblue3", POS)));abline(v=0)
abline(h=mean(diffList[[4]] < 0), lty=2, lwd=2)
abline(h=1 - mean(diffList[[4]] > 0), lty=2, lwd=2)

dev.off()

XX = diffList[[4]]
XX = XX[XX > 0]
mean(XX)


mean(diffList[[1]] < 0)
mean(diffList[[1]] > 0)

mean(diffList[[2]] < 0)
mean(diffList[[2]] > 0)

mean(diffList[[3]] < 0)
mean(diffList[[3]] > 0)

mean(diffList[[4]] < 0)
mean(diffList[[4]] > 0)


XX = c(mean(diffList[[1]] < 0) /mean(diffList[[1]] > 0), 
mean(diffList[[2]] < 0) /mean(diffList[[2]] > 0),
mean(diffList[[3]] < 0) /mean(diffList[[3]] > 0),
mean(diffList[[4]] < 0) /mean(diffList[[4]] > 0))

sd(XX) / mean(XX)

