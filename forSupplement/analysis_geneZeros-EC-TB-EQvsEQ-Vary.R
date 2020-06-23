setwd("EqualizationManuscript")

load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(ggsci)
library("scales")
library(SingleCellExperiment)


sub637.ec <- counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "EQ")]
sub637.tb <- counts(filtered)[,which(filtered$CellType == "TB" & filtered$Experiment == "EQ")]
sub644.ec <- counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "EQ-Vary")]
sub644.tb <- counts(filtered)[,which(filtered$CellType == "TB" & filtered$Experiment == "EQ-Vary")]


## Sub637 versus Sub644, EC
## Use only the cells in common.
ec.cells.qc <- intersect(colnames(sub637.ec), colnames(sub644.ec))

set.seed(444)
sum(colSums(sub637.ec[,ec.cells.qc]) - colSums(sub644.ec[,ec.cells.qc]) > 0)
sum(colSums(sub637.ec[,ec.cells.qc]) - colSums(sub644.ec[,ec.cells.qc]) < 0)
use.ec <- sample(names(which((colSums(sub637.ec[,ec.cells.qc]) - colSums(sub644.ec[,ec.cells.qc]) > 0))), 13)
use.ec <- c(use.ec, names(which(colSums(sub637.ec[,ec.cells.qc]) - colSums(sub644.ec[,ec.cells.qc]) < 0)))


useG <- intersect(rownames(sub637.ec), rownames(sub644.ec))

MedExpr <- log(apply(counts(filtered)[,which(filtered$Experiment == "Unequalized")], 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr)
grps <- length(splitby) / 4
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

countsNEW = sub637.ec[useG,use.ec]
countsOLD = sub644.ec[useG,use.ec]

zeroNEW <- apply(countsNEW, 1, function(x) mean(x!=0))
zeroOLD <- apply(countsOLD, 1, function(x) mean(x!=0))

meanDiff <- sort(zeroNEW - zeroOLD)


diffList <- list()
for (i in length(sreg_orig1):1) {
  topgenes <- names(sreg_orig1[[i]])
  diffList[[i]] <- meanDiff[topgenes]
}


pdf("PLOTS/SuppFig3A_changeInZeros_proportion_EQvsEQ-Vary_EC.pdf", height=4, width=12, useDingbats = F)
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


detectG1 <- rowMeans(sub637.ec[useG,use.ec]==0)
sum(detectG1 == 0)

detectG2 <- rowMeans(sub644.ec[useG,use.ec]==0)
sum(detectG2 == 0)

(sum(detectG2 == 0) - sum(detectG1 == 0)) / sum(detectG1 == 0)



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


XX = c(diffList[[2]], diffList[[3]], diffList[[4]])
XX = XX[XX > 0]
mean(XX)


XX = c(mean(diffList[[1]] < 0) /mean(diffList[[1]] > 0), 
       mean(diffList[[2]] < 0) /mean(diffList[[2]] > 0),
       mean(diffList[[3]] < 0) /mean(diffList[[3]] > 0),
       mean(diffList[[4]] < 0) /mean(diffList[[4]] > 0))

sd(XX) / mean(XX)


########################################################################################
########################################################################################

## Sub637 versus Sub644, EC
sum(colSums(sub644.tb[,]) - colSums(sub637.tb[,]) < 0)
sum(colSums(sub644.tb[,]) - colSums(sub637.tb[,]) > 0)


set.seed(77122)
use.tb <- sample(names(which((colSums(sub644.tb[,]) - colSums(sub637.tb[,]) < 0))), 22)
use.tb <- c(use.tb, names(which(colSums(sub644.tb[,]) - colSums(sub637.tb[,]) > 0)))

MedExpr <- log(apply(counts(filtered)[,which(filtered$Experiment == "Unequalized")], 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr)
grps <- length(splitby) / 4
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

countsNEW = sub637.tb[useG,use.tb]
countsOLD = sub644.tb[useG,use.tb]


zeroNEW <- apply(countsNEW, 1, function(x) mean(x!=0))
zeroOLD <- apply(countsOLD, 1, function(x) mean(x!=0))

meanDiff <- sort(zeroNEW - zeroOLD)


diffList <- list()
for (i in length(sreg_orig1):1) {
  topgenes <- names(sreg_orig1[[i]])
  diffList[[i]] <- meanDiff[topgenes]
}


pdf("PLOTS/SuppFig3B_changeInZeros_proportion_EQvsEQ-Vary_TB.pdf", height=4, width=12, useDingbats = F)
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

XX = c(diffList[[2]], diffList[[3]], diffList[[4]])
XX = XX[XX > 0]
mean(XX)


mean(diffList[[1]] == 0)


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
