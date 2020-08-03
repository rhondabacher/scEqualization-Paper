setwd("EqualizationManuscript")

load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(ggsci)
library("scales")
library(SingleCellExperiment)


## gene-level detection

sub637.ec <- counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "EQ")]
sub637.tb <- counts(filtered)[,which(filtered$CellType == "TB" & filtered$Experiment == "EQ")]
orig.all.ec <- counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "Unequalized")]
orig.all.tb <- counts(filtered)[,which(filtered$CellType == "TB" & filtered$Experiment == "Unequalized")]

set.seed(6771)

## Use only the cells in common.
ec.cells.qc <- intersect(colnames(sub637.ec), colnames(orig.all.ec))

## Make sure that any difference in depth is not the cause...so pick an equal number that is bigger/smaller. Balance!
sum(colSums(sub637.ec[,ec.cells.qc]) - colSums(orig.all.ec[,ec.cells.qc]) > 0)
sum(colSums(sub637.ec[,ec.cells.qc]) - colSums(orig.all.ec[,ec.cells.qc]) < 0)
use.ec <- sample(names(which((colSums(sub637.ec[,ec.cells.qc]) - colSums(orig.all.ec[,ec.cells.qc]) > 0))), 16)
use.ec <- c(use.ec, names(which(colSums(sub637.ec[,ec.cells.qc]) - colSums(orig.all.ec[,ec.cells.qc]) < 0)))


## Sub637 versus Orig, EC
useG <- intersect(rownames(sub637.ec), rownames(orig.all.ec))

MedExpr <- log(apply(counts(filtered)[,which(filtered$Experiment == "Unequalized")], 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr)
grps <- length(splitby) / 4
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

countsNEW = sub637.ec[useG,use.ec]
countsOLD = orig.all.ec[useG,use.ec]


zeroNEW <- apply(countsNEW, 1, function(x) mean(x!=0))
zeroOLD <- apply(countsOLD, 1, function(x) mean(x!=0))

meanDiff <- sort(zeroNEW - zeroOLD)

## The signs here are flipped. More zeros in 637 means the meanDiff is negative.
sub637.more.zeros.ec <- meanDiff[which(meanDiff < 0)]
sub637.less.zeros.ec <- rev(meanDiff[which(meanDiff > 0)])


write.table(sub637.more.zeros.ec, file="OUT/genes_moreZerosin637_EC.csv", sep=",",col.names=F, quote=F)
write.table(sub637.less.zeros.ec, file="OUT/genes_lessZerosin637_EC.csv", sep=",",col.names=F, quote=F)

diffList <- list()
for (i in length(sreg_orig1):1) {
	topgenes <- names(sreg_orig1[[i]])
  diffList[[i]] <- meanDiff[topgenes]
}

pdf("PLOTS/changeInZeros_proportion_Sub637_EC_QCscater.pdf", height=4, width=12, useDingbats = FALSE)
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


detectG1 <- rowMeans(orig.all.ec[useG,ec.cells.qc]!=0)
sum(detectG1 == 1)

detectG2 <- rowMeans(sub637.ec[useG,ec.cells.qc]!=0)
sum(detectG2 == 1)

(sum(detectG2 == 1) - sum(detectG1 == 1)) / sum(detectG1 == 1)
(sum(detectG2 == 1) - sum(detectG1 == 1))
(sum(detectG2 == 0) - sum(detectG1 == 0)) / sum(detectG1 == 0)
(sum(detectG2 == 0) - sum(detectG1 == 0))

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

length(intersect(names(sub637.more.zeros.ec)[1:200], names(diffList[[1]])))
########################################################################################
########################################################################################
tb.cells.qc <- intersect(colnames(sub637.tb), colnames(orig.all.tb))


sum(colSums(sub637.tb[,tb.cells.qc]) - colSums(orig.all.tb[,tb.cells.qc]) > 0)
sum(colSums(sub637.tb[,tb.cells.qc]) - colSums(orig.all.tb[,tb.cells.qc]) < 0)


set.seed(6800)
use.tb <- sample(names(which((colSums(sub637.tb[,tb.cells.qc]) - colSums(orig.all.tb[,tb.cells.qc]) > 0))), 21)
use.tb <- c(use.tb, names(which(colSums(sub637.tb[,tb.cells.qc]) - colSums(orig.all.tb[,tb.cells.qc]) < 0)))

sum(colSums(sub637.tb[,use.tb]) - colSums(orig.all.tb[,use.tb]) > 0)
sum(colSums(sub637.tb[,use.tb]) - colSums(orig.all.tb[,use.tb]) < 0)


## Sub637 versus Orig, TB
MedExpr <- log(apply(orig.all.tb[,use.tb], 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr)
grps <- length(splitby) / 4
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

countsNEW = sub637.tb[useG,use.tb]
countsOLD = orig.all.tb[useG,use.tb]

zeroNEW <- apply(countsNEW, 1, function(x) mean(x!=0))
zeroOLD <- apply(countsOLD, 1, function(x) mean(x!=0))

meanDiff <- sort(zeroNEW - zeroOLD)

## The signs here are flipped more zeros in 637 means the meanDiff is negative.
sub637.more.zeros.tb <- meanDiff[which(meanDiff < 0)]
sub637.less.zeros.tb <- rev(meanDiff[which(meanDiff > 0)])

write.table(sub637.more.zeros.tb, file="OUT/genes_moreZerosin637_TB.csv", sep=",",col.names=F, quote=F)
write.table(sub637.less.zeros.tb, file="OUT/genes_lessZerosin637_TB.csv", sep=",",col.names=F, quote=F)


diffList <- list()
for (i in length(sreg_orig1):1) {
	topgenes <- names(sreg_orig1[[i]])
  diffList[[i]] <- meanDiff[topgenes]
}

length(intersect(names(sub637.more.zeros.tb)[1:200], names(diffList[[1]])))

pdf("PLOTS/changeInZeros_proportion_Sub637_TB_QCscater.pdf", height=4, width=12, useDingbats = FALSE)
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


detectG1 <- rowMeans(orig.all.tb[useG,tb.cells.qc]!=0)
sum(detectG1 == 1)

detectG2 <- rowMeans(sub637.tb[useG,tb.cells.qc]!=0)
sum(detectG2 == 1)

(sum(detectG2 == 1) - sum(detectG1 == 1)) / sum(detectG1 == 1)
(sum(detectG2 == 0) - sum(detectG1 == 0)) / sum(detectG1 == 0)
(sum(detectG2 == 1) - sum(detectG1 == 1))
(sum(detectG2 == 0) - sum(detectG1 == 0))


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


