## Calculate seq depth and detection rates for experiments:

setwd("EqualizationManuscript")

load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(SingleCellExperiment)


## Difference in the proportion of zeros per cell across experiments
## cell-specific detection rates
(mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ")]!=0)) - mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))) /  
  mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))

(mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Vary")]!=0)) - mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))) /  
  mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))

(mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Half")]!=0)) - mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ")]!=0))) /  
  mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ")]!=0))

(mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Half")]!=0)) - mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))) /  
  mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))

mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ")]!=0)) - mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))
mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Vary")]!=0)) - mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))
mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Half")]!=0)) - mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ")]!=0))
mean(colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Half")]!=0)) - mean(colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0))


## These are the difference in genes that are totally detected versus totally undetected:
zeroNEW <- apply(counts(filtered)[,which(filtered$CellType =="EC" & filtered$Experiment == "EQ")], 1, function(x) mean(x!=0))
zeroOLD <- apply(counts(filtered)[,which(filtered$CellType =="EC" & filtered$Experiment == "Unequalized")], 1, function(x) mean(x!=0))

sum(zeroNEW == 1) - sum(zeroOLD == 1)

sum(zeroNEW == 0) - sum(zeroOLD == 0)

(sum(zeroNEW == 1) - sum(zeroOLD == 1)) / sum(zeroOLD == 1)

(sum(zeroNEW == 0) - sum(zeroOLD == 0)) / sum(zeroOLD == 0)


zeroNEW <- apply(counts(filtered)[,which(filtered$CellType =="TB" & filtered$Experiment == "EQ")], 1, function(x) mean(x!=0))
zeroOLD <- apply(counts(filtered)[,which(filtered$CellType =="TB" & filtered$Experiment == "Unequalized")], 1, function(x) mean(x!=0))

sum(zeroNEW == 1) - sum(zeroOLD == 1)

sum(zeroNEW == 0) - sum(zeroOLD == 0)

(sum(zeroNEW == 1) - sum(zeroOLD == 1)) / sum(zeroOLD == 1)

(sum(zeroNEW == 0) - sum(zeroOLD == 0)) / sum(zeroOLD == 0)


# Detection rates:
X1 <- data.frame( Depth = (colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]!=0)) , Dataset = "Unequalized")
X2 <- data.frame( Depth = (colSums(counts(filtered)[,which(filtered$Experiment == "EQ")]!=0)) , Dataset = "EQ")
X3 <- data.frame( Depth = (colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Vary")]!=0)), Dataset = "EQ-Vary")
X4 <- data.frame( Depth = (colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Half")]!=0)), Dataset = "EQ-Half")

library(ggsci)
library("scales")

useCols <- pal_npg("nrc", alpha = 1)(4)

pdf("PLOTS/Figure_DetectionTotalCell.pdf", height=3.5, width=5, useDingbats=F)
longdata <- rbind(X1,X2, X3, X4)
par(mfrow=c(1,1))
par(mar=c(2,6,1,1), mgp = c(4.1, .5, 0))
yarrr::pirateplot(formula = Depth ~ Dataset,
           data = longdata,
           xlab = "", inf.method = "iqr", avg.line.fun = median, ylim=c(6000, 12000),
           ylab = "# Genes", pal=useCols[1:4],inf.f.o=.3,
           main = "", point.cex=1, bar.lwd=2, cex.lab=1.2, cex.axis=1.2,cex.names=1.2)
dev.off()


# Depths rates:
X1 <- data.frame( Depth = (colSums(counts(filtered)[,which(filtered$Experiment == "Unequalized")]))/1e6 , Dataset = "Unequalized")
X2 <- data.frame( Depth = (colSums(counts(filtered)[,which(filtered$Experiment == "EQ")]))/1e6 , Dataset = "EQ")
X3 <- data.frame( Depth = (colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Vary")]))/1e6, Dataset = "EQ-Vary")
X4 <- data.frame( Depth = (colSums(counts(filtered)[,which(filtered$Experiment == "EQ-Half")]))/1e6, Dataset = "EQ-Half")

library(ggsci)
library("scales")

# useCols.Main <- pal_npg("nrc", alpha = .5)(3)
# useCols <- rep(c(useCols.Main[1], useCols.Main[2], useCols.Main[3],
#                  useCols.Main[2]), 2)

pdf("PLOTS/forSupp_DepthTotalCell.pdf", height=3, width=4, useDingbats=F)
longdata <- rbind(X1,X2, X3, X4)
# longdata[,1] <- log(longdata[,1])
par(mfrow=c(1,1))
par(mar=c(2,3,1,1), mgp = c(2, .5, 0))
yarrr::pirateplot(formula = Depth ~ Dataset,
           data = longdata,
           xlab = "", inf.method = "iqr", avg.line.fun = median, ylim=c(0, 15),
           ylab = "Reads (per million)", pal=useCols[1:6],inf.f.o=.3,
           main = "Sequencing Depth per Cell", point.cex=1, bar.lwd=2, cex.lab=1.2, cex.axis=1.2,cex.names=1.2)
dev.off()


sd(X1$Depth) / mean(X1$Depth)
sd(X2$Depth) / mean(X2$Depth)
sd(X3$Depth) / mean(X3$Depth)


median(abs(X1$Depth - median(X1$Depth)))
median(abs(X2$Depth - median(X2$Depth)))
median(abs(X3$Depth - median(X3$Depth)))

