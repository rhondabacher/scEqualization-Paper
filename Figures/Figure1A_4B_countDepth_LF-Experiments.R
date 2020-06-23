## This makes a Figure 1A, Figure 4B and a Supplemental Figure.

setwd("EqualizationManuscript")

load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(SCnorm)

BiocParallel::register(BiocParallel::SerialParam())
source("CODE/functionsForCountDepthPlots.R")

library(SingleCellExperiment)

## For figure partial:
pdf("PLOTS/countDepthPlots_EC_OrigAndEQ.pdf", height=7, width=6, useDingbats=F)

ModeStat_SquaredError <- list()

par(mfrow=c(2,2), mar=c(5,5,2,1))
useCounts <- round(counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "Unequalized")])
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))
i = 1
gslopes0 <- getSlopes(Data = useCounts)
head(gslopes0)
makePlot(gslopes0, "EC2 - Original", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes0, sreg_reseq1, "EC2 - Original")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")

useCounts <- round(counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "EQ")])
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0])))
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))
i = i + 1
gslopes1 <- getSlopes(Data = useCounts) 
makePlot(gslopes1, "EC2 - Sub637", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "EC2 - Sub637")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")
sapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1)))
## [1] 0.6536204 0.3307241
dev.off()





BiocParallel::register(BiocParallel::SerialParam())
# EC first
useCounts <- counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "Unequalized")]
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))

ModeStat_SquaredError <- list()

pdf("PLOTS/countDepthPlots_EC_allSets.pdf", height=12, width=8, useDingbats=F)
par(mfrow=c(4,2), mar=c(5,5,2,1))
i = 1
gslopes1 <- getSlopes(Data = useCounts)
makePlot(gslopes1, "Unequalized", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "Unequalized")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")

useCounts <- counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "EQ")]
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))

i = i + 1
gslopes1 <- getSlopes(Data = useCounts) 
makePlot(gslopes1, "EQ", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "EQ")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")

useCounts <- counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "EQ-Half")]
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0])))
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))

i = i + 1
gslopes1 <- getSlopes(Data = useCounts) 
makePlot(gslopes1, "EQ-Half", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "EQ-Half")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")

useCounts <- counts(filtered)[,which(filtered$CellType == "EC" & filtered$Experiment == "EQ-Vary")]
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))

i = i + 1
gslopes1 <- getSlopes(Data = useCounts)
makePlot(gslopes1, "EQ-Vary", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "EQ-Vary")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")



dev.off()




sapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1)))

ModeStatistic_EC <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))
ModeStatistic_EC



# now TB
useCounts <- counts(filtered)[,which(filtered$CellType == "TB" & filtered$Experiment == "Unequalized")]
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))

ModeStat_SquaredError <- list()

pdf("PLOTS/countDepthPlots_TB_allSets.pdf", height=12, width=8, useDingbats=F)
par(mfrow=c(4,2), mar=c(5,5,2,1))
i = 1
gslopes1 <- getSlopes(Data = useCounts)
makePlot(gslopes1, "TB2 - Original", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "TB2 - Original")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")

useCounts <- counts(filtered)[,which(filtered$CellType == "TB" & filtered$Experiment == "EQ")]
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))


i = i + 1
gslopes1 <- getSlopes(Data = useCounts) 
makePlot(gslopes1, "EQ", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "EQ")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")

useCounts <- counts(filtered)[,which(filtered$CellType == "TB" & filtered$Experiment == "EQ-Half")]
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))

i = i + 1
gslopes1 <- getSlopes(Data = useCounts) 
makePlot(gslopes1, "EQ-Half", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "EQ-Half")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")

useCounts <- counts(filtered)[,which(filtered$CellType == "TB" & filtered$Experiment == "EQ-Vary")]
MedExpr <- log(apply(useCounts, 1, function(x) median(x[x!=0]))) 
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_reseq1 <- split(splitby, ceiling(seq_along(splitby) / grps))


i = i + 1
gslopes1 <- getSlopes(Data = useCounts)
makePlot(gslopes1, "EQ-Vary", sreg_reseq1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_reseq1, "EQ-Vary")
abline(v=c(0,.5), lty=3, lwd=1.5, col="gray6")

dev.off()


ModeStatistic_TB <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))
ModeStatistic_TB

save(ModeStatistic_TB, ModeStatistic_EC, file="RDATA/countDepth_Stats_Expr_QCscater_trimmed_genes.Rdata")

## Now make a plots with the other datasets!
load("RDATA/countDepth_Stats_Expr_QCscater_trimmed_genes.Rdata")

library(ggsci)
library("scales")

useCols.Main <- pal_npg("nrc", alpha = .5)(3)
useCols <- rep(c(useCols.Main[1], useCols.Main[2], useCols.Main[3],
                 useCols.Main[2]), 2)
exprStat <- c(ModeStatistic_EC[c(1,2,4,3)], ModeStatistic_TB[c(1,2,4,3)])

pdf("PLOTS/MAD_LF-Experiment.pdf", height=2, width=2, useDingbats=FALSE)
par(mar=c(2,2,.1,.1))
plot(c(1.4,1.5,1.6,1.7, 2.5,2.6,2.7,2.8), exprStat, 
        pch=rep(c(17,17,17,19), 2), 
        col=useCols, xaxt='n', yaxt='n', bty='n',
        xlab="", ylab = "", cex=1, cex.axis=1, 
        ylim=c(0, .9), xlim=c(1.1,3.1))
axis(1, at=c(1.5, 2.5), cex.axis=.7, label=c("EC2", "TB2"))
axis(2, at=c(0, .3, .6, .9), cex.axis=.7, las=2)
dev.off()

