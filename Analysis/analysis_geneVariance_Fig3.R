setwd("EqualizationManuscript")

load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

sub637 <- counts(filtered)[,which(filtered$Experiment == "EQ")]
orig.all <- counts(filtered)[,which(filtered$Experiment == "Unequalized")]

library(SingleCellExperiment)
genes <- intersect(rownames(sub637), rownames(orig.all))

# Make into SCE objects
chu.1 <- SingleCellExperiment(assays = list(counts = orig.all))
chu.2 <- SingleCellExperiment(assays = list(counts = sub637))


# library(scran)
# logcounts(chu.1) <- log(counts(chu.1) + 1)
# logcounts(chu.2) <- log(counts(chu.2) + 1)

# Normalize first
library(SCnorm)
Conditions <- c(rep("EC", length(grep("EC", colnames(orig.all)))), rep("TB", length(grep("TB", colnames(orig.all)))))
norm.chu <- SCnorm(Data = chu.1 , Conditions=Conditions, reportSF=T)

BiocParallel::register(BiocParallel::SerialParam())

Conditions <- c(rep("EC", length(grep("EC", colnames(sub637)))), rep("TB", length(grep("TB", colnames(sub637)))))
norm.chu.eq <- SCnorm(Data = chu.2, Conditions=Conditions, reportSF=T)

save.image("RDATA/diffVariance_readyNormalized_trimmed.RData")

load("RDATA/diffVariance_readyNormalized_trimmed.RData")

norm.chu.mat <- normcounts(norm.chu)
norm.chu.mat.log <- log(norm.chu.mat + 1)

norm.chu.eq.mat <- normcounts(norm.chu.eq)
norm.chu.mat.eq.log <- log(norm.chu.eq.mat + 1)

logcounts(chu.1) <- norm.chu.mat.log
logcounts(chu.2) <- norm.chu.mat.eq.log

# Decompose variance
library(scran)
alt.fit.1 <- trendVar(chu.1, use.spikes=FALSE, parametric=FALSE, 
                      method="spline")
alt.decomp.1 <- decomposeVar(chu.1, alt.fit.1)

alt.fit.2 <- trendVar(chu.2, use.spikes=FALSE, parametric=FALSE, method="spline")
alt.decomp.2 <- decomposeVar(chu.2, alt.fit.2)

orig.hg <- alt.decomp.1
sub637.hg <- alt.decomp.2


pdf("PLOTS/meanVar-decompose_OrigVsEq.pdf", height=5, width=10)
par(mfrow=c(1,2))
plot(alt.decomp.1$mean, alt.decomp.1$total, pch=16, cex=0.6, ylim=c(0,25),
     xlab="Mean log-expression", ylab="Variance of log-expression", main="Unequalized")
curve(alt.fit.1$trend(x), col="dodgerblue", lwd=2, add=TRUE)

plot(alt.decomp.2$mean, alt.decomp.2$total, pch=16, cex=0.6, ylim=c(0,25),
     xlab="Mean log-expression", ylab="Variance of log-expression", main="EQ")
curve(alt.fit.2$trend(x), col="dodgerblue", lwd=2, add=TRUE)
dev.off()


##########################################################################
### Try difference to top HVG:

X <- setdiff(rownames(subset(sub637.hg, FDR < .1)), rownames(subset(orig.hg, FDR < .1)))
Y <- setdiff(rownames(subset(orig.hg, FDR < .1)), rownames(subset(sub637.hg, FDR < .1)))

length(X); length(Y)
write.table(X, file="OUT/genes_moreVarInEQexpr_topDiff.txt", row.names=F, quote=F, col.names=F)
write.table(Y, file="OUT/genes_lessVarInEQexpr_topDiff.txt", row.names=F, quote=F, col.names=F)

# More genes are more variable in Orig experiment. 
# Put these into enrichment.
##########################################################################

## Make Figure
library(ggsci)
library("scales")

useCols.Main <- pal_aaas("default", alpha = 1)(4)

useCols.Main <- pal_npg("nrc", alpha = 1)(2)

useg <- intersect(c(rownames(subset(orig.hg, FDR < .1))), rownames(subset(sub637.hg, FDR < .1)))
length(useg)

pdf("PLOTS/meanVar_rankDiffPlot-ZOOM_FDR.1.pdf", height=6, width=4, useDingbats=FALSE)
par(mfrow=c(2,1), mar=c(5,5,2,1))
plot(alt.decomp.1$mean, alt.decomp.1$total, pch=16, cex=0.2, col="gray85", ylim=c(0,20),
     main="Unequalized",
     xlab="Mean of log-expression", ylab="Variance of log-expression", cex.axis=1.2, cex.lab=1.2)
points(alt.decomp.1[useg, "mean"], alt.decomp.1[useg, "total"],  pch=19,cex=.5, col=alpha("gray55",.8))
curve(alt.fit.1$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(alt.decomp.1[Y, "mean"], alt.decomp.1[Y, "total"],  pch=19,cex=.8, col=useCols.Main[1])
points(alt.decomp.1[X, "mean"], alt.decomp.1[X, "total"],  pch=19,cex=.8, col=useCols.Main[2])


plot(alt.decomp.2$mean, alt.decomp.2$total, pch=16, cex=0.2, col="gray85", ylim=c(0,20),
     main="EQ",
     xlab="Mean of log-expression", ylab="Variance of log-expression", cex.axis=1.2, cex.lab=1.2)
points(alt.decomp.2[useg, "mean"], alt.decomp.2[useg, "total"],  pch=19,cex=.5, col=alpha("gray55",.8))
curve(alt.fit.1$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(alt.decomp.2[Y, "mean"], alt.decomp.2[Y, "total"],  pch=19,cex=.8, col=useCols.Main[1])
points(alt.decomp.2[X, "mean"], alt.decomp.2[X, "total"],  pch=19,cex=.8, col=useCols.Main[2])
dev.off()




## Check housekeeping genes??
# hk: https://www.tau.ac.il/~elieis/HKG/HK_genes.txt

## Rank the Biological variable for each dataset:
orig.hg$varRank <- rank(orig.hg$bio)
sub637.hg$varRank <- rank(sub637.hg$bio)

# Difference in ranks
rankDiff <- orig.hg$varRank - sub637.hg$varRank
names(rankDiff) <- rownames(orig.hg)

# Only consider genes with smaller FDR; HVG in either dataset
useg <- intersect(rownames(orig.hg), unique(c(rownames(subset(orig.hg, FDR < .1)), rownames(subset(sub637.hg, FDR < .1)))))
length(useg)

rankDiff <- sort(rankDiff[useg])
mean(rankDiff < 0)
mean(rankDiff > 0) # 58%

sum(rankDiff < 0)
sum(rankDiff > 0) 
length(rankDiff)
# 767 / 1320


hk.genes <- read.table("DATA/hskp-genes.txt", header=T, stringsAsFactors=F)[,1]

hk.genes <- intersect(useg, hk.genes)
length(hk.genes)
Q <- orig.hg[hk.genes,"varRank"] - sub637.hg[hk.genes,"varRank"]
sum(Q > 0, na.rm=T)
mean(Q > 0, na.rm=T)
mean(Q < 0, na.rm=T)
### 0.78
### 290 / 370
# Compared to XX% of genes overall!

sum(orig.hg[hk.genes,"varRank"] - sub637.hg[hk.genes,"varRank"] < 0, na.rm=T)
sum(orig.hg[hk.genes,"varRank"] - sub637.hg[hk.genes,"varRank"] > 0, na.rm=T)
## A lot more genes have higher variation in Original than in Sub637!

Q <- (sub637.hg[hk.genes,"varRank"] - orig.hg[hk.genes,"varRank"]) / orig.hg[hk.genes,"varRank"]
mean(Q)
mean(Q[Q < 0], na.rm=T)
mean(Q[Q > 0], na.rm=T)


