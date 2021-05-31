
library(quantreg)
library(SCnorm)
library(parallel)
library(scaffold) 

ss3.data <- read.table("DATA/Smartseq3/Smartseq3.Fibroblasts.NovaSeq.readcounts.txt", header=TRUE)
mean(colSums(ss3.data)) / 1e6

## check QC
library(scater)
ss3data.sce <- SingleCellExperiment(assay=list(counts = data.matrix(ss3.data)))
example_sce <- addPerCellQC(ss3data.sce)
plotColData(example_sce, x = "total", y="detected")
keep.total.1 <- isOutlier(example_sce$sum, nmads=5, type="both")
keep.total.1 <- isOutlier(example_sce$detected, nmads=5, type="both")
sum(keep.total.1)
# no outliers

dim(ss3.data)
mean(colSums(ss3.data)) / 1e6
mean(colMeans(ss3.data!=0))

source("CODE/functionsForCountDepthPlots.R")
source("CODE/forRevisions/fasterSloperCalc_Helper.R")

gslopes0 <- fasterSlope(Data = data.matrix(ss3.data), ditherCounts=FALSE)

save.image("RDATA/BestMatchingSim_SS3_Fibroblast-reads.RData")

set.seed(5444)
RcppZiggurat::zsetseed(87)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = ss3.data))
params <- estimateScaffoldParameters(sce, equalizationAmount = -1, sceUMI=FALSE, useUMI = TRUE)
simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce)

save.image("RDATA/Sim_SS3_Fibroblast-reads.RData")

load("RDATA/Sim_SS3_Fibroblast-reads.RData")

# Let's first compare the counts
countsSim <- SingleCellExperiment::counts(simdata)
compareData <- ss3.data

gslopes1 <- fasterSlope(Data = data.matrix(countsSim), ditherCounts = FALSE)

save.image("RDATA/Sim_SS3_Fibroblast-reads.RData")

pdf("PLOTS/simMatching_CountDepth_SS3-Fibroblast-reads.pdf", height=8, width=8, useDingbats=F)
ModeStat_SquaredError<-list()
MedExpr <- (apply(ss3.data, 1, function(x) mean(x)))
splitby <- sort(MedExpr[MedExpr > 1])
grps <- length(splitby) / 10
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

par(mfrow=c(2,2), mar=c(5,5,2,1))

i = 1
makePlot(gslopes0, " ", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes0, sreg_orig1, " ")

i = i + 1
makePlot(gslopes1, " ", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_orig1, " ")

ModeStatistic_EC <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))

print(paste0(ModeStatistic_EC))
## "0.0958904109589043" "0.236790606653621" 
dev.off()


library(yarrr)


pdf("PLOTS/simMatching_SeqDepth_SS3-Fibroblast-reads.pdf", height=4, width=8, useDingbats=F)
X <- data.frame( Depth = colSums(countsSim) , Species = "Simulated")
Y <- data.frame( Depth = colSums(compareData) , Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata, 
           xlab = "", 
           ylab = "Sequencing Depth", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Sequencing Depth (millions)", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()




pdf("PLOTS/simMatching_CellDetectRate_SS3-Fibroblast-reads.pdf", height=4, width=8, useDingbats=F)
X <- data.frame( Depth = colSums(countsSim!=0) / nrow(countsSim), Species = "Simulated")
Y <- data.frame( Depth = colSums(compareData!=0) / nrow(countsSim), Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species, 
           data = longdata, ylim=c(0,1),
           xlab = "", 
           ylab = "Detection Rate", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Detection Rate", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("PLOTS/simMatching_GeneDetectRate_SS3-Fibroblast-reads.pdf", height=4, width=8, useDingbats=F)
set.seed(111)
Genes = rownames(countsSim)
XX <- sample(Genes, 200)
X1 <- apply(countsSim[XX,], 1, function(x) sum(x!=0, na.rm=T)) / dim(countsSim)[2]
X2 <- apply(compareData[XX,], 1, function(x) sum(x!=0, na.rm=T)) / dim(compareData)[2]
X <- data.frame( Depth =X1, Species = "Simulated")
Y <- data.frame( Depth = X2, Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "Detection Rate", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Detection Rate", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("PLOTS/simMatching_GeneMean_SS3-Fibroblast-reads.pdf", height=4, width=8, useDingbats=F)
X1 <- log(apply(countsSim[XX,], 1, function(x) mean(x, na.rm=T))+1)
X2 <- log(apply(compareData[XX,], 1, function(x) mean(x, na.rm=T))+1)
useg <- names(which(X1<Inf & X2 < Inf))
X <- data.frame( Depth =X1[useg], Species = "Simulated")
Y <- data.frame( Depth = X2[useg], Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "log (mean+1)", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="log (mean+1)", 
     cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("PLOTS/simMatching_GeneSD_SS3-Fibroblast-reads.pdf", height=4, width=8, useDingbats=F)
X1 <- log(apply(countsSim[XX,], 1, function(x) sd(x, na.rm=T))+1)
X2 <- log(apply(compareData[XX,], 1, function(x) sd(x, na.rm=T))+1)
useg <- names(which(X1<Inf & X2 < Inf))
X <- data.frame( Depth =X1[useg], Species = "Simulated")
Y <- data.frame( Depth = X2[useg], Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "log (sd+1)", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="log (sd+1)", 
     cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()







