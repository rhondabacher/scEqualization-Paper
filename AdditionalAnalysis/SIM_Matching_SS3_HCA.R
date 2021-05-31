library(quantreg)
library(SCnorm)
library(parallel)
library(scaffold) 

setwd("~/OneDrive - University of Florida/EqualizationManuscript")

ss3.data <- read.table("DATA/Smartseq3/HCA.readcounts.PBMC.txt", header=TRUE)

# ## check QC
library(scater)
ss3.sce <- SingleCellExperiment(assay=list(counts = data.matrix(ss3.data)))
example_sce <- addPerCellQC(ss3.sce)
plotColData(example_sce, x = "total", y="detected")

keep.total.1 <- isOutlier(example_sce$sum, nmads=5, type="both")
keep.total.2 <- isOutlier(example_sce$total, nmads=5, type="both")
keepcells <- keep.total.1 + keep.total.2
ss3.data <- ss3.data[,which(keepcells == 0)]
dim(ss3.data)

save(ss3.data, file="RDATA/dataReady_SS3-HCA-reads.RData")

dim(ss3.data)
mean(colSums(ss3.data)) / 1e6
mean(colMeans(ss3.data!=0))
# # # 

load("RDATA/dataReady_SS3-HCA-reads.RData")

source("CODE/functionsForCountDepthPlots.R")
source("CODE/fasterSloperCalc_Helper.R")

set.seed(1313)
RcppZiggurat::zsetseed(22)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = ss3.data))
params <- estimateScaffoldParameters(sce, equalizationAmount = .67, sceUMI = FALSE, 
                                     useUMI=TRUE)
simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce)

save.image("RDATA/Sim_SS3_HCA-reads.RData")

# Let's first compare the counts
countsSim <- SingleCellExperiment::counts(simdata)
compareData <- ss3.data

## Too many zeros to calclulate slopes accurately.

library(yarrr)

pdf("simMatching_SeqDepth_SS3-HCA-reads.pdf", height=4, width=8, useDingbats=F)
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
legend('bottomright', c("Simulated", "HCA Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()




pdf("simMatching_CellDetectRate_SS3-HCA-reads.pdf", height=4, width=8, useDingbats=F)
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
legend('bottomright', c("Simulated", "HCA Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("simMatching_GeneDetectRate_SS3-HCA-reads.pdf", height=4, width=8, useDingbats=F)
set.seed(111)
Genes = rownames(countsSim)
XX <- sample(Genes, 200)
X1 <- apply(countsSim[XX,], 1, function(x) sum(x!=0)) / dim(countsSim)[2]
X2 <- apply(compareData[XX,], 1, function(x) sum(x!=0)) / dim(compareData)[2]
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
legend('bottomright', c("Simulated", "HCA Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("simMatching_GeneMean_SS3-HCA-reads.pdf", height=4, width=8, useDingbats=F)

X1 <- log(apply(countsSim[XX,], 1, function(x) mean(x))+1)
X2 <- log(apply(compareData[XX,], 1, function(x) mean(x))+1)
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
legend('bottomright', c("Simulated", "HCA Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("simMatching_GeneSD_SS3-HCA-reads.pdf", height=4, width=8, useDingbats=F)
X1 <- log(apply(countsSim[XX,], 1, function(x) sd(x))+1)
X2 <- log(apply(compareData[XX,], 1, function(x) sd(x))+1)
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
legend('bottomright', c("Simulated", "HCA Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()







