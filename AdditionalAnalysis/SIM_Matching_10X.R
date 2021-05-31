# Sim matching 10X data


library(DropletUtils)

#import
data10x <- read10xCounts("DATA/10X_data/")

## qc
data10x.ed <- emptyDrops(data10x@assays@data$counts)
is.cell.1 <- data10x.ed$FDR <= 0.01
sum(is.cell.1, na.rm=TRUE)

data10x <- data10x[,which(is.cell.1==TRUE)]

library(ggplot2)
library(scater)
data10x <- addPerCellQC(data10x, 
                            subsets=list(Mito=grep("MT-", rownames(data10x))))

plotColData(data10x, x = "sum", y="detected") 
is.mito <- grep("^MT-", rowData(data10x)$Symbol)
pbmc.qc <- perCellQCMetrics(data10x, subsets=list(MT=is.mito))

plot(pbmc.qc$sum, pbmc.qc$subsets_MT_percent, log="x",
     xlab="Total count", ylab='Mitochondrial %')

keep.total.1 <- isOutlier(pbmc.qc$sum, nmads=3, type="both")
keep.total.2 <- isOutlier(pbmc.qc$detected, nmads=3, type="both")
keep.total.3 <- isOutlier(pbmc.qc$subsets_MT_percent, nmads=3, type="higher")

toKeep <- keep.total.1 + keep.total.2 + keep.total.3

data10x <- data10x[,which(toKeep == 0)]

save(data10x, file="RDATA/ready10Xdata.RData")

dim(data10x)

dim(data10x)
mean(colSums(as.matrix(counts(data10x)))) / 1e6
mean(colMeans(as.matrix(counts(data10x))!=0))

### 

library(scaffold)


load("RDATA/ready10Xdata.RData")
set.seed(554488)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(counts(data10x))))
params <- estimateScaffoldParameters(sce, sceUMI = TRUE, protocol = "10X", useUMI=TRUE) 

simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce)

save.image("RDATA/Sim_10Xdata.RData")

countsSim <- simdata@assays@data$umi_counts
compareData <- as.matrix(counts(data10x))


library(yarrr)


pdf("PLOTS/simMatching_SeqDepth_10X.pdf", height=4, width=8, useDingbats=F)
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
legend('bottomright', c("Simulated", "10X Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("PLOTS/simMatching_CellDetectRate_10X.pdf", height=4, width=8, useDingbats=F)
X <- data.frame( Depth = colSums(countsSim!=0) / nrow(countsSim), Species = "Simulated")
Y <- data.frame( Depth = colSums(compareData!=0) / nrow(countsSim), Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species, 
           data = longdata, ylim=c(0,.2),
           xlab = "", 
           ylab = "Detection Rate", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Detection Rate", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "10X Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("PLOTS/simMatching_GeneDetectRate_10X.pdf", height=4, width=8, useDingbats=F)

set.seed(6655)
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
legend('bottomright', c("Simulated", "10X Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("PLOTS/simMatching_GeneMean_10X.pdf", height=4, width=8, useDingbats=F)

X1 <- (apply(countsSim[XX,], 1, function(x) mean(x)))
X2 <- (apply(compareData[XX,], 1, function(x) mean(x)))
useg <- names(which(X1<Inf & X2 < Inf))
X <- data.frame( Depth =X1[useg], Species = "Simulated")
Y <- data.frame( Depth = X2[useg], Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "mean", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="mean", 
     cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "10X Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("PLOTS/simMatching_GeneSD_10X.pdf", height=4, width=8, useDingbats=F)
X1 <- (apply(countsSim[XX,], 1, function(x) sd(x)))
X2 <- (apply(compareData[XX,], 1, function(x) sd(x)))
useg <- names(which(X1<Inf & X2 < Inf & X1>-Inf &X2>-Inf))
X <- data.frame( Depth =X1[useg], Species = "Simulated")
Y <- data.frame( Depth = X2[useg], Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "sd", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="sd", 
     cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "10X Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


