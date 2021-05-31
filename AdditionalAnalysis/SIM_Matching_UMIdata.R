## UMI data: 
##### GSE81076-GPL16791.rds from conquer


library(SingleCellExperiment)
library(MultiAssayExperiment)
umidata <- readRDS("DATA/GSE81076-GPL16791.rds")
umidata = experiments(umidata)[[1]] 
umidata <- assay(umidata, "count") 
# 
# 
# ## Need to do some QC on the UMI data:
library(scater)
umidata.sce <- SingleCellExperiment(assay=list(counts = umidata))
example_sce <- addPerCellQC(umidata.sce)
plotColData(example_sce, x = "total", y="detected")

keep.total.1 <- isOutlier(example_sce$detected, nmads=5, type="higher")
umidata <- umidata[,!keep.total.1]
# 
save(umidata, file="RDATA/dataReady_UMI.RData")

dim(umidata)
mean(colSums(umidata)) / 1e6
mean(colMeans(umidata!=0))
### 

library(scaffold)

load("RDATA/dataReady_UMI.RData")
set.seed(99119)
RcppZiggurat::zsetseed(99)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = umidata))
params <- estimateScaffoldParameters(sce, equalizationAmount = 1, useUMI = TRUE, sceUMI=TRUE)
simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce)

save.image("RDATA/Sim_UMIdata.RData")

countsSim <- simdata@assays@data$umi_counts 
compareData <- umidata


dim(countsSim)
dim(compareData)
library(SCnorm)

source("CODE/forRevisions/fasterSloperCalc_Helper.R")
gslopes0 <- fasterSlope(compareData, ditherCounts=TRUE)
gslopes1 <- fasterSlope(countsSim, ditherCounts = TRUE)


save.image("RDATA/Sim_UMIdata.RData")

source("CODE/functionsForCountDepthPlots.R")
ModeStat_SquaredError<-list()


pdf("PLOTS/simMatching_CountDepth_UMI.pdf", height=8, width=8, useDingbats=F)
MedExpr <- log(apply(umidata, 1, function(x) mean(x)))
splitby <- sort(MedExpr[MedExpr > 0])
grps <- length(splitby) / 10
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

par(mfrow=c(2,2), mar=c(5,5,2,1))

i = 1
makePlot(gslopes0, "UMI", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes0, sreg_orig1, "UMI")

i = i + 1
makePlot(gslopes1, "Simulated", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_orig1, "Simulated")

ModeStatistic_EC <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))

print(paste0(ModeStatistic_EC))
## "0.307240704500979" "0.512720156555773"
dev.off()




library(yarrr)


pdf("PLOTS/simMatching_SeqDepth_UMI.pdf", height=4, width=8, useDingbats=F)
X <- data.frame( Depth = colSums(countsSim)/1e3 , Species = "Simulated")
Y <- data.frame( Depth = colSums(compareData)/1e3 , Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata, 
           xlab = "", 
           ylab = "Sequencing Depth (thousands)", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Sequencing Depth (thousands)", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "UMI Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("PLOTS/simMatching_CellDetectRate_UMI.pdf", height=4, width=8, useDingbats=F)
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
legend('bottomright', c("Simulated", "UMI Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("PLOTS/simMatching_GeneDetectRate_UMI.pdf", height=4, width=8, useDingbats=F)

set.seed(545)
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
legend('bottomright', c("Simulated", "UMI Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("PLOTS/simMatching_GeneMean_UMI.pdf", height=4, width=8, useDingbats=F)

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
legend('bottomright', c("Simulated", "UMI Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("PLOTS/simMatching_GeneSD_UMI.pdf", height=4, width=8, useDingbats=F)
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
legend('bottomright', c("Simulated", "UMI Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()
