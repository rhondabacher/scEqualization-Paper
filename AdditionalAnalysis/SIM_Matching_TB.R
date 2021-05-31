
library(quantreg)
library(SCnorm)
library(parallel)
library(scaffold) 



load(file="~/OneDrive - University of Florida/EqualizationManuscript/RDATA/QC-cells_usingScater_trimmed_genes.RData")


library(SingleCellExperiment)
tbcounts <- counts(filtered)[,which(filtered$Experiment == "Unequalized" & filtered$CellType == "TB")]
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = tbcounts))

source("CODE/functionsForCountDepthPlots.R")
ModeStat_SquaredError <- list()

source("CODE/forRevisions/fasterSloperCalc_Helper.R")
gslopes0 <- fasterSlope(Data = tbcounts, ditherCounts=FALSE)
MedExpr <- log(apply(tbcounts, 1, function(x) median(x[x!=0])))
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))


set.seed(1364)
  RcppZiggurat::zsetseed(88)
  
  params <- estimateScaffoldParameters(sce, equalizationAmount = .67)
  simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce)
  
  simcounts <- counts(simdata)
  rownames(simcounts)<- rownames(tbcounts)
  gslopes1 <- fasterSlope(Data = data.matrix(simcounts), ditherCounts=FALSE)
  ModeStat_SquaredError<-list()
  
pdf(paste0("PLOTS/countDepthPlots_TB_Simulation_NOEQ.pdf"), height=8, width=8, useDingbats = F)
  par(mfrow=c(2,2), mar=c(5,5,2,1))
  
  i = 1
  makePlot(gslopes0, "unEQ TB", sreg_orig1)
  ModeStat_SquaredError[[i]] <- dotPlot(gslopes0, sreg_orig1, "MAD = 0.695")

  i = i + 1
  makePlot(gslopes1, "Simulated", sreg_orig1)
  ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_orig1, "MAD = 0.542")
  
  ModeStatistic_TB <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))
  
  print(paste0(ModeStatistic_TB))

  ## "0.694716242661448" "0.542074363992173"
dev.off()


ModeStatistic_TB_BestSim <- ModeStatistic_TB

originalCounts <- tbcounts
save(ModeStatistic_TB_BestSim, simdata, params, originalCounts,
 file="RDATA/Sim_LiFangTB_trimmed.RData")




##################

countsSim <- counts(simdata)
compareData <- originalCounts

library(yarrr)
pdf("PLOTS/sim_SeqDepth_TB.pdf", height=4, width=8, useDingbats=F)
X <- data.frame( Depth = colSums(countsSim) / 1000000, Species = "Simulated")
Y <- data.frame( Depth = colSums(compareData) / 1000000, Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "Sequencing Depth (millions)", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Sequencing Depth (millions)", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "TB Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("PLOTS/sim_CellDetectRate_TB.pdf", height=4, width=8, useDingbats=F)

X <- data.frame( Depth = colSums(countsSim!=0) / nrow(countsSim), Species = "Simulated")
Y <- data.frame( Depth = colSums(compareData!=0) / nrow(countsSim), Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species, ylim=c(0,1),
           data = longdata,
           xlab = "", 
           ylab = "Detection Rate", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Detection Rate", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('topleft', c("Simulated", "TB Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("PLOTS/sim_GeneDetectRate_TB.pdf", height=4, width=8, useDingbats=F)
set.seed(8276)
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
legend('topleft', c("Simulated", "TB Data"), col=c("brown1","cornflowerblue"), lwd=3)

dev.off()



pdf("PLOTS/sim_GeneMean_TB.pdf", height=4, width=8, useDingbats=F)

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
legend('bottomright', c("Simulated", "TB Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("PLOTS/sim_GeneSD_TB.pdf", height=4, width=8, useDingbats=F)
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
legend('bottomright', c("Simulated", "TB Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()




