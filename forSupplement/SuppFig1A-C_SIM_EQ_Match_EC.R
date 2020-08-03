# Code to generate simulated set that best matches the original dataset:
library(quantreg)
library(SCnorm)
library(parallel)

library(devtools)
# install("~/Desktop/SOFTWARE_DEVEL/Splash")
library(Splash)

setwd("~/OneDrive - University of Florida/EqualizationManuscript")
source("CODE/functionsForCountDepthPlots.R")

load(file="RDATA/BestMatchingSim_LiFangEC2_trimmed.RData")

sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = originalCounts))

gslopes0 <- getSlopes(Data = originalCounts)
MedExpr <- log(apply(originalCounts, 1, function(x) median(x[x!=0])))
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

set.seed(2479)
params <- estimateSplashParameters(sce, 
                              captureEfficiency = simdata@metadata$capEfficiency, 
                                   degree = 3, percentRange = 0)
alterSimdata <- simulateSplash(params, sce, model = "p", 
                            inputInitial = simdata@metadata$initialSimCounts)

simcounts.eq <- counts(alterSimdata)
BiocParallel::register(BiocParallel::SerialParam())
gslopes1 <- getSlopes(Data = simcounts.eq)
ModeStat_SquaredError<-list()

pdf("PLOTS/countDepthPlots_EC_MatchingSimulation_EQ.pdf", height=8, width=4, useDingbats = F)
par(mfrow=c(2,1), mar=c(5,5,2,1))
i = 1
makePlot(gslopes1, "EC2 - Sim", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_orig1, "EC2 - Sim")
ModeStatistic_EC <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))
print(ModeStatistic_EC)
dev.off()

## Get statistics for paper:
# 

simcounts <- counts(simdata)
(mean(colSums(simcounts.eq!=0)) - mean(colSums(simcounts!=0))) /
  mean(colSums(simcounts!=0))

zeroNEW <- apply(simcounts.eq, 1, function(x) mean(x!=0))
zeroOLD <- apply(simcounts, 1, function(x) mean(x!=0))

(sum(zeroNEW == 1) - sum(zeroOLD == 1)) / sum(zeroOLD == 1)

(sum(zeroNEW == 0) - sum(zeroOLD == 0)) / sum(zeroOLD == 0)
