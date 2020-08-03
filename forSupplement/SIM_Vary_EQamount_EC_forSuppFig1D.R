# Code to generate simulated set that best matches the original dataset:
library(quantreg)
library(SCnorm)
library(parallel)

library(Splash)

setwd("/ufrc/rbacher/rbacher/HPC/EQUALIZATION_PAPER")
source("RCODE/functionsForCountDepthPlots.R")

load(file="RDATA/BestMatchingSim_LiFangEC2_trimmed.RData")

sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = originalCounts))

MedExpr <- log(apply(originalCounts, 1, function(x) median(x[x!=0])))
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

params <- estimateSplashParameters(sce, degree = 3,
                                   percentRange = -1)
alterParams <- params
alterParams@captureEfficiency = simdata@metadata$capEfficiency
pcntRange <- c(0,.125,.25,.50,-1)
trackMAD <- list()
repMAD <- c()

set.seed(3222)

for(j in 1:length(pcntRange)) {
  repMAD <- c()
  alterParams@percentRange <- pcntRange[j]
  for(i in 1:20) {
    alterSimdata <- simulateSplash(alterParams, sce, model = "p", 
                                   inputInitial = simdata@metadata$initialSimCounts)
    simcounts <- counts(alterSimdata)
    BiocParallel::register(BiocParallel::SerialParam())
    gslopes1 <- getSlopes(Data = simcounts)
    
    MedExpr <- log(apply(simcounts, 1, function(x) median(x[x!=0])))
    splitby <- sort(MedExpr[MedExpr >= log(2)])
    grps <- length(splitby) / 10
    sreg_orig2 <- split(splitby, ceiling(seq_along(splitby) / grps))
    
    estMAD <- modeCalc(gslopes1, sreg_orig2)
    
    repMAD <- c(repMAD, median(abs(estMAD[1:10] - 1)))

  }
  trackMAD[[j]] <- repMAD
  names(trackMAD)[j] <- pcntRange[j]
}
trackMAD

save(trackMAD, file="RDATA/simulateFromMatching_calcMAD_varyEQ_trimmed.RData")

