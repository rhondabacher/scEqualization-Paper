# Code to generate simulated set that best matches the original dataset:
library(Splash)

load(file="RDATA/BestMatchingSim_LiFangEC2_trimmed.RData")

sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = originalCounts))

set.seed(58893)
propGenesPerCell <- c()
propGenesDetectedALL <- c()
propGenesDetectedNONE <- c()
numGenesDetectedALL <- c()
meanDEV <- c()
meanMEAN <- c()
meanDEPTH <- c()
meanCLD <- c()
meanGLD <- c()
for(i in 1:25) {
  params <- estimateSplashParameters(sce, 
                                     degree = 3, percentRange = -1)
  mainSimdata <- simulateSplash(splashParams = params, sce, model = "p", SD=.02)
  INIT <- mainSimdata@metadata$initialSimCounts
  params@captureEfficiency = mainSimdata@metadata$capEfficiency 
  
  params@percentRange = -1
  uneqSimdata <- simulateSplash(splashParams = params, sce, model = "p", 
                                inputInitial = INIT)
  simcountsUN <- round(counts(uneqSimdata))
  
  params@percentRange <- 0
  eqSimdata <- simulateSplash(params, sce, model = "p", 
                              inputInitial = INIT)
  simcountsEQ <- round(counts(eqSimdata))
  
  ## Check Zeros:
  
  propGenesPerCell[i] <- (mean(colSums(simcountsEQ!=0)) - mean(colSums(simcountsUN!=0))) / mean(colSums(simcountsUN!=0)) 
  
  zeroNEW <- apply(simcountsEQ, 1, function(x) mean(x!=0))
  zeroOLD <- apply(simcountsUN, 1, function(x) mean(x!=0))
  
  propGenesDetectedALL[i] <- (mean(zeroNEW == 1) - mean(zeroOLD == 1)) / mean(zeroOLD == 1)
  
  propGenesDetectedNONE[i] <- (mean(zeroNEW == 0) - mean(zeroOLD == 0)) / mean(zeroOLD == 0)
  
  numGenesDetectedALL[i] <- (sum(zeroNEW == 1) - sum(zeroOLD == 1))
  
  X1 <- (apply(simcountsEQ, 1, function(x) sd(x)))
  X2 <- (apply(simcountsUN, 1, function(x) sd(x)))
  prop1 <- ((X1 - X2) / X2)
  prop1 <- prop1[which(X2 > 0)]
  meanDEV[i] <- median(prop1, na.rm=T)
  
  X1 <- (apply(simcountsEQ, 1, function(x) mean(x)))
  X2 <- (apply(simcountsUN, 1, function(x) mean(x)))
  prop1 <- ((X1 - X2) / X2)
  prop1 <- prop1[which(X2 > 0)]
  meanMEAN[i] <- median(prop1, na.rm=T)
  
  X1 <- (apply(simcountsEQ, 2, function(x) sum(x!=0)))
  X2 <- (apply(simcountsUN, 2, function(x) sum(x!=0)))
  prop1 <- ((X1 - X2) / X2)
  prop1 <- prop1[which(X2 > 0)]
  meanCLD[i] <- median(prop1, na.rm=T)
  
  
  X1 <- (apply(simcountsEQ, 1, function(x) sum(x!=0)))
  X2 <- (apply(simcountsUN, 1, function(x) sum(x!=0)))
  prop1 <- ((X1 - X2) / X2)
  prop1 <- prop1[which(X2 > 0)]
  meanGLD[i] <- median(prop1, na.rm=T)
  
  
  X1 <- colSums(simcountsEQ) / 1000000
  X2 <- colSums(simcountsUN) / 1000000
  meanDEPTH[i] <- (sd(X1) - sd(X2)) / sd(X2)
}


mean(propGenesPerCell)
mean(propGenesDetectedALL)
mean(propGenesDetectedNONE)
mean(numGenesDetectedALL)

mean(meanDEPTH)
mean(meanGLD)
mean(meanCLD)
mean(meanMEAN)
mean(meanDEV)


save.image("RDATA/compare_SIM-properties_unEQvsEQ_v2.RData")


