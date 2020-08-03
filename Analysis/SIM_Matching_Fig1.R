library(devtools)
# install("~/Desktop/SOFTWARE_DEVEL/scaffold/")
install("~/Desktop/SOFTWARE_DEVEL/Splash")

# Code to generate simulated set that best matches the original dataset:
library(splatter)
library(quantreg)
library(SCnorm)
library(parallel)
library(Splash) ## Change to scaffold LATER!

setwd("~/OneDrive - University of Florida/EqualizationManuscript")

## Don't want the data normalized yet.
load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

counts <- counts(filtered)[,which(filtered$Experiment == "Unequalized")]

counts <- counts[,grepl("EC",colnames(counts))] # Only the EC cells!
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

source("CODE/functionsForCountDepthPlots.R")
ModeStat_SquaredError <- list()

gslopes0 <- getSlopes(Data = counts)
MedExpr <- log(apply(counts, 1, function(x) median(x[x!=0])))
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))


set.seed(4869)
params <- estimateSplashParameters(sce, degree = 3, 
                                   percentRange = -1)
simdata <- simulateSplash(splashParams=params, originalSCE=sce, model = "p")

simcounts <- counts(simdata)
BiocParallel::register(BiocParallel::SerialParam())
gslopes1 <- getSlopes(Data = simcounts)
ModeStat_SquaredError<-list()

pdf(paste0("PLOTS/countDepthPlots_EC_MatchingSimulation_NOEQ.pdf"), height=8, width=8, useDingbats = F)
par(mfrow=c(2,2), mar=c(5,5,2,1))

i = 1
makePlot(gslopes0, "EC2 - Original", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes0, sreg_orig1, "EC2 - Original")
MedExpr <- log(apply(simcounts, 1, function(x) median(x[x!=0])))
splitby <- sort(MedExpr[MedExpr >= log(2)])
grps <- length(splitby) / 10
sreg_orig2 <- split(splitby, ceiling(seq_along(splitby) / grps))

i = i + 1
makePlot(gslopes1, "EC2 - Sim", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_orig1, "EC2 - Sim")

ModeStatistic_EC <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))
print(paste0(ModeStatistic_EC))

dev.off()



ModeStatistic_EC <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))
print(ModeStatistic_EC)

ModeStatistic_EC_BestSim <- ModeStatistic_EC

originalCounts <- counts
save(ModeStatistic_EC_BestSim, simdata, params, originalCounts, 
     file="RDATA/BestMatchingSim_LiFangEC2_trimmed.RData")


