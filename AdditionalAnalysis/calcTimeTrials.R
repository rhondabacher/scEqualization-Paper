library(scaffold) 


load(file="~/OneDrive - University of Florida/EqualizationManuscript/RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(SingleCellExperiment)
tbcounts <- counts(filtered)[,which(filtered$Experiment == "Unequalized" & filtered$CellType == "TB")] 
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = tbcounts))

set.seed(1364)
RcppZiggurat::zsetseed(88)

getTIME1 <- list()
for(k in 1:10) {
  getTIME1[[k]] <- system.time(simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce))
  print(k)
}
save(getTIME1, file="timeTrials-nonUMI.RData")


params <- estimateScaffoldParameters(sce, numCells = 1000)
getTIME2 <- list()
for(k in 1:10) {
  getTIME2[[k]] <- system.time(simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce))
}

save(getTIME1,getTIME2, file="timeTrials-nonUMI.RData")

mean(unlist(sapply(getTIME1, function(x) x[3]))) / 60
mean(unlist(sapply(getTIME2, function(x) x[3]))) / 60

params <- estimateScaffoldParameters(sce, numCells = 5000)
getTIME3 <- list()
for(k in 1:10) {
  getTIME3[[k]] <- system.time(simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce))
}

save(getTIME1,getTIME2,getTIME3, file="timeTrials-nonUMI.RData")

load("RDATA/timeTrials-nonUMI.RData")
mean(unlist(sapply(getTIME1, function(x) x[3]))) / 60 # 0.14
mean(unlist(sapply(getTIME2, function(x) x[3]))) / 60 # 1.52
mean(unlist(sapply(getTIME3, function(x) x[3]))) / 60 # 9.41

##############################################################################################################3

load("RDATA/forRevision/Sim_10Xdata.RData")
set.seed(99119)
RcppZiggurat::zsetseed(99)
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data.matrix(counts(data10x))))

params <- estimateScaffoldParameters(sce, useUMI = TRUE, fromUMI=TRUE, numCells=1000)
getTIME1 <- list()
for(k in 6:10) {
  getTIME1[[k]] <- system.time(simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce))
  print(k)
}
save(getTIME1, file="timeTrials-10X.RData")


params <- estimateScaffoldParameters(sce, useUMI = TRUE, fromUMI=TRUE, numCells=5000)
getTIME2 <- list()
for(k in 1:10) {
  getTIME2[[k]] <- system.time(simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce))
}

save(getTIME1,getTIME2, file="timeTrials-10X.RData")

params <- estimateScaffoldParameters(sce, useUMI = TRUE, fromUMI=TRUE, numCells=10000)
getTIME3 <- list()
for(k in 1:10) {
  getTIME3[[k]] <- system.time(simdata <- simulateScaffold(scaffoldParams=params, originalSCE=sce))
  print(k)
}

save(getTIME1,getTIME2,getTIME3, file="timeTrials-10X.RData")


load("timeTrials-10X.RData")

mean(unlist(sapply(getTIME1, function(x) x[3]))) / 60
mean(unlist(sapply(getTIME2, function(x) x[3]))) / 60
mean(unlist(sapply(getTIME3, function(x) x[3]))) / 60

