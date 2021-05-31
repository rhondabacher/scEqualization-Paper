

load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(SingleCellExperiment)
counts <- counts(filtered)[,which(filtered$Experiment == "Unequalized" & filtered$CellType == "EC")]
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

library(scaffold)
library(Rcpp)
set.seed(6432)
RcppZiggurat::zsetseed(9921)

scaffoldParams <- estimateScaffoldParams(sce, numCells = c(50,50),
                                             equalizationAmount = 1,
                                             popHet = c(1,1), 
                                             usePops = list(fc_mean = c(0,2),
                                                            fc_sd = c(0, .4),
                                                            propGenes = c(0,.6)))

newsce <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce)

library(scater)
simcounts.sce <- newsce
simcounts.sce <- logNormCounts(simcounts.sce) 
simcounts.sce <- runPCA(simcounts.sce, name="PCA",
                      ncomponents=15)
plotReducedDim(simcounts.sce, dimred = "PCA", colour_by = 'cellPopulation')
simcounts.sce <- runTSNE(simcounts.sce, perplexity=10, 
                       dimred="PCA", n_dimred=10)
Q <- plotReducedDim(simcounts.sce, dimred = "TSNE", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)

pdf(file = "multipop_2pops50-50.pdf", height = 5, width = 5)
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")
dev.off()
 


set.seed(3222)
RcppZiggurat::zsetseed(21)

scaffoldParams <- estimateScaffoldParams(sce, numCells = c(97, 3),
                                             equalizationAmount = 1,
                                             popHet = c(1,1),
                                             usePops = list(fc_mean = c(0,2), 
                                                            fc_sd = c(0, .4),
                                                            propGenes = c(0,.6)))
newsce <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce)

library(scater)
simcounts.sce <- newsce
simcounts.sce <- logNormCounts(simcounts.sce) 
simcounts.sce <- runPCA(simcounts.sce, name="PCA",
                        ncomponents=15)
simcounts.sce <- runTSNE(simcounts.sce, perplexity=10, 
                         dimred="PCA", n_dimred=10)
Q <- plotReducedDim(simcounts.sce, dimred = "TSNE", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)

pdf(file = "multipop_Rarepops97-3.pdf", height = 5, width = 5)
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")
dev.off()


set.seed(3222)
RcppZiggurat::zsetseed(21)

scaffoldParams <- estimateScaffoldParams(sce, numCells = c(50,50,50),
                                             equalizationAmount = 1,
                                             popHet = c(1,1),
                                             usePops = list(fc_mean = c(0,2,1.5), 
                                                            fc_sd = c(0, .4,.4),
                                                            propGenes = c(0,.6,.4)))
newsce <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce)

library(scater)
simcounts.sce <- newsce
simcounts.sce <- logNormCounts(simcounts.sce) 
simcounts.sce <- runPCA(simcounts.sce, name="PCA",
                        ncomponents=15)
simcounts.sce <- runTSNE(simcounts.sce, perplexity=10, 
                         dimred="PCA", n_dimred=10)
Q <- plotReducedDim(simcounts.sce, dimred = "TSNE", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)

pdf(file = "multipop_pops3.pdf", height = 5, width = 5)
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")
dev.off()


