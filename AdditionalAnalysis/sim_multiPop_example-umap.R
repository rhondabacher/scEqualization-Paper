load(file="QC-cells_usingScater_trimmed_genes.RData")

library(SingleCellExperiment)
counts <- counts(filtered)[,which(filtered$Experiment == "Unequalized" & filtered$CellType == "EC")]
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

library(scaffold)
library(Rcpp)
library(scater)
library(cluster)
require(gridExtra)
library(scran)

set.seed(8112)
RcppZiggurat::zsetseed(9898)

scaffoldParams <- estimateScaffoldParameters(sce, numCells = c(50,40),
                                             usePops = list(fc_mean = c(0,1.5), 
                                                            fc_sd = c(0, .5),
                                                            propGenes = c(0,.1)))

scaffoldParams@equalizationAmount <- 1
newsce <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce)

scaffoldParams@equalizationAmount <- 0
newsce.eq <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce, 
                              inputInitial = newsce@metadata$initialSimCounts,
                              amplifiedMolecules = newsce@metadata$amplifiedMolecules)

simcounts.sce <- newsce
simcounts.sce <- logNormCounts(simcounts.sce, size.factors=rep(1,sum(scaffoldParams@numCells)))  
simcounts.sce <- runPCA(simcounts.sce, name="PCA", scale=TRUE, 
                        ncomponents=15)


simcounts.sce.nrm <- newsce
simcounts.sce.nrm <- computeSumFactors(simcounts.sce.nrm)
simcounts.sce.nrm <- logNormCounts(simcounts.sce.nrm)  
simcounts.sce.nrm <- runPCA(simcounts.sce.nrm, name="PCA", scale=TRUE, 
                            ncomponents=15)


simcounts.sce.eq <- newsce.eq
simcounts.sce.eq <- logNormCounts(simcounts.sce.eq, size.factors=rep(1,sum(scaffoldParams@numCells)))  
simcounts.sce.eq <- runPCA(simcounts.sce.eq, name="PCA",scale=TRUE, 
                           ncomponents=15)

simcounts.sce.sd <- newsce.eq
simcounts.sce.sd <- computeSumFactors(simcounts.sce.sd)
simcounts.sce.sd <- logNormCounts(simcounts.sce.sd)
simcounts.sce.sd <- runPCA(simcounts.sce.sd, name="PCA",scale=TRUE, 
                           ncomponents=15)

simcounts.sce <- runUMAP(simcounts.sce, name="UMAP", n_neighbors=10,
                         dimred="PCA", n_dimred=10)
simcounts.sce.nrm <- runUMAP(simcounts.sce.nrm, name="UMAP", n_neighbors=10,
                             dimred="PCA", n_dimred=10)
simcounts.sce.eq <- runUMAP(simcounts.sce.eq, name="UMAP", n_neighbors=10,
                            dimred="PCA", n_dimred=10)
simcounts.sce.sd <- runUMAP(simcounts.sce.sd, name="UMAP", n_neighbors=10,
                            dimred="PCA", n_dimred=10)

Q <- plotReducedDim(simcounts.sce, dimred = "UMAP", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)
p1 <- Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")

Q <- plotReducedDim(simcounts.sce.nrm, dimred = "UMAP", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)
p2 <- Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")

Q <- plotReducedDim(simcounts.sce.eq, dimred = "UMAP", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)
p3 <- Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")

Q <- plotReducedDim(simcounts.sce.sd, dimred = "UMAP", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)
p4 <- Q + theme(text = element_text(size=20)) + scale_colour_discrete(name=i)


grid.arrange(p1,p2,p3,p4, nrow=2)


pdf("umap_example.pdf", height=6, width=4)
grid.arrange(p2,p4, nrow=2)
dev.off()

