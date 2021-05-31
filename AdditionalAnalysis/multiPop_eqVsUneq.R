
load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(SingleCellExperiment)
counts <- counts(filtered)[,which(filtered$Experiment == "Unequalized" & filtered$CellType == "EC")]
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))


library(scaffold)
library(Rcpp)
library(scater)
library(cluster)

set.seed(2424)
RcppZiggurat::zsetseed(99)

sil_eq <- list()
sil_uneq <- list()
sil_diff <- list()
for (k in 1:25) {
scaffoldParams <- estimateScaffoldParameters(sce, numCells = c(50,40),
                                             equalizationAmount = 1,
                                             popHet = c(.5,2),
                                             usePops = list(fc_mean = c(0,1.25), 
                                                             fc_sd=c(0,.4),
                                                            propGenes = c(0,.15)))
newsce <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce)

simcounts.sce <- newsce
simcounts.sce <- logNormCounts(simcounts.sce) 
simcounts.sce <- runPCA(simcounts.sce, name="PCA",
                        ncomponents=15)
# plotReducedDim(simcounts.sce, dimred = "PCA", colour_by = 'cellPopulation')
simcounts.sce <- runTSNE(simcounts.sce, perplexity=10,
                         dimred="PCA", n_dimred=10)
Q <- plotReducedDim(simcounts.sce, dimred = "TSNE", colour_by = "cellPopulation")
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")

scaffoldParams@captureEfficiency <- colData(newsce)$capEfficiency
scaffoldParams@equalizationAmount <- 0
newsce.eq <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce, inputInitial = newsce@metadata$initialSimCounts)

simcounts.sce.eq <- newsce.eq
simcounts.sce.eq <- logNormCounts(simcounts.sce.eq) 
simcounts.sce.eq <- runPCA(simcounts.sce.eq, name="PCA",
                        ncomponents=15)
# plotReducedDim(simcounts.sce.eq, dimred = "PCA", colour_by = 'cellPopulation')
simcounts.sce.eq <- runTSNE(simcounts.sce.eq, perplexity=10,
                         dimred="PCA", n_dimred=10)
Q <- plotReducedDim(simcounts.sce.eq, dimred = "TSNE", colour_by = "cellPopulation")
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")


uneq.dist <- dist(simcounts.sce@int_colData$reducedDims$TSNE[,1:2])
uneq.sil <- silhouette(x=as.numeric( colData(simcounts.sce)$cellPopulation), dist=uneq.dist)
sil_uneq[[k]] <- summary(uneq.sil)$avg.width
print(summary(uneq.sil)$avg.width)

eq.dist <- dist(simcounts.sce.eq@int_colData$reducedDims$TSNE[,1:2])
eq.sil <- silhouette(x=as.numeric( colData(simcounts.sce.eq)$cellPopulation), dist=eq.dist)
sil_eq[[k]] <- summary(eq.sil)$avg.width
print(summary(eq.sil)$avg.width)

sil_diff[[k]] <- sil_eq[[k]] - sil_uneq[[k]]
}

print(unlist(sil_diff))
mean(unlist(sil_diff)) #0.07
median(unlist(sil_diff)) #0.11
boxplot(unlist(sil_uneq), unlist(sil_eq))
boxplot(unlist(sil_diff))

plot(unlist(sil_uneq), unlist(sil_eq))

save.image(file="RDATA/multiPop_SilDistanceCompv2.RData")


load(file="RDATA/multiPop_SilDistanceCompv2.RData")


pdf("PLOTS/comparisonMultiPop.pdf", height=5, width=5)
par(mfrow=c(1,1), mar=c(3,5,2,1))
boxplot(unlist(sil_uneq), unlist(sil_eq), ylab="Average Silhouette Width", cex.axis=1.3,cex.lab=1.5,
        main="",
        names=c("unEQ", "EQ"))
dev.off()

# Example:

set.seed(390)
RcppZiggurat::zsetseed(31)
scaffoldParams <- estimateScaffoldParameters(sce, numCells = c(50,40),
                                             percentRange = -1,
                                             degree = c(.5,2),
                                             usePops = list(fc = c(0,1.25), 
                                                            propGenes = c(0,.15)))
newsce <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce)

simcounts.sce <- newsce
simcounts.sce <- logNormCounts(simcounts.sce) 
simcounts.sce <- runPCA(simcounts.sce, name="PCA",
                        ncomponents=15)
# plotReducedDim(simcounts.sce, dimred = "PCA", colour_by = 'cellPopulation')
simcounts.sce <- runTSNE(simcounts.sce, perplexity=10,
                         dimred="PCA", n_dimred=10)
Q <- plotReducedDim(simcounts.sce, dimred = "TSNE", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")

scaffoldParams@captureEfficiency <- colData(newsce)$capEfficiency
scaffoldParams@percentRange <- 0
newsce.eq <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce, inputInitial = newsce@metadata$initialSimCounts)

simcounts.sce.eq <- newsce.eq
simcounts.sce.eq <- logNormCounts(simcounts.sce.eq) 
simcounts.sce.eq <- runPCA(simcounts.sce.eq, name="PCA",
                           ncomponents=15)
simcounts.sce.eq <- runTSNE(simcounts.sce.eq, perplexity=10,
                            dimred="PCA", n_dimred=10)
Q <- plotReducedDim(simcounts.sce.eq, dimred = "TSNE", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")


uneq.dist <- dist(simcounts.sce@int_colData$reducedDims$TSNE[,1:2])
uneq.sil <- silhouette(x=as.numeric( colData(simcounts.sce)$cellPopulation), dist=uneq.dist)
print(summary(uneq.sil)$avg.width)

eq.dist <- dist(simcounts.sce.eq@int_colData$reducedDims$TSNE[,1:2])
eq.sil <- silhouette(x=as.numeric( colData(simcounts.sce.eq)$cellPopulation), dist=eq.dist)
print(summary(eq.sil)$avg.width)


pdf("PLOTS/comparisonMultiPop_unEQ.pdf", height=5, width=5)
Q <- plotReducedDim(simcounts.sce, dimred = "TSNE", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")
dev.off()


pdf("PLOTS/comparisonMultiPop_EQ.pdf", height=5, width=5)
Q <- plotReducedDim(simcounts.sce.eq, dimred = "TSNE", colour_by = "cellPopulation",
                    point_alpha = .8, point_size=4)
Q + theme(text = element_text(size=20)) + scale_colour_discrete(name="Population")
dev.off()





