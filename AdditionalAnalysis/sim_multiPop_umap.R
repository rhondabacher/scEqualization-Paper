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

set.seed(1515)
RcppZiggurat::zsetseed(1221)

sil_eq_euc <- list()
sil_uneq_euc <- list()
sil_diff_euc <- list()
sil_sd_euc <- list()
sil_nrm_euc <- list()

sil_eq_euc_all <- list()
sil_uneq_euc_all <- list()
sil_diff_euc_all <- list()
sil_sd_euc_all <- list()
sil_nrm_euc_all <- list()


for (j in 1:250) {

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
  simcounts.sce <- runPCA(simcounts.sce, name="PCA", scale=TRUE, #subset_row=gvar.uneq,
                          ncomponents=15)
       
  simcounts.sce.nrm <- newsce
  simcounts.sce.nrm <- computeSumFactors(simcounts.sce.nrm)
  simcounts.sce.nrm <- logNormCounts(simcounts.sce.nrm)  
  simcounts.sce.nrm <- runPCA(simcounts.sce.nrm, name="PCA", scale=TRUE, 
                          ncomponents=15)
   
  simcounts.sce.eq <- newsce.eq
  simcounts.sce.eq <- logNormCounts(simcounts.sce.eq, size.factors=rep(1,sum(scaffoldParams@numCells)))  
  simcounts.sce.eq <- runPCA(simcounts.sce.eq, name="PCA",scale=TRUE, #subset_row=gvar.eq,
                             ncomponents=15)

  simcounts.sce.sd <- newsce.eq
  simcounts.sce.sd <- computeSumFactors(simcounts.sce.sd)
  simcounts.sce.sd <- logNormCounts(simcounts.sce.sd)
  simcounts.sce.sd <- runPCA(simcounts.sce.sd, name="PCA",scale=TRUE, 
                             ncomponents=15)
	
	 for(k in 1:25) {

	     simcounts.sce <- runUMAP(simcounts.sce, name="UMAP", n_neighbors=10,
	                              dimred="PCA", n_dimred=10)
	     simcounts.sce.nrm <- runUMAP(simcounts.sce.nrm, name="UMAP", n_neighbors=10,
	                                  dimred="PCA", n_dimred=10)
	    simcounts.sce.eq <- runUMAP(simcounts.sce.eq, name="UMAP", n_neighbors=10,
	                                 dimred="PCA", n_dimred=10)
	     simcounts.sce.sd <- runUMAP(simcounts.sce.sd, name="UMAP", n_neighbors=10,
	                                 dimred="PCA", n_dimred=10)

    uneq.dist <- dist((simcounts.sce@int_colData$reducedDims$UMAP[,1:2]), method = "e")
    uneq.sil <- silhouette(x=as.numeric( colData(simcounts.sce)$cellPopulation), dist=uneq.dist)
    sil_uneq_euc[[k]] <- summary(uneq.sil)$avg.width
  
    nrm.dist <- dist((simcounts.sce.nrm@int_colData$reducedDims$UMAP[,1:2]), method = "e")
    nrm.sil <- silhouette(x=as.numeric( colData(simcounts.sce.nrm)$cellPopulation), dist=nrm.dist)
    sil_nrm_euc[[k]] <- summary(nrm.sil)$avg.width
    
    eq.dist <- dist((simcounts.sce.eq@int_colData$reducedDims$UMAP[,1:2]), method = "e")
    eq.sil <- silhouette(x=as.numeric( colData(simcounts.sce.eq)$cellPopulation), dist=eq.dist)
    sil_eq_euc[[k]] <- summary(eq.sil)$avg.width
    
    sd.dist <- dist((simcounts.sce.sd@int_colData$reducedDims$UMAP[,1:2]), method = "e")
    sd.sil <- silhouette(x=as.numeric( colData(simcounts.sce.sd)$cellPopulation), dist=sd.dist)
    sil_sd_euc[[k]] <- summary(sd.sil)$avg.width
		
  
	}
	
	sil_uneq_euc_all[[j]] <- sil_uneq_euc
	sil_nrm_euc_all[[j]] <- sil_nrm_euc
	sil_eq_euc_all[[j]] <- sil_eq_euc
	sil_sd_euc_all[[j]] <- sil_sd_euc
    
}

save.image("multiPop_umap.RData")
 