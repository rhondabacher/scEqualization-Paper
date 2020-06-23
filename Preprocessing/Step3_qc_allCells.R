## QC using scater R package:
################
library(scater)

setwd("EqualizationManuscript")
load(file="RDATA/ALL_exprData_trimmed_genes.RData")

pd <- data.frame(CellType = rep(c(rep("EC", length(ec.cells)), rep("TB", length(tb.cells))), 4),
                 Experiment = rep(c("Unequalized", "EQ", "EQ-Vary", "EQ-Half"), 
                                  each=length(c(ec.cells, tb.cells))))

genez <- rownames(unEQ_EC)
exprdat <- SingleCellExperiment(assays = list(counts = data.matrix(round(cbind(unEQ_EC[genez,ec.cells], 
                                                           unEQ_TB[genez,tb.cells],
                                                           EQ[genez,c(ec.cells, tb.cells)],
                                                           EQ_Vary[genez,c(ec.cells, tb.cells)], 
                                                           EQ_75pcnt[genez,c(ec.cells, tb.cells)])))), colData = pd)

assay(exprdat, "logcounts") <- log(counts(exprdat) + 1)


## PCA of all data. 
# plotPCA(exprdat, colour_by = "Experiment", shape_by="CellType")

# plotPCA(exprdat, colour_by = "CellType", shape_by="Experiment")
# Separated by cell type and some minor separation by experiment.

## Do some QC checks
exprdat <- calculateQCMetrics(exprdat)

colnames(colData(exprdat))
colnames(rowData(exprdat))

# Plot log10 depth versus %counts in top 100 genes:
p1 <- plotColData(exprdat, x = "log10_total_counts", 
    y = "pct_counts_in_top_50_features", colour_by="Experiment", shape_by="CellType")
    
plot(p1)

pdf("PLOTS/qcPlot_scater_genes.pdf", height=3, width=4)
print(p1)
dev.off()
# Definitely want to remove outlier looking cells.


median(exprdat$log10_total_counts) + 2*sd(exprdat$log10_total_counts)
median(exprdat$log10_total_counts) - 2*sd(exprdat$log10_total_counts)
## use lower bound

mean(exprdat$pct_counts_in_top_50_features) + 2*sd(exprdat$pct_counts_in_top_50_features)
mean(exprdat$pct_counts_in_top_50_features) - 2*sd(exprdat$pct_counts_in_top_50_features)
## use upper bound

keep.total <- exprdat$log10_total_counts > 5.4
keep.n <- exprdat$pct_counts_in_top_50_features < 31
filtered <- exprdat[,keep.total & keep.n]
dim(filtered)

save(filtered, file="RDATA/QC-cells_usingScater_trimmed_genes.RData")


