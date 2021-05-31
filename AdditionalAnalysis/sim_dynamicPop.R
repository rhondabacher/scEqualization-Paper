

load(file="RDATA/QC-cells_usingScater_trimmed_genes.RData")

library(SingleCellExperiment)
eccounts <- counts(filtered)[,which(filtered$Experiment == "Unequalized" & filtered$CellType == "EC")]
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = eccounts))

library(scaffold)
library(Rcpp)

library(ggplot2)


set.seed(2123)
RcppZiggurat::zsetseed(772)

# Simulate non-dynamic population as a comparison:
scaffoldParams0 <- estimateScaffoldParams(sce, numCells = c(100),
                                             equalizationAmount = 1,
                                             popHet = c(1,1))
newsce0 <- simulateScaffold(scaffoldParams=scaffoldParams0, originalSCE=sce)

## Show trajectory of the data:
library(SCORPIUS)
space0 <- reduce_dimensionality(t(counts(newsce0)), "pearson")
traj0 <- infer_trajectory(space0)
pdf("trajNonDyn.pdf", height=4, width=4)
draw_trajectory_plot(space0, contour = F)+ theme(text = element_text(size=20))
dev.off()

# Simulate dynamic population 20% of genes dynamic:
scaffoldParams <- estimateScaffoldParameters(sce, numCells = c(100),
                                             equalizationAmount = 1,
                                             popHet = c(1,1), 
                                             useDynamic = list(propGenes=.2))
newsce <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce)

## Show trajectory of the data:
library(SCORPIUS)
space <- reduce_dimensionality(t(counts(newsce)), "pearson")
traj <- infer_trajectory(space)
uneq.time <- rev(traj$time)
pdf("trajDynUnEQ.pdf", height=4, width=4)
draw_trajectory_plot(space, path=traj$path, contour = F) + theme(text = element_text(size=20))
dev.off()

scaffoldParams@captureEfficiency <- newsce$capEfficiency
scaffoldParams@percentRange <- 0
newsce.eq <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce, 
                               inputInitial = metadata(newsce)[[1]])
space <- reduce_dimensionality(t(counts(newsce.eq)), "pearson")
traj <- infer_trajectory(space)
eq.time <- (traj$time)
pdf("PLOTS/forRevision/trajDynEQ.pdf", height=4, width=4)
draw_trajectory_plot(space, path=traj$path, contour = F)+ theme(text = element_text(size=20))
dev.off()

pdf("trajDyn_compareTime.pdf", height=4, width=4)
plot(uneq.time, eq.time, xlab="UnEQ Pseudotime", ylab="EQ Pseudotime")
dev.off()



set.seed(6733)
RcppZiggurat::zsetseed(3211)

Rescale <- function(Data){
  InEC=Data
  ECNo0=InEC[which(rowMeans(InEC)>0),]
  ECSC=t(apply(ECNo0,1,scale))
  rownames(ECSC)=rownames(ECNo0)
  colnames(ECSC)=colnames(ECNo0)
  ECSC
}

auc.eq <- c()
auc.uneq <- c()
library(data.table)
for(i in 1:25) {
  scaffoldParams <- estimateScaffoldParameters(sce, numCells = c(100),
                                          equalizationAmount = 1,
                                          popHet = c(1,1), 
                                          useDynamic = list(propGenes=.2))
  newsce <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce)
  scaffoldParams@captureEfficiency <- newsce$capEfficiency
  scaffoldParams@percentRange <- 0
  newsce.eq <- simulateScaffold(scaffoldParams=scaffoldParams, originalSCE=sce, 
                                inputInitial = metadata(newsce)[[1]])
  
  library(SCORPIUS)
  space <- reduce_dimensionality(t(counts(newsce)), "pearson")
  traj <- infer_trajectory(space)
  uneq.time <- rev(traj$time)
  
  space <- reduce_dimensionality(t(counts(newsce.eq)), "pearson")
  traj <- infer_trajectory(space)
  eq.time <- (traj$time)
  
  if(cor(uneq.time, eq.time) < 0) {eq.time <- rev(traj$time)} 
  
  data.use <- counts(newsce.eq)
  data.sc <- Rescale(data.use[,order(eq.time)])
  fitall.all.eq <- sapply(1:nrow(data.sc),function(j){
    print(j)
    tt=lm(data.sc[j,]~poly(1:ncol(data.sc),2))
    t2 <- anova(tt)$'Pr(>F)'[1]
    t2  })
  names(fitall.all.eq) <- rownames(data.sc)
  fitall.all.eq.adj <- p.adjust(fitall.all.eq, 'fdr')

  data.use <- counts(newsce)
  data.sc <- Rescale(data.use[,order(uneq.time)])
  fitall.all <- sapply(1:nrow(data.sc),function(j){
    print(j)
    tt=lm(data.sc[j,]~poly(1:ncol(data.sc),2))
    t2 <- anova(tt)$'Pr(>F)'[1]
    t2  })
  names(fitall.all) <- rownames(data.sc)
  fitall.all.adj <- p.adjust(fitall.all, 'fdr')

  XX <- rownames(counts(sce))[1:13992]
  XXX <- rownames(counts(sce))[13993:17490]

  YY <- intersect(XX, names(fitall.all.adj))
  YYY <- intersect(XXX, names(fitall.all.adj))
  sum(fitall.all.adj[YY] < .05) / (sum(fitall.all.adj[YYY] < .05)+sum(fitall.all.adj[YY] < .05))
  sum(fitall.all.adj[YYY] < .05) / (sum(fitall.all.adj[YYY] < .05)+sum(fitall.all.adj[YY] < .05))

  YY <- intersect(XX, names(fitall.all.eq))
  YYY <- intersect(XXX, names(fitall.all.eq))
  sum(fitall.all.eq.adj[YY] < .05) / (sum(fitall.all.eq.adj[YYY] < .05)+sum(fitall.all.eq.adj[YY] < .05))
  sum(fitall.all.eq.adj[YYY] < .05) / (sum(fitall.all.eq.adj[YYY] < .05)+sum(fitall.all.eq.adj[YY] < .05))

  
  ## Calculate the AUC
YY <- intersect(XX, names(fitall.all.eq))
YYY <- intersect(XXX, names(fitall.all.eq))
ordered.eq <- fitall.all.eq.adj[c(YY,YYY)]
ordered.eq[ordered.eq < .05] <- 1
ordered.eq[ordered.eq != 1] <- 0
df <- data.frame(predictions = ordered.eq, labels=c(rep(0, length(YY)),rep(1, length(YYY))))
library(pROC)
pROC_obj.eq <- roc(df$labels,df$predictions,
                smoothed = TRUE, percent=TRUE)

YY <- intersect(XX, names(fitall.all.adj))
YYY <- intersect(XXX, names(fitall.all.adj))
ordered.og <- fitall.all.adj[c(YY,YYY)]
ordered.og[ordered.og < .05] <- 1
ordered.og[ordered.og != 1] <- 0
df <- data.frame(predictions = ordered.og, labels=c(rep(0, length(YY)),rep(1, length(YYY))))
library(pROC)
pROC_obj <- roc(df$labels,df$predictions,percent=TRUE,
                smoothed = TRUE)

auc.eq[i] <- pROC_obj.eq$auc[1]
auc.uneq[i] <- pROC_obj$auc[1]
}

save.image("RDATA/sim_dynPop.RData")

load("RDATA/sim_dynPop.RData")
mean(auc.eq) # 83.02555
mean(auc.uneq) # 81.74244


