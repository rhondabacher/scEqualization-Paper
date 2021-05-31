load("RDATA/forRevision/Sim_SS3_HCA-reads.RData")

ss3.umi <- read.table("DATA/Smartseq3/HCA.UMIcounts.PBMC.txt", header=TRUE)

ss3.umi <- ss3.umi[,colnames(ss3.data)]
countsSim <- simdata@assays@data$umi_counts
compareData <- ss3.umi

library(yarrr)

pdf("simMatching_SeqDepth_SS3-HCA-umi.pdf", height=4, width=8, useDingbats=F)
X <- data.frame( Depth = colSums(countsSim) , Species = "Simulated")
Y <- data.frame( Depth = colSums(compareData) , Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata, 
           xlab = "", 
           ylab = "Sequencing Depth", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Sequencing Depth (millions)", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()




pdf("simMatching_CellDetectRate_SS3-HCA-umi.pdf", height=4, width=8, useDingbats=F)
X <- data.frame( Depth = colSums(countsSim!=0) / nrow(countsSim), Species = "Simulated")
Y <- data.frame( Depth = colSums(compareData!=0) / nrow(countsSim), Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species, 
           data = longdata, ylim=c(0,1),
           xlab = "", 
           ylab = "Detection Rate", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Detection Rate", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()


pdf("simMatching_GeneDetectRate_SS3-HCA-umi.pdf", height=4, width=8, useDingbats=F)
set.seed(777)
Genes = rownames(countsSim)
XX <- sample(Genes, 200)
X1 <- apply(countsSim[XX,], 1, function(x) sum(x!=0, na.rm=T)) / dim(countsSim)[2]
X2 <- apply(compareData[XX,], 1, function(x) sum(x!=0, na.rm=T)) / dim(compareData)[2]
X <- data.frame( Depth =X1, Species = "Simulated")
Y <- data.frame( Depth = X2, Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "Detection Rate", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="Detection Rate", 
     cex.axis=1.3, cex.lab=1.5, cex.main=1.5, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "EC Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("simMatching_GeneMean_SS3-HCA-umi.pdf", height=4, width=8, useDingbats=F)

X1 <- log(apply(countsSim[XX,], 1, function(x) mean(x, na.rm=T))+1)
X2 <- log(apply(compareData[XX,], 1, function(x) mean(x, na.rm=T))+1)
useg <- names(which(X1<Inf & X2 < Inf))
X <- data.frame( Depth =X1[useg], Species = "Simulated")
Y <- data.frame( Depth = X2[useg], Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "log (mean+1)", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="log (mean+1)", 
     cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()



pdf("simMatching_GeneSD_SS3-HCA-umi.pdf", height=4, width=8, useDingbats=F)
X1 <- log(apply(countsSim[XX,], 1, function(x) sd(x, na.rm=T))+1)
X2 <- log(apply(compareData[XX,], 1, function(x) sd(x, na.rm=T))+1)
useg <- names(which(X1<Inf & X2 < Inf))
X <- data.frame( Depth =X1[useg], Species = "Simulated")
Y <- data.frame( Depth = X2[useg], Species = "Original")
longdata <- rbind(Y, X)
par(mfrow=c(1,2))
par(mar=c(4,5,2,1), mgp = c(3, .5, 0))
pirateplot(formula = Depth ~ Species,
           data = longdata,
           xlab = "", 
           ylab = "log (sd+1)", pal=c("cornflowerblue", "brown1"),
           main = "", point.cex=1.1, bar.lwd=1, cex.lab=1.5, cex.axis=1.3,cex.names=1.5)
par(mar=c(5.5,4,2,1), mgp = c(2.5, 1, 0))
plot(ecdf(X$Depth), col="brown1", main="", xlab="log (sd+1)", 
     cex.axis=2, cex.lab=2, cex.main=2, xlim=c(0, max(X$Depth, Y$Depth)))
plot(ecdf(Y$Depth), add=T, col="cornflowerblue")
legend('bottomright', c("Simulated", "H1 Data"), col=c("brown1","cornflowerblue"), lwd=3)
dev.off()







