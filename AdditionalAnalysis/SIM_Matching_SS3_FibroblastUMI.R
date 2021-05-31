
load("RDATA/Sim_SS3_Fibroblast-reads.RData")


ss3.umi <- read.table("DATA/Smartseq3/Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt", header=TRUE)

gslopes0 <- fasterSlope(Data = data.matrix(ss3.umi), ditherCounts=TRUE)

countsSim <- simdata@assays@data$umi_counts
compareData <- ss3.umi

gslopes1 <- fasterSlope(Data = data.matrix(countsSim), ditherCounts = TRUE)

save.image("RDATA/BestSim_SS3_Fibroblast-umi.RData")


makePlot <- function(SLOPES, TYPE, sreg) {
        colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), bias=2)(n = length(sreg))
        
        Height <- rep(NA, length(sreg))
        for (i in 1:length(sreg)) {
                useg <- names(sreg[[i]])
                rqdens <- density(na.omit(SLOPES[useg]))
                peak <- which.max(rqdens$y)
                Height[i] <- rqdens$y[peak]
        }
        YMAX <- max(Height)
        YMAX <- ceiling(YMAX/0.5)*0.5 + .1
        YMAX <- max(YMAX, 1)
        plot(density(na.omit(SLOPES), from=min(SLOPES, na.rm=T), to=max(SLOPES, na.rm=T), adjust=1), 
             xlab="Slope", ylab="Density",  main=TYPE, 
             cex.lab=2, cex.main=1.5, cex.axis=2,xlim=c(-2,3), 
             ylim=c(0,10), lwd=1, col="white", xaxt='n', yaxt='n', bty='n')
        
        
        for (i in length(sreg):1) {
                topgenes<-names(sreg[[i]])
                lines(density(na.omit(SLOPES[topgenes]), from=-3, to=3, adjust=1), lwd=1, col=colors[i])
        }
        abline(v=0, lwd=1, col="black")
        axis(side = 1, at=seq(-2,3, by=1), lwd.ticks=1, lwd = 1, cex.axis=1.6,font=1)
        axis(side = 2, at=seq(0,YMAX,by=.5), lwd = 1, lwd.ticks=1, cex.axis=1.6,font=1)
        abline(v=1, lwd=1, col="black")
        
}


pdf("PLOTS/sim_CountDepth_SS3-Fibroblast-umi.pdf", height=8, width=8, useDingbats=F)
ModeStat_SquaredError<-list()
MedExpr <- (apply(ss3.data, 1, function(x) mean(x)))
splitby <- sort(MedExpr[MedExpr > 1])
grps <- length(splitby) / 10
sreg_orig1 <- split(splitby, ceiling(seq_along(splitby) / grps))

par(mfrow=c(2,2), mar=c(5,5,2,1))

i = 1
makePlot(gslopes0, " ", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes0, sreg_orig1, "Original")

i = i + 1
makePlot(gslopes1, " ", sreg_orig1)
ModeStat_SquaredError[[i]] <- dotPlot(gslopes1, sreg_orig1, "SIM")

ModeStatistic_EC <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))

print(paste0(ModeStatistic_EC))
## "0.135029354207436" "0.26027397260274" 
dev.off()
# 

library(yarrr)


pdf("PLOTS/sim_SeqDepth_SS3-Fibroblast-umi.pdf", height=4, width=8, useDingbats=F)
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




pdf("PLOTS/sim_CellDetectRate_SS3-Fibroblast-umi.pdf", height=4, width=8, useDingbats=F)
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


pdf("PLOTS/sim_GeneDetectRate_SS3-Fibroblast-umi.pdf", height=4, width=8, useDingbats=F)
set.seed(111)
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



pdf("PLOTS/sim_GeneMean_SS3-Fibroblast-umi.pdf", height=4, width=8, useDingbats=F)

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



pdf("PLOTS/sim_GeneSD_SS3-Fibroblast-umi.pdf", height=4, width=8, useDingbats=F)
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







