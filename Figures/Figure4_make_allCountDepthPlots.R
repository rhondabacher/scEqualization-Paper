# For all my public datasets make plots:
setwd("EqualizationManuscript")

source("CODE/functionsForCountDepthPlots.R")

#############
colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), bias=2)(n = 10)

datasets <- list.files(path="RDATA/Conquer_Processed/", pattern=".Rdata")

TotalCells_Initial <- c()
TotalCells_Filter <- c()
AverageDepth <- c()
AverageCDR <- c()
ModeStat_SquaredError <- list()

for(i in 1:length(datasets)) {
	
	load(paste0("RDATA/Conquer_Processed/",datasets[i]))
	counts <- round(counts.sub)
	counts <- counts[,which(colSums(counts) > 10000)]
	
  geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
  counts <- counts[which(geneMeans >= 2),]
  
	TotalCells_Initial[i] <- ncol(counts.sub)
	TotalCells_Filter[i] <- ncol(counts)
	AverageDepth[i] <- mean(colSums(counts))
  AverageCDR[i] <- mean(colMeans(counts != 0))
	
  if(ncol(counts) > 10) {
    gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

    MedExpr <- log(geneMeans[names(gslopes)])
    splitby <- sort(MedExpr)
    grps <- length(splitby) / 10
    sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

    pdf(paste0("PLOTS/countDepth_",gsub(".Rdata", "", datasets[i]),".pdf"), height=10, width=5, useDingbats = F)
      par(mar=c(5,5,2,1), mfrow=c(2,1))

      makePlot(SLOPES = gslopes, TYPE = "Quantile", sreg)
      ModeStat_SquaredError[[i]] <- dotPlot(SLOPES = gslopes, sreg, NAME = "Quantile")

    dev.off()

    print(i)
  }
}
names(TotalCells_Initial) <- gsub(".Rdata", "", datasets)
names(TotalCells_Filter) <- gsub(".Rdata", "", datasets)
names(AverageDepth) <- gsub(".Rdata", "", datasets)
names(AverageCDR) <- gsub(".Rdata", "", datasets)
names(ModeStat_SquaredError) <- gsub(".Rdata", "", datasets)

save(TotalCells_Initial, 
TotalCells_Filter,
AverageDepth,
AverageCDR,
ModeStat_SquaredError, file="RDATA/countDepth_Stats_PublicData.Rdata")

##############################################################################

# Plot one example from each dataset:

ModeStat_SquaredError <- NULL
## "countDepth_GSE45719-Deng_earlyblast_Embryo2.pdf"
load(paste0("RDATA/Conquer_Processed/",datasets[4]))
counts <- round(counts.sub)
counts <- counts[,which(colSums(counts) > 10000)]

geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
counts <- counts[which(geneMeans >= 2),]

gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

MedExpr <- log(geneMeans[names(gslopes)])
splitby <- sort(MedExpr)
grps <- length(splitby) / 10
sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

pdf(paste0("PLOTS/countDepth_",gsub(".Rdata", "", datasets[4]),".pdf"), height=7, width=3, useDingbats = F)
  par(mar=c(5,5,2,1), mfrow=c(2,1))

  makePlot(SLOPES = gslopes, TYPE = "", sreg)
  ModeStat_SquaredErrorX <- dotPlot(SLOPES = gslopes, sreg, NAME = "")

dev.off()
  
# countDepth_GSE63818-GPL16791-Guo_M11W_Embryo2.pdf
load(paste0("RDATA/Conquer_Processed/",datasets[32]))
counts <- round(counts.sub)
counts <- counts[,which(colSums(counts) > 10000)]

geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
counts <- counts[which(geneMeans >= 2),]

gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

MedExpr <- log(geneMeans[names(gslopes)])
splitby <- sort(MedExpr)
grps <- length(splitby) / 10
sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

pdf(paste0("PLOTS/countDepth_",gsub(".Rdata", "", datasets[32]),".pdf"), height=7, width=3, useDingbats = F)
  par(mar=c(5,5,2,1), mfrow=c(2,1))

  makePlot(SLOPES = gslopes, TYPE = "", sreg)
  ModeStat_SquaredErrorX <- dotPlot(SLOPES = gslopes, sreg, NAME = "")

dev.off()

  
# "GSE48968-GPL13112-Shalek_UNSTIM-REP1.Rdata"
load(paste0("RDATA/Conquer_Processed/",datasets[25]))
counts <- round(counts.sub)
counts <- counts[,which(colSums(counts) > 10000)]

geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
counts <- counts[which(geneMeans >= 2),]

gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

MedExpr <- log(geneMeans[names(gslopes)])
splitby <- sort(MedExpr)
grps <- length(splitby) / 10
sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

pdf(paste0("PLOTS/countDepth_",gsub(".Rdata", "", datasets[25]),".pdf"), height=7, width=3, useDingbats = F)
  par(mar=c(5,5,2,1), mfrow=c(2,1))

  makePlot(SLOPES = gslopes, TYPE = "", sreg)
  ModeStat_SquaredErrorX <- dotPlot(SLOPES = gslopes, sreg, NAME = "")

dev.off()
##############################################################################
##############################################################################

load("RDATA/countDepth_Stats_PublicData.Rdata")
# i = 37
i = i + 1
load("RDATA/picelli_data.Rdata")
	counts <- round(picellidata.orig)
	TotalCells_Initial[i] <- ncol(counts)
	counts <- counts[,which(colSums(counts) > 10000)]
	TotalCells_Filter[i] <- ncol(counts)
	AverageDepth[i] <- mean(colSums(counts))
	AverageCDR[i] <- mean(colMeans(counts != 0))
	
  geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
  counts <- counts[which(geneMeans >= 2),]
	gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

	MedExpr <- log(geneMeans[names(gslopes)])
	splitby <- sort(MedExpr) 
	grps <- length(splitby) / 10
	sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

	colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), bias=2)(n = 10)

	pdf(paste0("PLOTS/countDepth_PicelliOrig.pdf"), height=7, width=3, useDingbats = F)
		par(mar=c(5,5,2,1), mfrow=c(2,1))
		makePlot(SLOPES = gslopes, TYPE = "", sreg)
		
		ModeStat_SquaredError[[i]] <- dotPlot(SLOPES = gslopes, sreg, NAME = "")
		
	dev.off()


i = i + 1	
load("RDATA/islam_data.Rdata")
	counts <- round(islamES)
	TotalCells_Initial[i] <- ncol(counts)
	counts <- counts[,which(colSums(counts) > 10000)]
	TotalCells_Filter[i] <- ncol(counts)
	AverageDepth[i] <- mean(colSums(counts))
	AverageCDR[i] <- mean(colMeans(counts != 0))
	
  geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
  counts <- counts[which(geneMeans >= 2),]
	gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

	MedExpr <- log(geneMeans[names(gslopes)])
	splitby <- sort(MedExpr) 
	grps <- length(splitby) / 10
	sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

	pdf(paste0("PLOTS/countDepth_IslamES.pdf"), height=10, width=5, useDingbats = F)
		par(mar=c(5,5,2,1), mfrow=c(2,1))
		makePlot(SLOPES = gslopes, TYPE = "Quantile", sreg)

		ModeStat_SquaredError[[i]] <- dotPlot(SLOPES = gslopes, sreg, NAME = "Quantile")
dev.off()	

i = i + 1	
	counts <- round(islamEF)
	TotalCells_Initial[i] <- ncol(counts)
	counts <- counts[,which(colSums(counts) > 10000)]
	TotalCells_Filter[i] <- ncol(counts)
	AverageDepth[i] <- mean(colSums(counts))
	AverageCDR[i] <- mean(colMeans(counts != 0))
	
  geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
  counts <- counts[which(geneMeans >= 2),]
	gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

	MedExpr <- log(geneMeans[names(gslopes)])
	splitby <- sort(MedExpr) 
	grps <- length(splitby) / 10
	sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

	pdf(paste0("PLOTS/countDepth_IslamEF.pdf"), height=7, width=3, useDingbats = F)
		par(mar=c(5,5,2,1), mfrow=c(2,1))
		makePlot(SLOPES = gslopes, TYPE = "", sreg)
		
		ModeStat_SquaredError[[i]] <- dotPlot(SLOPES = gslopes, sreg, NAME = "")
	dev.off()

	
i = i + 1	
load("RDATA/thomsonLab_bulkdata_USE.Rdata")
	counts <- round(bulkH1data)
	TotalCells_Initial[i] <- ncol(counts)
	counts <- counts[,which(colSums(counts) > 10000)]
	TotalCells_Filter[i] <- ncol(counts)
	AverageDepth[i] <- mean(colSums(counts))
  AverageCDR[i] <- mean(colMeans(counts != 0))
  
  geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
  counts <- counts[which(geneMeans >= 2),]
	gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

	MedExpr <- log(geneMeans[names(gslopes)])
	splitby <- sort(MedExpr) 
	grps <- length(splitby) / 10
	sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

	pdf(paste0("PLOTS/countDepth_BulkH1.pdf"), height=7, width=3, useDingbats = F)
		par(mar=c(5,5,2,1), mfrow=c(2,1))
		makePlot(SLOPES = gslopes, TYPE = "", sreg)
		
		ModeStat_SquaredError[[i]] <- dotPlot(SLOPES = gslopes, sreg, NAME = "")
dev.off()


load("RDATA/thomsonLab_scdata_fromGEO.Rdata")
for (j in 1:15) {
  print(i)
  i = i + 1	
	counts <- round(data.matrix(thomsondata[,batches$start[j]:batches$end[j]]))
  
  TotalCells_Initial[i] <- ncol(counts)
	counts <- counts[,which(colSums(counts) > 10000)]
  TotalCells_Filter[i] <- ncol(counts)
  AverageDepth[i] <- mean(colSums(counts))
  AverageCDR[i] <- mean(colMeans(counts != 0))

  geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
  counts <- counts[which(geneMeans >= 2),]
	gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

	MedExpr <- log(geneMeans[names(gslopes)]) 
	splitby <- sort(MedExpr) 
	grps <- length(splitby) / 10
	sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

  pdf(paste0("PLOTS/countDepth_",batches$batch[j],".pdf"), height=10, width=5, useDingbats = F)
    par(mar=c(5,5,2,1), mfrow=c(2,1))
    makePlot(SLOPES = gslopes, TYPE = "Quantile", sreg)
		ModeStat_SquaredError[[i]] <- dotPlot(SLOPES = gslopes, sreg, NAME = "Quantile")
  dev.off()
  print(j)
}

# Plot one specific example:
j=15
counts <- round(data.matrix(thomsondata[,batches$start[j]:batches$end[j]]))

geneMeans <- apply(counts, 1, function(x) median(x[x!=0]))
counts <- counts[which(geneMeans >= 2),]
gslopes <- getSlopes(Data = counts, ditherCounts=FALSE)

MedExpr <- log(geneMeans[names(gslopes)]) 
splitby <- sort(MedExpr) 
grps <- length(splitby) / 10
sreg <- split(splitby, ceiling(seq_along(splitby) / grps))

pdf(paste0("PLOTS/countDepth_",batches$batch[j],".pdf"), height=7, width=3, useDingbats = F)
  par(mar=c(5,5,2,1), mfrow=c(2,1))
  makePlot(SLOPES = gslopes, TYPE = "", sreg)
	ModeStat_SquaredErrorX <- dotPlot(SLOPES = gslopes, sreg, NAME = "")
dev.off()

################################################################################################################

#### Now plot the statistics together

names(ModeStat_SquaredError)[38:41] <- c("Picelli1",  "IslamES", "IslamEF", "BulkH1")
names(TotalCells_Initial)[38:41] <- c("Picelli1", "IslamES", "IslamEF", "BulkH1")
names(TotalCells_Filter)[38:41] <- c("Picelli1", "IslamES", "IslamEF", "BulkH1")
names(AverageDepth)[38:41] <- c("Picelli1",  "IslamES", "IslamEF", "BulkH1") 
names(AverageCDR)[38:41] <- c("Picelli1",  "IslamES", "IslamEF", "BulkH1")

names(TotalCells_Initial)[42:56] <- as.character(batches$batch)
names(TotalCells_Filter)[42:56] <- as.character(batches$batch)
names(AverageDepth)[42:56] <- as.character(batches$batch)
names(AverageCDR)[42:56] <- as.character(batches$batch)
names(ModeStat_SquaredError)[42:56] <- as.character(batches$batch)


save(TotalCells_Initial, 
TotalCells_Filter,
AverageDepth,
AverageCDR,
ModeStat_SquaredError, file="RDATA/countDepth_Stats_ALLPublicData.Rdata")

############################################################################################################


load("RDATA/countDepth_Stats_ALLPublicData.Rdata")

library(ggsci)
library("scales")
show_col(pal_lancet("lanonc", alpha = .5)(7))


ModeStatistic <- unlist(lapply(ModeStat_SquaredError, function(x) median(abs(x[1:10] - 1))))
cbind(1:length(ModeStatistic), ModeStatistic)


min(ModeStatistic)
max(ModeStatistic)
sd(ModeStatistic)

sd(ModeStatistic[1:11]) 
sd(ModeStatistic[27:37])
sd(ModeStatistic[12:26])
sd(ModeStatistic[42:56])
sd(ModeStatistic[39:40])

median(abs(ModeStatistic - median(ModeStatistic)))

median(abs(ModeStatistic[1:11] - median(ModeStatistic[1:11])))
median(abs(ModeStatistic[27:37] - median(ModeStatistic[27:37])))
median(abs(ModeStatistic[12:26] - median(ModeStatistic[12:26])))
median(abs(ModeStatistic[42:56] - median(ModeStatistic[42:56])))


ModeStatistic <- ModeStatistic[c(41, 38, 1:11, 27:37,  12:26, 39:40, 42:56)]
useCols.Main <- pal_lancet("lanonc", alpha = 1)(7)
useCols.Dots <- pal_lancet("lanonc", alpha = .5)(7)
useCols.Dots <- c(rep(useCols.Dots[1], 1), 
                  rep(useCols.Dots[2], 1), 
                  rep(useCols.Dots[3], 11), 
                  rep(useCols.Dots[4], 11), 
                  rep(useCols.Dots[5], 15), 
                  rep(useCols.Dots[6], 2),
                  rep(useCols.Dots[7], 15))


pdf("PLOTS/MAD_AllDatasets_Statistic.pdf", height=2, width=6.6)
par(mar=c(2,2,1,1))
plot(1:length(ModeStatistic), (ModeStatistic), pch=16, col=useCols.Dots, xaxt='n', yaxt='n',
     xlab="", ylab = "", cex=1.5, cex.axis=1, ylim=c(0, 1.2))
axis(2, at=c(0, .3, .6, .9, 1.2), cex.axis=1)
lines(c(.5, 1.5), rep(mean(ModeStatistic[1]), 2), pch=16, col=useCols.Main[1], lwd=2)
lines(c(1.5:2.5), rep(mean(ModeStatistic[2]), 2), pch=16, col=useCols.Main[2], lwd=2)
lines(3:13, rep(mean(ModeStatistic[3:13]), 11), pch=16, col=useCols.Main[3], lwd=2)
lines(3:13, rep(mean(ModeStatistic[3:13])+sd(ModeStatistic[3:13]), 11), pch=16, col=useCols.Main[3], lwd=1, lty=2)
lines(3:13, rep(mean(ModeStatistic[3:13])-sd(ModeStatistic[3:13]), 11), pch=16, col=useCols.Main[3], lwd=1, lty=2)
lines(14:24, rep(mean(ModeStatistic[14:24]), 11), pch=16, col=useCols.Main[4], lwd=2)
lines(14:24, rep(mean(ModeStatistic[14:24]), 11), pch=16, col=useCols.Main[4], lwd=1, lty=2)
lines(14:24, rep(mean(ModeStatistic[14:24])+sd(ModeStatistic[14:24]), 11), pch=16, col=useCols.Main[4], lwd=1, lty=2)
lines(14:24, rep(mean(ModeStatistic[14:24])-sd(ModeStatistic[14:24]), 11), pch=16, col=useCols.Main[4], lwd=1, lty=2)

lines(25:39, rep(mean(ModeStatistic[25:39]), 15), pch=16, col=useCols.Main[5], lwd=2)
lines(25:39, rep(mean(ModeStatistic[25:39])+sd(ModeStatistic[25:39]), 15), pch=16, col=useCols.Main[5], lwd=1, lty=2)
lines(25:39, rep(mean(ModeStatistic[25:39])-sd(ModeStatistic[25:39]), 15), pch=16, col=useCols.Main[5], lwd=1, lty=2)

lines(40:41, rep(mean(ModeStatistic[40:41]), 2), pch=16, col=useCols.Main[6], lwd=2)
lines(40:41, rep(mean(ModeStatistic[40:41])+sd(ModeStatistic[40:41]), 2), pch=16, col=useCols.Main[6], lwd=1, lty=2)
lines(40:41, rep(mean(ModeStatistic[40:41])-sd(ModeStatistic[40:41]), 2), pch=16, col=useCols.Main[6], lwd=1, lty=2)

lines(42:56, rep(mean(ModeStatistic[42:56]), 15), pch=16, col=useCols.Main[7], lwd=2)
lines(42:56, rep(mean(ModeStatistic[42:56])+sd(ModeStatistic[42:56]), 15), pch=16, col=useCols.Main[7], lwd=1, lty=2)
lines(42:56, rep(mean(ModeStatistic[42:56])-sd(ModeStatistic[42:56]), 15), pch=16, col=useCols.Main[7], lwd=1, lty=2)

dev.off()




mean(ModeStatistic[1])
mean(ModeStatistic[2]) 
mean(ModeStatistic[3:13])
mean(ModeStatistic[14:24])
mean(ModeStatistic[25:39])
mean(ModeStatistic[40:41])
mean(ModeStatistic[42:56])

summary(TotalCells_Initial[1:11])
summary(TotalCells_Initial[12:26])
summary(TotalCells_Initial[27:37])
summary(TotalCells_Initial[39:40])
summary(TotalCells_Initial[42:56])

mean(AverageDepth[1:11]) / 1e6
mean(AverageDepth[12:26]) / 1e6
mean(AverageDepth[27:37]) / 1e6
(AverageDepth[c(38,41)])  / 1e6
mean(AverageDepth[39:40])  / 1e6
mean(AverageDepth[42:56])  / 1e6


mean(AverageCDR[1:11]) 
mean(AverageCDR[12:26])
mean(AverageCDR[27:37])
(AverageCDR[c(38,41)]) 
mean(AverageCDR[39:40]) 
mean(AverageCDR[42:56]) 

