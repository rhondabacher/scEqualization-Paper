setwd("~/OneDrive - University of Florida/EqualizationManuscript/")
load("RDATA/compare_SIM-properties_unEQvsunEQ_v2.RData")

meanDEV.0 <- meanDEV
meanDEPTH.0 <- meanDEPTH
meanMEAN.0 <- meanMEAN
meanGLD.0 <- meanGLD
meanCLD.0 <- meanCLD

load("RDATA/compare_SIM-properties_unEQvsEQ_v2.RData")


pdf("PLOTS/comparisonSIM_percents.pdf", height=5, width=10)
par(mfrow=c(1,2), mar=c(3,5,2,1))
boxplot(meanDEV.0*100, meanDEV*100, ylab="Percent", cex.axis=1.3,cex.lab=1.5,
        main="Reduction in Gene-Specific Expression Variability",
          names=c("unEQ vs. unEQ", "unEQ vs. EQ"))
boxplot(meanDEPTH.0*100, meanDEPTH*100, ylab="Percent", cex.axis=1.3,cex.lab=1.5,
        main="Reduction in Sequencing Depth Variability",
        names=c("unEQ vs. unEQ", "unEQ vs. EQ"))
dev.off()
## Make supplemental figure out of this. 



boxplot(meanMEAN.0, meanMEAN, ylab="Difference in Sequencing Depths", 
        main="Difference in Sequencing Depths",
        names=c("unEQ vs. unEQ", "unEQ vs. EQ"))

boxplot(meanCLD.0, meanCLD, ylab="Difference in Sequencing Depths", 
        main="Difference in Sequencing Depths",
        names=c("unEQ vs. unEQ", "unEQ vs. EQ"))

boxplot(meanLD.0, meanGLD, ylab="Difference in Sequencing Depths", 
        main="Difference in Sequencing Depths",
        names=c("unEQ vs. unEQ", "unEQ vs. EQ"))

