setwd("EqualizationManuscript")
load(file="RDATA/ALL_exprData_trimmed_genes.RData")

library(ggcorrplot)

XX <- cbind(sub147[genez,ec.cells], sub637[genez,ec.cells])
colnames(XX) <- c(paste0("EQ",1:45), paste0("EQ-Vary",1:45))
corr <- round(cor(XX), 1)
corr <- corr[1:45,46:90]
G1 <- ggcorrplot(corr,outline.col = "white", title="EC: unEQ - EQ",
           ggtheme = ggplot2::theme_gray,show.legend=F,
           tl.cex=8,tl.srt = 90,
           colors = c("#6D9EC1", "white", "#E46726"))


XX <- cbind(sub147[genez,ec.cells], sub644[genez,ec.cells])
colnames(XX) <- c(paste0("EQ",1:45), paste0("EQ-Vary",1:45))
corr <- round(cor(XX), 2)
corr <- corr[1:45,46:90]
G2 <- ggcorrplot(corr,outline.col = "white", title="EC: unEQ - EQ-Vary",
           ggtheme = ggplot2::theme_gray,show.legend=F,
           tl.cex=8,tl.srt = 90,
           colors = c("#6D9EC1", "white", "#E46726"))
           
XX <- cbind(sub147[genez,ec.cells], sub671[genez,ec.cells])
colnames(XX) <- c(paste0("EQ",1:45), paste0("EQ-Vary",1:45))
corr <- round(cor(XX), 1)
corr <- corr[1:45,46:90]
G3 <- ggcorrplot(corr,outline.col = "white", title="EC: unEQ - EQ-75%",
           ggtheme = ggplot2::theme_gray,show.legend=F,
           tl.cex=8,tl.srt = 90,
           colors = c("#6D9EC1", "white", "#E46726"))








XX <- cbind(sub150[genez,tb.cells], sub637[genez,tb.cells])
colnames(XX) <- c(paste0("EQ",1:51), paste0("EQ-Vary",1:51))
corr <- round(cor(XX), 1)
corr <- corr[1:51,52:102]
G4 <- ggcorrplot(corr,outline.col = "white", title="TB: unEQ - EQ",
           ggtheme = ggplot2::theme_gray,show.legend=F,
           tl.cex=8,tl.srt = 90,
           colors = c("#6D9EC1", "white", "#E46726"))


XX <- cbind(sub150[genez,tb.cells], sub644[genez,tb.cells])
colnames(XX) <- c(paste0("EQ",1:51), paste0("EQ-Vary",1:51))
corr <- round(cor(XX), 2)
corr <- corr[1:51,52:102]
G5 <- ggcorrplot(corr,outline.col = "white", title="TB: unEQ - EQ-Vary",
           ggtheme = ggplot2::theme_gray,show.legend=T,
           tl.cex=8,tl.srt = 90,
           colors = c("#6D9EC1", "white", "#E46726"))

XX <- cbind(sub150[genez,tb.cells], sub671[genez,tb.cells])
colnames(XX) <- c(paste0("EQ",1:51), paste0("EQ-Vary",1:51))
corr <- round(cor(XX), 1)
corr <- corr[1:51,52:102]
G6 <- ggcorrplot(corr,outline.col = "white", title="TB: unEQ - EQ-Vary",
           ggtheme = ggplot2::theme_gray,show.legend=F,
           tl.cex=8,tl.srt = 90,
           colors = c("#6D9EC1", "white", "#E46726"))

library(gridExtra)
pdf("PLOTS/cellCorrelations_1.pdf", height=6, width=8, useDingbats = F)
print(G1)
dev.off()
pdf("PLOTS/cellCorrelations_2.pdf", height=6, width=8, useDingbats = F)
print(G2)
dev.off()
pdf("PLOTS/cellCorrelations_3.pdf", height=6, width=8, useDingbats = F)
print(G3)
dev.off()
pdf("PLOTS/cellCorrelations_4.pdf", height=6, width=8, useDingbats = F)
print(G4)
dev.off()
pdf("PLOTS/cellCorrelations_5.pdf", height=6, width=8, useDingbats = F)
print(G5)
dev.off()
pdf("PLOTS/cellCorrelations_6.pdf", height=6, width=8, useDingbats = F)
print(G6)
dev.off()
