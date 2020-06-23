# Process data from Shalek, Guo, and Deng from conquer database:

setwd("EqualizationManuscript")

# Deng data, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
library("MultiAssayExperiment")
XX <- readRDS(paste0("RDATA/Conquer_Downloaded/GSE45719.rds"))
Y = experiments(XX)[[1]] #this gives us the GENE data
counts <- assay(Y, "count") # get the count matrix

#Break into 'conditions/groups'
counts.sub <- counts[,paste0("GSM11",seq(12490,12503, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_16cell_Embryo1.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12515,12527, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_16cell_Embryo5.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12528,12539, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_16cell_Embryo6.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12611,12625, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_earlyblast_Embryo2.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12626,12640, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_earlyblast_Embryo3.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12641,12653, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_earlyblast_Embryo4.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12664,12682, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_lateblast_Embryo1.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12683,12693, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_lateblast_Embryo2.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12706,12727, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_midblast_Embryo1.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12728,12747, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_midblast_Embryo2.Rdata")

counts.sub <- counts[,paste0("GSM11",seq(12748,12765, 1))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE45719-Deng_midblast_Embryo3.Rdata")

###########################

# Shalek data
library("MultiAssayExperiment")
XX <- readRDS(paste0("RDATA/Conquer_Downloaded/GSE48968-GPL13112.rds"))
Y = experiments(XX)[[1]]
counts <- assay(Y, "count")

# Break into groups

counts.sub <- counts[,paste0("GSM1189",seq(322,417, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_Shalek_UNSTIM-REP1.Rdata")

counts.sub <- counts[,paste0("GSM1189",seq(418,513, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_Shalek_LPS1h.Rdata")

counts.sub <- counts[,paste0("GSM1189",seq(514,609, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_Shalek_LPS2h.Rdata")

counts.sub <- counts[,paste0("GSM1189",seq(610,704, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_Shalek_LPS4h.Rdata")

counts.sub <- counts[,paste0("GSM1189",seq(705,800, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_Shalek_LPS6h.Rdata")

counts.sub <- counts[,paste0("GSM1189",seq(801,896, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_Shalek_UNSTIM-REP2.Rdata")

counts.sub <- counts[,paste0("GSM1189",seq(897,992, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_LPS4h-REP2.Rdata")

counts.sub <- counts[,c(paste0("GSM1189",seq(993,999, 1)), paste0("GSM119000",seq(1,9, 1)), paste0("GSM11900",seq(10,88, 1)))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_PAM1h.Rdata")

counts.sub <- counts[,c(paste0("GSM11900",seq(89,99, 1)), paste0("GSM1190",seq(100,184, 1)))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_PAM2h.Rdata")

counts.sub <- counts[,paste0("GSM1190",seq(185,269, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_PAM4h.Rdata")

counts.sub <- counts[,paste0("GSM1190",seq(270,333, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_PAM6h.Rdata")

counts.sub <- counts[,paste0("GSM1190",seq(334,429, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_PIC1h.Rdata")

counts.sub <- counts[,paste0("GSM1190",seq(430,513, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_PIC2h.Rdata")

counts.sub <- counts[,paste0("GSM1190",seq(514,603, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_PIC4h.Rdata")

counts.sub <- counts[,paste0("GSM1190",seq(604,699, 1))]
save(counts.sub, file="RDATA/Conquer_Processed/GSE48968-GPL13112-Shalek_PIC6h.Rdata")




#Guo dataset:

library("MultiAssayExperiment")
XX <- readRDS(paste0("RDATA/Conquer_Downloaded/GSE63818-GPL16791.rds"))
Y = experiments(XX)[[1]]
counts <- assay(Y, "count")

counts.sub <- counts[,c(paste0("GSM1677",seq(479,494, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_M7W_Embryo1.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(505,517, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_M7W_Embryo3.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(518,537, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_M10W_Embryo1.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(538,549, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_M11W_Embryo1.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(550,564, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_M11W_Embryo2.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(565,595, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_M19W_Embryo1.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(596,621, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_M19W_Embryo2.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(634,651, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_F8W_Embryo1.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(652,664, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_F10W_Embryo2.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(675,705, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_F17W_Embryo1.Rdata")

counts.sub <- counts[,c(paste0("GSM1677",seq(706,718, 1)))]; dim(counts.sub)
save(counts.sub, file="RDATA/Conquer_Processed/GSE63818-GPL16791-Guo_MSoma7W_Embryo1.Rdata")

