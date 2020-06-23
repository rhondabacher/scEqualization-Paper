# Process other publica datasets for easy access (save as .RData)

setwd("EqualizationManuscript")

# Thomson lab batches
thomsondata <- read.csv("DATA/GSE75748_sc_cell_type_ec.csv", header=T, row.names=1)
colnames(thomsondata)
thomsondata <- thomsondata[,sort(colnames(thomsondata))]

# thomsondata.more <- read.csv("DATA/GSE75748_sc_time_course_ec.csv", header=T, row.names=1)
# colnames(thomsondata.more)

cellNames <- table(unlist(lapply(strsplit(colnames(thomsondata), ".", fixed=TRUE), function(x) x[1])))
head(cellNames)
batches <- NULL
batches$batch <- names(cellNames)
batches$cells <- cellNames
batches$start <- c(1, 65, 139, 203, 244, 304, 382, 456, 543, 618, 667, 746, 777, 863, 950)
batches$end <- c(64, 138, 202, 243, 303, 381, 455, 542, 617, 666, 745, 776, 862, 949, 1018)


save(batches, thomsondata, file="RDATA/thomsonLab_scdata_fromGEO.Rdata")


# load("DATA/scdata.Rdata") #scdata
# load("DATA/cells_perbatch.Rdata") #batches
# 
# save.image(file="RDATA/thomsonLab_scdata.Rdata")


# Picelli data

txt_files = list.files(path="DATA/GSE49321", pattern="*.txt")
data_list = lapply(txt_files, function(x) read.table(paste0("DATA/GSE49321/", x), header=F)[,c(4)])
picellidata <- data.frame(do.call(cbind, data_list))
picellidata.orig <- picellidata
rownames(picellidata.orig) <- paste0("Gene_", 1:nrow(picellidata))
picellidata.orig <- data.matrix(picellidata.orig)


picellidata <- cbind(Gene = read.table(paste0("DATA/GSE49321/", txt_files[1]), header=F, stringsAsFactors=F)[,c(1)], 
						picellidata, stringsAsFactors=F)
dupg <- unique(picellidata[,1][which(duplicated(picellidata[,1]))])
for(i in 1:length(dupg)) {
	toRM <- which(picellidata$Gene == dupg[i])
	newG <- c(Gene = i, colSums(picellidata[toRM, 2:ncol(picellidata)]))
	picellidata <- picellidata[-toRM,]
	picellidata <- rbind(picellidata, newG)
}
dim(picellidata)
picellidata[1:5,1:5]
which(picellidata[,1] == "1")

picellidata[which(picellidata[,1] == "1"):nrow(picellidata),1] <- dupg

rownames(picellidata) <- picellidata[,1]
picellidata = picellidata[,-1]

picellidata <- data.matrix(picellidata)

save(picellidata, picellidata.orig, file="RDATA/picelli_data.Rdata")

##################################################################################################################
# Islam data

alldata<- read.table("DATA/GSE29087_L139_expression_tab_fix.txt", 
                 stringsAsFactors=F, header=T, fill=TRUE, row.names=1)
alldata[1:4,1:10]
alldata <- alldata[,-c(93:96)] #remove negative controls, 93,94,95,96
dim(alldata)

islamES <- data.matrix(alldata[,1:48])
islamEF <- data.matrix(alldata[,49:92])

save(islamES, islamEF, file="RDATA/islam_data.Rdata")

######################################################################################################################
# Thomson lab bulk data

library(readxl)

inData <- read_excel("DATA/GSE85917_Bacher.RSEM.xlsx", sheet=5)
head(inData); dim(inData)
bulkH1data <- data.matrix(inData[,-1])
rownames(bulkH1data) <- inData$Gene.ID
save(bulkH1data, file="RDATA/thomsonLab_bulkdata_USE.Rdata")

######################################################################################################################

