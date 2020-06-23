## All data is on hipergator. Run this script to save locally.
setwd("/Volumes/HPC/EQUALIZATION_PAPER/")

## Formatting EQ
mymat <- read.table("RSEMDATA/EQ/exprmat_EQ_trimmed_codingRef.txt", header=T, stringsAsFactors=F)
colnames(mymat) <- lapply(strsplit(colnames(mymat), "_"), function(x) paste0(x[2], "_", x[3]))
keymat = read.csv("FASTQ/Ln_1_d030c6e4175f81f9_EQ.csv", header=T, skip=1, stringsAsFactors=F)
keymat <- keymat[match(colnames(mymat), keymat$Sample_ID),]
colnames(mymat) <- keymat$Description
EQ <- mymat
save(EQ, file="RDATA/trimmedExpCountMat_genes_fromHPC_EQ.RDATA")

## Formatting EQ-Vary
mymat <- read.table("RSEMDATA/EQ-Vary/exprmat_EQ-Vary_trimmed_codingRef.txt", header=T, stringsAsFactors=F)
colnames(mymat) <- lapply(strsplit(colnames(mymat), "_"), function(x) paste0(x[2], "_", x[3]))
dim(mymat)
keymat = read.csv("FASTQ/Ln_2_15796e6ed1e36030_EQ-Vary.csv", header=T, skip=1, stringsAsFactors=F)
keymat <- keymat[match(colnames(mymat), keymat$Sample_ID),]
colnames(mymat) <- keymat$Description
EQ-Vary <- mymat
save(EQ-Vary, file="RDATA/trimmedExpCountMat_genes_fromHPC_EQ-Vary.RDATA")


## Formatting EQ-75pcnt
mymat <- read.table("RSEMDATA/EQ-75pcnt/exprmat_EQ-75pcnt_trimmed_codingRef.txt", header=T, stringsAsFactors=F)
colnames(mymat) <- lapply(strsplit(colnames(mymat), "_"), function(x) paste0(x[2], "_", x[3]))
dim(mymat)
keymat = read.csv("FASTQ/Ln_1_325587f4c3df848b_EQ-75pcnt.csv", header=T, skip=1, stringsAsFactors=F)
keymat <- keymat[match(colnames(mymat), keymat$Sample_ID),]
colnames(mymat) <- keymat$Description
EQ_75pcnt <- mymat
save(EQ_75pcnt, file="RDATA/trimmedExpCountMat_genes_fromHPC_EQ-75pcnt.RDATA")



## Formatting unEQ_TB
mymat <- read.table("RSEMDATA/unEQ_TB/exprmat_unEQ_TB_trimmed_codingRef.txt", header=T, stringsAsFactors=F)
colnames(mymat) <- lapply(strsplit(colnames(mymat), "_"), function(x) paste0(x[2], "_", x[3]))
colnames(mymat) <- gsub(".", "-", colnames(mymat), fixed=T)
dim(mymat)

keymat1 = read.csv("FASTQ/Ln_2_f4943063e178af60_unEQ_TB.csv", header=T, skip=0, stringsAsFactors=F)[-1,]
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_3_22ac9d93ee033462_unEQ_TB.csv", header=T, skip=0, stringsAsFactors=F)[-1,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_4_9a9fe938674362b7_unEQ_TB.csv", header=T, skip=0, stringsAsFactors=F)[-1,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_5_8f2215ea0bb500ba_unEQ_TB.csv", header=T, skip=0, stringsAsFactors=F)[-1,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_6_a59e5b531cc73c81_unEQ_TB.csv", header=T, skip=0, stringsAsFactors=F)[-1,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_7_0d1cbc4bb2a61dfc_unEQ_TB.csv", header=T, skip=0, stringsAsFactors=F)[-1,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_1_23d115d0864e9ea2_unEQ_TB.csv", header=T, skip=0, stringsAsFactors=F))
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_1_2fd5c3ae1330c0bf_unEQ_TB.csv", header=T, skip=0, stringsAsFactors=F)[-1,])

keymat1$Sample_ID <- paste0(sapply(strsplit(keymat1$SampleID.FlowcellSampleID, "_"), function(x) x[2]), "_", keymat1$Index)
keymat2 <- keymat1[match(colnames(mymat), keymat1$Sample_ID),]

keymat2$Description <- sapply(strsplit(keymat2$Description, "_"), function(x) {
  if(nchar(x[2]) < 3 ) {Y=paste0(x[1],"_", "0", x[2])}
    else {Y = paste0(x[1],"_",x[2])}
    return(Y)
})
colnames(mymat) <- keymat2$Description
colnames(mymat) <- gsub("TB", "TBb2", colnames(mymat), fixed=T)
unEQ_TB <- mymat[,which(colnames(mymat) %in% colnames(EQ))]
save(unEQ_TB, file="RDATA/trimmedExpCountMat_genes_fromHPC_unEQ_TB.RDATA")


## Formatting unEQ_EC
mymat <- read.table("RSEMDATA/unEQ_EC/exprmat_unEQ_EC_trimmed_codingRef.txt", header=T, stringsAsFactors=F)
colnames(mymat) <- lapply(strsplit(colnames(mymat), "_"), function(x) paste0(x[2], "_", x[3]))
colnames(mymat) <- gsub(".", "-", colnames(mymat), fixed=T)
dim(mymat)

keymat1 = read.csv("FASTQ/Ln_8_15d88367f863f22d_unEQ_EC.csv", header=T, skip=0, stringsAsFactors=F)[,]
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_1_f0585f8ca4d1b042_unEQ_EC.csv", header=T, skip=0, stringsAsFactors=F)[,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_2_62c0f3e9ccc44cb8_unEQ_EC.csv", header=T, skip=0, stringsAsFactors=F)[,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_3_666ef10660879702_unEQ_EC.csv", header=T, skip=0, stringsAsFactors=F)[,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_4_8865fae2fb5d4e04_unEQ_EC.csv", header=T, skip=0, stringsAsFactors=F)[,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_5_da3a2b8f8de5d143_unEQ_EC.csv", header=T, skip=0, stringsAsFactors=F)[,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_6_44273364a93c7547_unEQ_EC.csv", header=T, skip=0, stringsAsFactors=F)[,])
keymat1 = rbind(keymat1, read.csv("FASTQ/Ln_7_1c70ef08150caa7e_unEQ_EC.csv", header=T, skip=0, stringsAsFactors=F)[,])

keymat1$Sample_ID <- paste0(sapply(strsplit(keymat1$SampleID.FlowcellSampleID, "_"), function(x) x[2]), "_", keymat1$Index)
keymat2 <- keymat1[match(colnames(mymat), keymat1$Sample_ID),]

keymat2$Description <- sapply(strsplit(keymat2$Description, "_"), function(x) {
  if(nchar(x[2]) < 3 ) {Y=paste0(x[1],"_", "0", x[2])}
    else {Y = paste0(x[1],"_",x[2])}
    return(Y)
})
colnames(mymat) <- keymat2$Description
colnames(mymat) <- gsub("EC", "ECb2", colnames(mymat), fixed=T)

unEQ_EC <- mymat[,which(colnames(mymat) %in% colnames(EQ))]
save(unEQ_EC, file="RDATA/trimmedExpCountMat_genes_fromHPC_unEQ_EC.RDATA")




## Formatting EQ_R2
mymat <- read.table("RSEMDATA/EQ_R2/exprmat_EQ_R2_trimmed_codingRef.txt", header=T, stringsAsFactors=F)
colnames(mymat) <- lapply(strsplit(colnames(mymat), "_"), function(x) paste0(x[2], "_", x[3]))
keymat = read.csv("FASTQ/Ln_1_d030c6e4175f81f9_EQ.csv", header=T, skip=1, stringsAsFactors=F)
keymat <- keymat[match(colnames(mymat), keymat$Sample_ID),]
colnames(mymat) <- keymat$Description
EQ_R2 <- mymat
save(EQ_R2, file="RDATA/trimmedExpCountMat_genes_fromHPC_EQ_R2.RDATA")

## Formatting EQ-Vary
mymat <- read.table("RSEMDATA/EQ-Vary_R2/exprmat_EQ-Vary_R2_trimmed_codingRef.txt", header=T, stringsAsFactors=F)
colnames(mymat) <- lapply(strsplit(colnames(mymat), "_"), function(x) paste0(x[2], "_", x[3]))
dim(mymat)
keymat = read.csv("FASTQ/Ln_2_15796e6ed1e36030_EQ-Vary.csv", header=T, skip=1, stringsAsFactors=F)
keymat <- keymat[match(colnames(mymat), keymat$Sample_ID),]
colnames(mymat) <- keymat$Description
EQ-Vary_R2 <- mymat
save(EQ-Vary_R2, file="RDATA/trimmedExpCountMat_genes_fromHPC_EQ-Vary_R2.RDATA")


## Formatting EQ-75pcnt
mymat <- read.table("RSEMDATA/EQ-75pcnt_R2/exprmat_EQ-75pcnt_R2_trimmed_codingRef.txt", header=T, stringsAsFactors=F)
colnames(mymat) <- lapply(strsplit(colnames(mymat), "_"), function(x) paste0(x[2], "_", x[3]))
dim(mymat)
keymat = read.csv("FASTQ/Ln_1_325587f4c3df848b_EQ-75pcnt.csv", header=T, skip=1, stringsAsFactors=F)
keymat <- keymat[match(colnames(mymat), keymat$Sample_ID),]
colnames(mymat) <- keymat$Description
EQ_75pcnt_R2 <- mymat
save(EQ_75pcnt_R2, file="RDATA/trimmedExpCountMat_genes_fromHPC_EQ-75pcnt_R2.RDATA")



