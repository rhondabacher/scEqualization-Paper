setwd("EqualizationManuscript")

load("RDATA/trimmedExpCountMat_genes_fromHPC_EQ_Vary.RDATA")
load("RDATA/trimmedExpCountMat_genes_fromHPC_EQ.RDATA")
load("RDATA/trimmedExpCountMat_genes_fromHPC_EQ_75pcnt.RDATA")

load("RDATA/trimmedExpCountMat_genes_fromHPC_EQ_Vary_R2.RDATA")
load("RDATA/trimmedExpCountMat_genes_fromHPC_EQ_R2.RDATA")
load("RDATA/trimmedExpCountMat_genes_fromHPC_EQ_75pcnt_R2.RDATA")

load("RDATA/trimmedExpCountMat_genes_fromHPC_unEQ_TB.RDATA")
load("RDATA/trimmedExpCountMat_genes_fromHPC_unEQ_EC.RDATA")


## Handle the rownames with gene names appended.
genes <- rownames(EQ_Vary_R2)
ensg <- genes
genes <- strsplit(genes, "_")
genes <- sapply(genes, function(x) {paste(x[-1], collapse="_")})

genes <- data.frame(ensg=ensg, use.g=genes, stringsAsFactors=F)

unEQ_EC <- data.frame(ensg = rownames(unEQ_EC), unEQ_EC, stringsAsFactors=F)
unEQ_EC <- merge(genes, unEQ_EC, by = "ensg")
unEQ_EC = aggregate(.~ use.g, data=unEQ_EC[,-1], FUN=sum)
rownames(unEQ_EC) <- unEQ_EC[,1]
unEQ_EC <- unEQ_EC[,-1]

unEQ_TB <- data.frame(ensg = rownames(unEQ_TB), unEQ_TB, stringsAsFactors=F)
unEQ_TB <- merge(genes, unEQ_TB, by = "ensg")
unEQ_TB = aggregate(.~ use.g, data=unEQ_TB[,-1], FUN=sum)
rownames(unEQ_TB) <- unEQ_TB[,1]
unEQ_TB <- unEQ_TB[,-1]

EQ <- data.frame(ensg = rownames(EQ), EQ, stringsAsFactors=F)
EQ <- merge(genes, EQ, by = "ensg")
EQ = aggregate(.~ use.g, data=EQ[,-1], FUN=sum)
rownames(EQ) <- EQ[,1]
EQ <- EQ[,-1]

EQ_R2 <- data.frame(ensg = rownames(EQ_R2), EQ_R2, stringsAsFactors=F)
EQ_R2 <- merge(genes, EQ_R2, by = "ensg")
EQ_R2 = aggregate(.~ use.g, data=EQ_R2[,-1], FUN=sum)
rownames(EQ_R2) <- EQ_R2[,1]
EQ_R2 <- EQ_R2[,-1]

EQ_Vary <- data.frame(ensg = rownames(EQ_Vary), EQ_Vary, stringsAsFactors=F)
EQ_Vary <- merge(genes, EQ_Vary, by = "ensg")
EQ_Vary = aggregate(.~ use.g, data=EQ_Vary[,-1], FUN=sum)
rownames(EQ_Vary) <- EQ_Vary[,1]
EQ_Vary <- EQ_Vary[,-1]

EQ_Vary_R2 <- data.frame(ensg = rownames(EQ_Vary_R2), EQ_Vary_R2, stringsAsFactors=F)
EQ_Vary_R2 <- merge(genes, EQ_Vary_R2, by = "ensg")
EQ_Vary_R2 = aggregate(.~ use.g, data=EQ_Vary_R2[,-1], FUN=sum)
rownames(EQ_Vary_R2) <- EQ_Vary_R2[,1]
EQ_Vary_R2 <- EQ_Vary_R2[,-1]

EQ_75pcnt <- data.frame(ensg = rownames(EQ_75pcnt), EQ_75pcnt, stringsAsFactors=F)
EQ_75pcnt <- merge(genes, EQ_75pcnt, by = "ensg")
EQ_75pcnt = aggregate(.~ use.g, data=EQ_75pcnt[,-1], FUN=sum)
rownames(EQ_75pcnt) <- EQ_75pcnt[,1]
EQ_75pcnt <- EQ_75pcnt[,-1]

EQ_75pcnt_R2 <- data.frame(ensg = rownames(EQ_75pcnt_R2), EQ_75pcnt_R2, stringsAsFactors=F)
EQ_75pcnt_R2 <- merge(genes, EQ_75pcnt_R2, by = "ensg")
EQ_75pcnt_R2 = aggregate(.~ use.g, data=EQ_75pcnt_R2[,-1], FUN=sum)
rownames(EQ_75pcnt_R2) <- EQ_75pcnt_R2[,1]
EQ_75pcnt_R2 <- EQ_75pcnt_R2[,-1]

# These should be the same
ec.cells <- intersect(colnames(unEQ_EC), colnames(EQ_Vary))
tb.cells <- intersect(colnames(unEQ_TB), colnames(EQ_Vary))

length(ec.cells)
length(tb.cells)

unEQ_TB <- data.matrix(unEQ_TB)
unEQ_EC <- data.matrix(unEQ_EC)

EQr1 <- data.matrix(EQ)
EQr2 <- data.matrix(EQ_R2)
EQ <- EQr1[rownames(EQr1), colnames(EQr1)] + EQ_R2[rownames(EQr1), colnames(EQr1)]

EQ_Varyr1 <- data.matrix(EQ_Vary)
EQ_Varyr2 <- data.matrix(EQ_Vary_R2)
EQ_Vary <- EQ_Varyr1[rownames(EQ_Varyr1), colnames(EQ_Varyr1)] + EQ_Vary_R2[rownames(EQ_Varyr1), colnames(EQ_Varyr1)]

EQ_75pcntr1 <- data.matrix(EQ_75pcnt)
EQ_75pcntr2 <- data.matrix(EQ_75pcnt_R2)
EQ_75pcnt <- EQ_75pcntr1[rownames(EQ_75pcntr1), colnames(EQ_75pcntr1)] + EQ_75pcnt_R2[rownames(EQ_75pcntr1), colnames(EQ_75pcntr1)]


rm(EQ_Vary_R2)
rm(EQ_R2)
rm(EQ_75pcnt_R2)

save.image(file="RDATA/ALL_exprData_trimmed_genes.RData")

