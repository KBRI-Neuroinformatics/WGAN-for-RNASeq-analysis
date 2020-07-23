source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("srnadiff")

sessionInfo()

library(DESeq2)
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)
library(biomaRt)
library(genefilter)

library(gplots)
library(RColorBrewer)
library(dplyr)

rm(list=ls())
ls()

dircsv <- '/fstorage_data/mkcheon/RNA-seq/GSE104775_5xFAD_TREM2_Cortex_72/Hisat2_GRCm38.93_with_Known_dta-cuff/DESeq2_htseq/'

countData <- as.matrix(read.csv(paste(dircsv,"Raw_counts_from_HTSeq_w_HumanGenes.csv",sep=""),row.names="gene_id"))
colData <- read.csv(paste(dircsv,"pheno.csv",sep=""),row.names=1)
tail(countData,3)

all(rownames(colData) %in% colnames(countData))
countData <- countData[,rownames(colData)]
all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ Pheno)
dds <- DESeq(dds)

############### Ensembl ID & Gene Name ######################

ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
lAe <- listAttributes(ensembl)

ensemblID <- dds@rowRanges@partitioning@NAMES
ensemblID_GeneName <- getBM(attributes=c("ensembl_gene_id","external_gene_name","description","gene_biotype"),
                            filters="ensembl_gene_id",values=ensemblID,mart= ensembl)

ensemblID_GeneName_w_HG <- rbind(ensemblID_GeneName,c("ENSG00000080815","hPsen1","",""),
                                 c("ENSG00000095970","hTrem2","",""),c("ENSG00000142192","hAPP","",""))


############## RLog Transform ########################

rld <- rlog(dds)
ddsN <- counts(dds,normalized=TRUE)
write.csv( assay(rld),file=paste(dircsv,"/GSE104775_RLD_w_HG.csv",sep=""))
write.csv( ddsN,file=paste(dircsv,"/GSE104775_DDSN_w_HG.csv",sep=""))

sampleDists <- stats::dist(t(assay(rld)),method="euclidean")
sampleDistMatrix <- as.matrix(sampleDists)


###### DEG Analysis for 7M #####################
res <- results(dds,contrast = c("Pheno","AD7M","WT7M"))
#res005 <- res[which(res$padj < 0.05),]
#res005M <- as.matrix(res005[order(res005$padj,decreasing=FALSE),])
#write.csv(res005M,file=paste(dircsv,"/AD7M_over_TA7M_q005.csv",sep=""))

###### Add basemean for AD, NC #################
baseMeanPerLvl <- sapply(c("AD7M","WT7M"), 
                 function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds@colData@listData$Pheno == lvl] ) )
baseSdsPerLvl <- sapply(c("AD7M","WT7M"), 
                         function(lvl) rowSds(counts(dds,normalized=TRUE)[,dds@colData@listData$Pheno == lvl] ) )

resv <- cbind(as.data.frame(baseMeanPerLvl), res)

resv005 <- resv[which(resv$padj <0.05),]
resv005$gene_name <- ensemblID_GeneName_w_HG$external_gene_name[
                      match(as.list(resv005@rownames),ensemblID_GeneName_w_HG$ensembl_gene_id)]
resv005 <- resv005[,c(9,1:8)]

resv005M <- as.matrix(resv005[order(resv005$padj,decreasing = FALSE),])
write.csv(resv005M,file=paste(dircsv,"/AD7M_over_WT7M_q005_PMean.csv",sep=""))

resv005_bM100 <- resv005[which(resv005$baseMean > 100),]
resv005_bM100M <- as.matrix(resv005_bM100[order(resv005_bM100$padj,decreasing = FALSE),])
write.csv(resv005_bM100M,file=paste(dircsv,"AD7M_over_WT7M_q005_baseMean100_PMean.csv",sep=""))


################################################


