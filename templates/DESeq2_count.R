#!/usr/bin/env Rscript
rm(list=ls())


#################### install packages ######################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

###########################################################

library("DESeq2")
#order columns
col_order=c("Geneid", "SRR628582.bam", "SRR628583.bam", "SRR628584.bam", "SRR628585.bam", "SRR628586.bam", "SRR628587.bam", "SRR628588.bam", "SRR628589.bam")
samples<-read.table(file ="${count}",header = TRUE, sep='\t', )
#head(samples)
newsample=subset(samples, select=-c(Chr, Start, End, Strand, Length))
newsample <- newsample[, col_order]
metaDat=read.csv("${des}", header=TRUE, sep=";")

#DESeq
dds <- DESeqDataSetFromMatrix(countData=newsample, 
                              colData=metaDat, 
                              design=~mutation, tidy = TRUE)

dds <- DESeq(dds)

res <- results(dds, name="mutation_WT_vs_R625C", alpha = 0.05)

#table(res\$padj<0.05)[2]


library("EnhancedVolcano")
pdf(paste0("volcano.pdf"))


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

dev.off()

nbgenes=length(res\$padj[which(res\$padj<0.05)])
cat(nbgenes, file="results.txt")
