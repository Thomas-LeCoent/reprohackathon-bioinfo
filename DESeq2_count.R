###installation de DESeq2, passila, Rsubread

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("pasilla")
# 
# 
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")
# 
rm(list=ls())

library("DESeq2")
#order columns
col_order=c("Geneid", "SRR628582.bam", "SRR628583.bam", "SRR628584.bam", "SRR628585.bam", "SRR628586.bam", "SRR628587.bam", "SRR628588.bam", "SRR628589.bam")
samples<-read.table(file ="counts.txt",header = TRUE, sep='\t', )
#head(samples)
newsample=subset(samples, select=-c(Chr, Start, End, Strand, Length))
newsample <- newsample[, col_order]
metaDat=read.csv("descriptionMutation.txt", header=TRUE, sep=";")

#DESeq
dds <- DESeqDataSetFromMatrix(countData=newsample, 
                              colData=metaDat, 
                              design=~mutation, tidy = TRUE)

dds <- DESeq(dds)

res <- results(dds, name="mutation_WT_vs_R625C", alpha = 0.05)

table(res$padj<0.05)

