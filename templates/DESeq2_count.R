#!/usr/bin/env Rscript
rm(list=ls())


#################### install packages ######################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
BiocManager::install("EnsDb.Hsapiens.v79")

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

library("EnsDb.Hsapiens.v79")

id_and_pvalues <- data.frame(cbind(rownames(res[which(res$padj < 0.05),]), res[which(res\$padj<0.05),]\$padj))
colnames(id_and_pvalues) = c("GENEID", "padj")

geneids <- rownames(res[which(res\$padj < 0.05),])
genes <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneids, keytype = "GENEID", column = c("GENENAME"))

results <- merge(genes, id_and_pvalues, all.y = T)

write.csv(results, file= "results.txt", row.names = F)

#table(res\$padj<0.05)[2]

###########################################################

library("EnhancedVolcano")
pdf(paste0("volcano.pdf"))


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

dev.off()
