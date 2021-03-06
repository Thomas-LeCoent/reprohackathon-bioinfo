#!/usr/bin/env Rscript
rm(list=ls())


#################### install packages ######################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
BiocManager::install("EnsDb.Hsapiens.v79")

###########################################################

library("DESeq2")
samples<-read.table(file ="${count}",header = TRUE, sep='\t', )
#head(samples)
newsample=subset(samples, select=-c(Chr, Start, End, Strand, Length))
newsample <- newsample[, order(colnames(newsample))]#order column bam files
metaDat=read.csv("${des}", header=TRUE, sep=";")

#DESeq
dds <- DESeqDataSetFromMatrix(countData=newsample, 
                              colData=metaDat, 
                              design=~mutation, tidy = TRUE)

dds <- DESeq(dds)

res <- results(dds, name="mutation_WT_vs_R625C", alpha = 0.05)




library("EnsDb.Hsapiens.v79")
geneids <- rownames(res[which(res\$padj < 0.05),])
genes <- ensembldb::select(EnsDb.Hsapiens.v79, keys = geneids, keytype = "GENEID", column = c("GENENAME"))

res.df <- data.frame(cbind(rownames(res[which(res\$padj < 0.05),]), res[which(res\$"padj" < 0.05),]))
names(res.df)[names(res.df) == colnames(res.df)[[1]]] <- "GENEID"
results <- merge(genes, res.df, all.y = T)


write.csv(results, file= "results.txt", row.names = F)

#table(res\$padj<0.05)[2]

###################Volcanoplot Graphic#############################

library("EnhancedVolcano")
pdf(paste0("volcano.pdf"))


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

dev.off()
