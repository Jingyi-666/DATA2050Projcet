#BiocManager::install("Rsubread")
#BiocManager::install("GenomicAlignments")
library("Rsubread")
#library("GenomicAlignments")
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(dplyr)

#?featureCounts

### read BAM files from bowtie2 (female) ###
fc1 <- featureCounts(files=c('data/SRR960177.bam', 
                            'data/SRR21065600.bam', 
                            'data/SRR21065599.bam'),
                     annot.ext = 'data/genomic.gtf', 
                     isGTFAnnotationFile = TRUE, 
                     isPairedEnd=TRUE)
fc2 <- featureCounts(files=c('data/SRR960178.bam'), 
                     annot.ext = 'data/genomic.gtf', 
                     isGTFAnnotationFile = TRUE)

countPaired <- as.data.frame(fc1$counts)
countSingle <- as.data.frame(fc2$counts)
fc <- merge(countPaired, countSingle, by = "row.names")
rownames(fc) <- fc[, 1]
fc <- fc[, -1]
head(fc)

sampleTable <- read.csv('data/sample_table.csv', row.names=1)
sampleTable
colnames(fc) <- sampleTable$Run

countData <- as.data.frame(fc)
colData <- as.data.frame(sampleTable)
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=colData, 
                              design=~Sex, 
                              tidy=TRUE)
dds

#filter out all 0 counts
keep <- rowSums(counts(dds) > 0) > 0
dds <- dds[keep,]

### factor levels ###
dds$Sex <- factor(dds$Sex, levels = c("female","male"))


################################################################
############### differential expression analysis ###############

dds <- DESeq(dds)
# specify the contrast
res <- results(dds, name="Sex_male_vs_female")
res


### log fold change shrinkage ###
# LFC useful for visualization
resLFC <- lfcShrink(dds, coef="Sex_male_vs_female", type="apeglm")
resLFC


### p-values ###
resOrdered <- res[order(res$pvalue),]

summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

# by default alpha=0.1
res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)


################################################################
########################### export ###########################

resSig <- subset(resOrdered, padj < 0.05)
resSig

write.csv(as.data.frame(resSig), file="Sex_male_vs_female.csv")


################################################################
######################## Quality Check ########################
library(pheatmap)
library(vsn)
library(hexbin)
library("RColorBrewer")

# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
# Regularized log transformation
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

#this gives log2(n + 1)
ntd <- ?normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))


################################################################
######################## visualization ########################

### MA-plot ###
png(filename="MAplot_male_vs_female.png", width = 720, height = 720)
MAPlot <- plotMA(resLFC, colNonSig = "gray80", colSig = "navy", 
                 colLine = "grey60", main = "MA plot: male vs female")
dev.off()

# interactively identify the row
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

### heatmap ###
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
heat.df <- as.data.frame(colData(dds)[,c("Sex")])

ntd_heatmap <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                        cluster_cols=FALSE, main="ntd")

vsd_heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                        cluster_cols=FALSE, main="vsd")

rld_heatmap <- pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                        cluster_cols=FALSE, main="rld")

combined_plot <- grid.arrange(grobs = list(ntd_heatmap$gtable, 
                                           vsd_heatmap$gtable, 
                                           rld_heatmap$gtable), 
                              nrow = 1)
ggsave("combined_heatmaps.png", plot = combined_plot)

### volcano plot ###
res_vol <- res %>% data.frame() %>% 
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 1)
summary(res_vol$threshold_OE)
volcano_plot <- 
  ggplot(res_vol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold_OE)) +
  ggtitle("Volcano plot: male vs female") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
ggsave("volcano_plot_male_vs_female.png", plot = volcano_plot)
