#if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("airway")
#BiocManager::install("ggplot2")
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("dplyr")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("AnnotationDbi")
BiocManager::install("genefilter")
install.packages("pheatmap")
install.packages("PoiClaClu")

library(DESeq2)
library(airway)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(ggplot2)
library(genefilter)
library(EnhancedVolcano)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)


## More info https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


# load data from "airway" dataset and create the DESeq object
data("airway")
se <- airway
colData(se)
dds <- DESeqDataSet(se, design = ~ cell + dex)


# plot heatmap of Poisson distances between samples
# use Poisson distance for raw (non-normalized) count data
# use Euclidean distance for data normalized by regularized-logarithm transformation (rlog) or variance stablization transfromation (vst)
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd)
rownames(samplePoisDistMatrix) <- paste(dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


# perform differential gene expression analysis
dds <- DESeq(dds)

# obtain the results
results <- results(dds)

# plot average expression versus log2 fold change - points are colored blue if Padj < 0.1
plotMA(results, ylim = c(-4, 4))


# plot histogram of P-values
hist(results$pvalue, breaks=20, col="grey50", border="white" )

# plot histogram of P-values - improved version by filtering out genes with very low expression levels
hist(results$pvalue[results$baseMean > 1], breaks = 20, col = "grey50", border = "white")



## add gene annotation to results table
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(results), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")
results = cbind(ENSEMBL = rownames(results), results)
anno_results <- left_join(as.data.frame(results), anno )
head(anno_results) 



# volcano plot

#Default cutoffs are log2FC > |2| and adjusted P-value < 0.05
EnhancedVolcano(anno_results, lab = anno_results$SYMBOL, x = "log2FoldChange", y="padj")

#Add custom log2FC and adjusted P-value cutoffs and size of points and labels
EnhancedVolcano(anno_results, lab = anno_results$SYMBOL, x = "log2FoldChange", y="padj", pCutoff = 10e-5, FCcutoff = 2, pointSize = 1.5, labSize = 3.0, title = "Untreated versus treated")

#Adjust axis limits
EnhancedVolcano(anno_results, lab = anno_results$SYMBOL, x = "log2FoldChange", y="padj", xlim = c(-5, 5), ylim = c(0, -log10(10e-10)), title = "Untreated versus treated")

#Modify border and remove gridlines
EnhancedVolcano(anno_results, lab = anno_results$SYMBOL, x = "log2FoldChange", y="padj", border = "full", borderWidth = 1.5, borderColour = "black", gridlines.major = FALSE, gridlines.minor = FALSE, title = "Untreated versus treated")


# perform regularized-logarithm transformation (rlog) on the data
rld <- rlog(dds)


# plot Principal Component Analysis
plotPCA(rld, intgroup = c("dex"))
plotPCA(rld, intgroup = c("dex", "cell"))



# subset genes 
resultsSig <- anno_results[which(anno_results$padj < 0.01 & abs(anno_results$log2FoldChange) >= 1 & anno_results$baseMean >= 20), ]

# print DE genes with strongest downregulation (head) and upregulation (tail)
head(resultsSig[order(resultsSig$log2FoldChange ), ] )
tail(resultsSig[order(resultsSig$log2FoldChange ), ] )

# plot expression of individual genes
# gene with largest positive log2 Fold Change 
plotCounts(dds, gene=which.max(anno_results$log2FoldChange), intgroup="dex")

# specific gene of interest
plotCounts(dds, gene="ENSG00000127954", intgroup="dex")



## EXERCISE 1 ##

# extract gene counts for a gene of interest
geneCounts <- plotCounts(dds, gene = "ENSG00000127954", intgroup = c("dex","cell"), returnData = TRUE)
geneCounts

# use ggplot2 to create a custom plot in which data points from the same cell line have the same color






# plot heatmap of all DE genes
mat <- assay(rld)
idx <- resultsSig$ENSEMBL
DEgenes <- mat[idx,]

annotation <- as.data.frame(colData(rld)[, c("cell","dex")])
pheatmap(DEgenes, scale = "row", show_rownames = FALSE, clustering_distance_rows = "correlation", annotation_col = annotation, main="Differentially Expressed genes")

# plot heatmap of top 20 DE genes
orderedSig <- resultsSig[order(resultsSig$padj), ]
id1 <- orderedSig$ENSEMBL
id2 <- orderedSig$SYMBOL
DE <- mat[id1,]
rownames(topDE) <- id2
top20DE <- head(topDE, n=20)

DEgenesTop <- mat[topDE,]
pheatmap(top20DE, scale = "row", clustering_distance_rows = "correlation", annotation_col = annotation, main="Top 20 Differentially Expressed genes")




## EXERCISE 2 ##

# generate these different plots for the selected RNA-seq data set (.rda file provided)

load(file = "ERP010786.rda")
# Rownames are the gene IDs. They contain a version. This should be removed.
rownames(se)
# Removes gene version
rownames(se) <- gsub(rownames(se), pattern = "\\..+", replacement = "")  
rownames(se)  # Gene ID is now clean



## HOMEWORK ##

# Find your own study on recount2 to analyze and create visualizations 