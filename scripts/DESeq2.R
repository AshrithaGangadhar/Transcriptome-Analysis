#loading libraries
library(ggplot2) #for creating plots and visualizations
library(pheatmap) #for generating heatmaps
library(DESeq2) #for differential expression analysis
library(dplyr) #for data manipulation
library(RColorBrewer) #for color palettes

#reading data
countdata  <- read.csv("counts.csv",sep=",",header=T, row.names = 1)
metadata <- read.csv("metadata.csv", sep = ",", header=T, row.names = 1)
metadata$Sample<- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$Sample), ]
head(metadata)
all(colnames(count) %in% rownames(metadata))
all(colnames(count) == rownames(metadata))

#creating DESeq2 dataset object from the count data matrix and metadata
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~Group)
ddsMat <- DESeq(ddsMat)

#results extraction
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.01)
summary(results)
mcols(results, use.names = T) 

#variance Stabilization and PCA Plot
ddsMat_rlog <- vst(ddsMat, blind = FALSE)
results_sig <- subset(results, padj < 0.05)
head(results_sig)
plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 500) +
  theme_bw() + 
  geom_point(size = 5) + 
  scale_x_continuous(limits = c(-150, 150)) +
  scale_y_continuous(limits = c(-100, 100)) + 
  ggtitle(label = "Principal Component Analysis (PCA)",
          subtitle = "PCA")

#Heatmap Generation
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:40, ]
mat
annotation_col = data.frame(
  Group = factor(colData(ddsMat_rlog)$Group),
  row.names = colData(ddsMat_rlog)$Sample
)

# Specify colors you want to annotate the columns by.
ann_colors = list(group = c(miR = "lightblue", Negative
 = "yellow"))

# Make Heatmap with pheatmap function.
## See more in documentation for customization
pheatmap(mat = mat,
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255),
         scale = "row", 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         fontsize = 6, 
         cellwidth = , 
)  
#Saving Results
write.table(mat, file = 'topgenes_30.txt',  sep = '\t', col.names = NA)
dev.off()

#volcano plot
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj),
                   lfc = results$log2FoldChange)

data <- na.omit(data)

pvalue = results_sig[which(results_sig[,5]<0.05), ]

upregulated = pvalue[which(pvalue[,2]>0), ]
upregulated
write.csv(upregulated, file = "upregulated.csv")
downregulated = pvalue[which(pvalue[,2]<0), ]
downregulated
write.csv(downregulated, file = "downregulated.csv")

# Color the points which are up or down
## If fold-change > 0 and pvalue > 1.3 (Increased significant)
## If fold-change < 0 and pvalue > 1.3 (Decreased significant)

data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol +  
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("FoldChange value"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values

#MA Plot and Dispersion Estimates
plotMA(results, ylim = c(-5, 5))
plotDispEsts(ddsMat)  

#Saving Significant Results
results_sig <- subset(results, padj < 0.05)
write.csv(results_sig, "results_sig.csv")