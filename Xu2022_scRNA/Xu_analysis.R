#!/bin/Rscript

library('Seurat')
library('SeuratData')
library('patchwork')
library('ggplot2')
library('ggrepel')
library('cowplot')
library('dplyr')
library('wesanderson')
library('yarrr')
library('tidyverse')

# Read data -----
combined <- readRDS('./combined.rds')
DefaultAssay(combined) <- "integrated"
combined <- FindClusters(combined, resolution = 0.2)

CA <- read.table('/PATH_TO_LIST/cellAdhesion.txt', header = F)
CA <- intersect(CA$V1, row.names(data))

NET <- read.table('/PATH_TO_LIST/NET.txt', header = F)
NET <- intersect(NET$V1, row.names(data))


# DEG analysis for cell adhesion genes -----
Neu <- subset(combined, idents = c('0', '2'))
i = 1
df <- as.data.frame(Neu@assays$RNA@data[rownames(Neu) %in% CA[i],])
colnames(df) <- 'gene'
df$type <- Neu@meta.data$type
test <- pairwise.wilcox.test(df$gene, df$type, p.adjust.method = "BH")
test.mtx <- as.data.frame(t(c(CA[i], mean(df[df$type == 'Healthy', 'gene']), 
                              mean(df[df$type == 'Mild', 'gene']),
                              mean(df[df$type == 'Severe', 'gene']),
                              test$p.value[1,1], test$p.value[2,1], test$p.value[2,2])))

colnames(test.mtx) <- c('Gene', 'mean Healthy', 'mean Mild', 'mean Severe', 
                        'adj.pvalue Healthy vs Mild', 'adj.pvalue Healthy vs Severe', 'adj.pvalue Mild vs Severe')

for (i in 2:length(CA)){
  
  df <- as.data.frame(Neu@assays$RNA@data[rownames(Neu) %in% CA[i],])
  colnames(df) <- 'gene'
  df$type <- Neu@meta.data$type
  test <- pairwise.wilcox.test(df$gene, df$type, p.adjust.method = "BH")
  tmp <- c(CA[i], mean(df[df$type == 'Healthy', 'gene']), 
           mean(df[df$type == 'Mild', 'gene']),
           mean(df[df$type == 'Severe', 'gene']),
           test$p.value[1,1], test$p.value[2,1], test$p.value[2,2])
  test.mtx <- rbind(test.mtx, tmp)
  
}

DEG <- FindMarkers(subset(combined, integrated_snn_res.0.2 == '0' | integrated_snn_res.0.2 == '2'), ident.1 = 'Severe', ident.2 = 'Healthy', logfc.threshold = 0.01)
DEG$BH.pvalue <- p.adjust(DEG$p_val, method = 'BH', n = length(DEG$p_val))

g = ggplot(NULL)
g <- g + geom_point(data = DEG[row.names(DEG) %in% CA, ], aes(x = avg_log2FC, y = -log10(BH.pvalue)), col = 'grey40', size = 3)
g <- g + geom_point(data = DEG[row.names(DEG) %in% CA & DEG$BH.pvalue < 0.05, ], 
                    aes(x = avg_log2FC, y = -log10(BH.pvalue)), col = 'purple3', size = 3)
g <- g + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")
#g <- g + geom_vline(xintercept = c(1, -1), linetype = "dashed", color = "gray")
g <- g + geom_text_repel(DEG[DEG$BH.pvalue < 0.05 & row.names(DEG) %in% CA,], 
                         mapping = aes(x = avg_log2FC, y = -log10(BH.pvalue), 
                                       label = row.names(DEG[DEG$BH.pvalue < 0.05 & row.names(DEG) %in% CA,]), size = 7), show.legend = FALSE)
g = g + theme_classic() + xlim(-2, 2)
plot(g)

VlnPlot(combined, features = c('CD44', 'SELL', 'ICAM3', 'CD93'), assay = 'RNA', slot = 'data', idents = c('0', '2'), group.by = 'type', cols = pal2)


# DEG analysis for NET genes -----
Neu <- subset(combined, idents = c('0', '2'))
i = 1
df <- as.data.frame(Neu@assays$RNA@data[NET[i],])
colnames(df) <- 'gene'
df$type <- Neu@meta.data$type
test <- pairwise.t.test(df$gene, df$type, p.adjust.method = "BH")
test.mtx <- as.data.frame(t(c(NET[i], mean(df[df$type == 'Healthy', 'gene']), 
                              mean(df[df$type == 'Mild', 'gene']),
                              mean(df[df$type == 'Severe', 'gene']),
                              test$p.value[1,1], test$p.value[2,1], test$p.value[2,2])))

colnames(test.mtx) <- c('Gene', 'mean Healthy', 'mean Mild', 'mean Severe', 
                        'adj.pvalue Healthy vs Mild', 'adj.pvalue Healthy vs Severe', 'adj.pvalue Mild vs Severe')

for (i in 2:length(NET)){
  
  df <- as.data.frame(Neu@assays$RNA@data[NET[i],])
  colnames(df) <- 'gene'
  df$type <- Neu@meta.data$type
  test <- pairwise.t.test(df$gene, df$type, p.adjust.method = 'BH')
  tmp <- c(NET[i], mean(df[df$type == 'Healthy', 'gene']), 
           mean(df[df$type == 'Mild', 'gene']),
           mean(df[df$type == 'Severe', 'gene']),
           test$p.value[1,1], test$p.value[2,1], test$p.value[2,2])
  test.mtx <- rbind(test.mtx, tmp)
  
}

VlnPlot(combined, features = c('PADI4', 'GSDMD'), assay = 'RNA', slot = 'data', idents = c('0', '2'), group.by = 'type', raster = FALSE)
