#!/bin/Rscript

library('Seurat')
library('SeuratData')
library('ggplot2')
library('ggrepel')

# Read data -----
data <- readRDS('blish_covid.seu.rds')

CA <- read.table('/PATH_TO_LIST/cellAdhesion.txt', header = F)
CA <- intersect(CA$V1, row.names(data))

NET <- read.table('/PATH_TO_LIST/NET.txt', header = F)
NET <- intersect(NET$V1, row.names(data))

# Subset neutrophils -----
Idents(data) <- 'cell.type'
Neu <- subset(data, cell.type == 'Neutrophil')

# DEG analysis for cell adhesion genes -----
i = 1
df <- as.data.frame(Neu@assays$RNA@data[rownames(Neu) %in% CA[i],])
colnames(df) <- 'gene'
df$Status <- Neu@meta.data$Status
df$Ventilated <- Neu@meta.data$Ventilated
test <- pairwise.wilcox.test(df$gene, df$Ventilated, p.adjust.method = "BH")
test.mtx <- as.data.frame(t(c(CA[i], mean(df[df$Ventilated == 'Healthy', 'gene']), 
                              mean(df[df$Ventilated == 'NonVent', 'gene']),
                              mean(df[df$Ventilated == 'Vent', 'gene']),
                              test$p.value[1,1], test$p.value[2,1], test$p.value[2,2])))

colnames(test.mtx) <- c('Gene', 'mean Healthy', 'mean NonVent', 'mean Vent', 
                        'adj.pvalue Healthy vs NonVent', 'adj.pvalue Healthy vs Vent', 'adj.pvalue NonVent vs Vent')

for (i in 2:length(CA)){
  
  df <- as.data.frame(Neu@assays$RNA@data[rownames(Neu) %in% CA[i],])
  colnames(df) <- 'gene'
  df$Status <- Neu@meta.data$Status
  df$Ventilated <- Neu@meta.data$Ventilated
  test <- pairwise.wilcox.test(df$gene, df$Ventilated, p.adjust.method = 'BH')
  tmp <- c(CA[i], mean(df[df$Ventilated == 'Healthy', 'gene']), 
           mean(df[df$Ventilated == 'NonVent', 'gene']),
           mean(df[df$Ventilated == 'Vent', 'gene']),
           test$p.value[1,1], test$p.value[2,1], test$p.value[2,2])
  test.mtx <- rbind(test.mtx, tmp)
  
}

Idents(data) <- 'Ventilated'
DEG <- FindMarkers(subset(data, cell.type == 'Neutrophil'), ident.1 = 'Vent', ident.2 = 'Healthy', logfc.threshold = 0.01)
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

VlnPlot(data, features = c('CD44', 'SELL', 'ICAM3', 'CD93'), assay = 'RNA', slot = 'data', idents = 'Neutrophil', group.by = 'Ventilated')


# DEG analysis for NET genes -----
Idents(data) <- 'cell.type'

Neu <- subset(data, idents = 'Neutrophil')
i = 1
df <- as.data.frame(Neu@assays$RNA@data[rownames(Neu) %in% NET[i],])
colnames(df) <- 'gene'
df$Status <- Neu@meta.data$Status
df$Ventilated <- Neu@meta.data$Ventilated
test <- pairwise.t.test(df$gene, df$Ventilated, p.adjust.method = "bonferroni")
test.mtx <- as.data.frame(t(c(NET[i], mean(df[df$Ventilated == 'Healthy', 'gene']), 
                              mean(df[df$Ventilated == 'NonVent', 'gene']),
                              mean(df[df$Ventilated == 'Vent', 'gene']),
                              test$p.value[1,1], test$p.value[2,1], test$p.value[2,2])))

colnames(test.mtx) <- c('Gene', 'mean Healthy', 'mean NonVent', 'mean Vent', 
                        'adj.pvalue Healthy vs NonVent', 'adj.pvalue Healthy vs Vent', 'adj.pvalue NonVent vs Vent')

for (i in 2:length(NET)){
  
  df <- as.data.frame(Neu@assays$RNA@data[rownames(Neu) %in% NET[i],])
  colnames(df) <- 'gene'
  df$Status <- Neu@meta.data$Status
  df$Ventilated <- Neu@meta.data$Ventilated
  test <- pairwise.t.test(df$gene, df$Ventilated, p.adjust.method = 'BH')
  tmp <- c(NET[i], mean(df[df$Ventilated == 'Healthy', 'gene']), 
           mean(df[df$Ventilated == 'NonVent', 'gene']),
           mean(df[df$Ventilated == 'Vent', 'gene']),
           test$p.value[1,1], test$p.value[2,1], test$p.value[2,2])
  test.mtx <- rbind(test.mtx, tmp)
  
}

VlnPlot(data, features = c('PADI4', 'GSDMD'),  assay = 'RNA', slot = 'data', idents = 'Neutrophil', group.by = 'Ventilated')
