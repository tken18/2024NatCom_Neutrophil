#!/bin/Rscript

library('Seurat')
library('SeuratData')
library('patchwork')
library('ggplot2')
library('cowplot')
library('dplyr')
library('tidyverse')
library('DoubletFinder')

setwd('~/GSE216020_RAW/')
lf = list.dirs(path = './', recursive = FALSE, full.names = FALSE) #The directory name of each sample should be set as the library name in the SRA run table.
meta.data <- read.csv('./SraRunTable.csv', header = T, sep = ',')
meta.data <- meta.data[!duplicated(meta.data$Library.Name),]

for(i in 2:length(lf)){
  
  count <- Read10X(paste0('./', lf[i]))
  data <- CreateSeuratObject(counts = count, min.cells = 5, min.features = 200, project = "Covid19")
  data$type <- meta.data[meta.data$Library.Name %in% lf[i], 'Genotype']
  data$time <- meta.data[meta.data$Library.Name %in% lf[i], 'Time']
  
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

  filtered.data <- subset(data, subset = percent.mt < 20)
  filtered.data <- NormalizeData(filtered.data)
  filtered.data <- FindVariableFeatures(filtered.data, selection.method = "vst", nfeatures = 2000)
  filtered.data <- ScaleData(filtered.data)
  filtered.data <- RunPCA(filtered.data)
  filtered.data <- RunUMAP(filtered.data, dims = 1:10)
  
  sweep.res.list <- paramSweep_v3(filtered.data, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  
  nExp_poi <- round(0.075*nrow(filtered.data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  filtered.data <- doubletFinder_v3(filtered.data, PCs = 1:10, pN = 0.25, 
                                    pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                                    nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  pdf(paste0('./Fig/', lf[i], ".dimplot.pdf"), width = 6, height = 5)
  print(DimPlot(filtered.data, group.by = paste0('DF.classifications_0.25_', as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), '_', nExp_poi)))
  dev.off()
  
  filtered.data <- filtered.data[,filtered.data@meta.data[[paste0('DF.classifications_0.25_', as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), '_', nExp_poi)]] %in% 'Singlet']

  filtered.data <- NormalizeData(filtered.data, normalization.method = 'RC')
  filtered.data <- ScaleData(filtered.data)
  filtered.data <- FindVariableFeatures(filtered.data, selection.method = "vst", nfeatures = 3000)
  assign(lf[i], filtered.data)

}

anchors <- FindIntegrationAnchors(object.list = list(GSM6656081, GSM6656082, GSM6656083, GSM6656084, GSM6656085, 
                                                     GSM6656086, GSM6656087, GSM6656088, GSM6656089, GSM6656090, 
                                                     GSM6656091, GSM6656092, GSM6656093, GSM6656094, GSM6656095, 
                                                     GSM6656096, GSM6656097, GSM6656098, GSM6656099, GSM6656100, 
                                                     GSM6656101, GSM6656102, GSM6656103, GSM6656104), dims = 1:20)

combined <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
saveRDS(combined, 'combined.rds')
