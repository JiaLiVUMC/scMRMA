##This code includes scMRMA annotation for all five datasets

library(Seurat)
library(scMRMA)

#mouse brain
#GSM3580745
#running time: ~2mins
brain <- readRDS("mouseBrain.rds")
brain.scMRMA <- scMRMA(brain,species = "Hs")
brain[["scMRMA"]] <- brain.scMRMA$multiR$annotationResult$Level4
DimPlot(brain,reduction = "umap",group.by = "scMRMA",repel = T,label = T)

#human pancreas
#GSE84133
#running time:~2mins
pancreas <- readRDS("humanPancreas.rds")
pancreas.scMRMA <- scMRMA(pancreas,species = "Hs")
pancreas[["scMRMA"]] <- pancreas.scMRMA$multiR$annotationResult$Level4
DimPlot(pancreas,reduction = "umap",group.by = "scMRMA",repel = T,label = T)

#human PBMC
#GSM2486333
#running time: ~2mins
pbmc <- readRDS("PBMC.rds")
pbmc.scMRMA <- scMRMA(pbmc,species = "Hs")
pbmc[["scMRMA"]] <- pbmc.scMRMA$multiR$annotationResult$Level4
DimPlot(pbmc,reduction = "umap",group.by = "scMRMA",repel = T,label = T)


#human SCC
#GSE123813
#running time: ~7mins
SCC <- readRDS("SCC.rds")
SCC.scMRMA <- scMRMA(SCC,species = "Hs")
SCC[["scMRMA"]] <- SCC.scMRMA$multiR$annotationResult$Level4
DimPlot(SCC,reduction = "umap",group.by = "scMRMA",repel = T,label = T)


#human lung
#GSE130148
#running time: ~4mins
lung <- readRDS("humanLung.rds")
lung.scMRMA <- scMRMA(lung,species = "Hs")
lung[["scMRMA"]] <- lung.scMRMA$multiR$annotationResult$Level4
DimPlot(lung,reduction = "umap",group.by = "scMRMA",repel = T,label = T)