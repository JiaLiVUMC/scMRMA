# scMRMA: single cell Multi-Resolution Marker-based Annotatation Algorithm

<p align="center">
  <img width="900"  src="https://github.com/JiaLiVUMC/scMRMA/blob/main/overview_scMRMA.png">
</p>

## Installation

scMRMA R package can be easily installed from Github using devtools:  

```
devtools::install_github("JiaLiVUMC/scMRMA")
```
some users might have issues when installing scMRMA package due to the version of C++, please check possible solution through [this website](https://teuder.github.io/rcpp4everyone_en/020_install.html)

## Usage

After installing scMRMA, use following codes to run examples

```R
library(scMRMA)
load(system.file("data", "colon1.Rdata", package = "scMRMA"))
result <- scMRMA(input = colon1,
                 species = "Mm",
                 db = "panglaodb",
                 p = 0.05,
                 normalizedData = F,
                 selfDB = NULL,
                 selfClusters = NULL)
```

`input` Count matrix with genes in row and cells in column. Formats of matrix, dgCMatrix, data.frame, Seurat and SingleCellExperiment object are all acceptable.

`species` Species of cell. Select `"Hs"` (default) or `"Mm"`.

`db` Hierarchical cell type reference database. Select `"panglaodb"` (default) or `"TcellAI"`.

`p` P value cutoff from fisher test for the significant cell type enrichment. Default is `0.05`.

`normalizedData` Use user-provided normalized data. Default is `F` to use default method for normalization.

`selfDB` Use user-provided or modified hierarchical cell type database.

`selfClusters` Use fixed clusters in each level. If provided cluster information, re-clustering will not be performed for intermediate nodes.

__Output__

`result` A list includes annotation results based on multi-resolution and uniform-resolution. 

`result$multiR$annotationResult` A data frame stores scMRMA annotation results for each cell in all reference levels. For example, totally four levels for database `panglaodb`.

`result$multiR$meta` A data frame contains scMRMA cluster, celltype activity score and p value information for each cell in each level.

`result$uniformR$annotationResult` A data frame stores uniform-resolution annotation results for each cell.

`result$uniformR$meta` A data frame contains uniform-resolution cluster, celltype activity score and p value information for each cell in each level.

## Example

__Self-defined database__

```R
# Note: please provide correct format of hierarchical database
# >leaf celltype, root celltype
# GeneA,GeneB,GeneC

CellType <- selfDefinedDatabase(file = system.file("data", "markerExample.txt", package = "scMRMA"))
```

__Add genes to existed database__

```R
# Note: provide the correct format for gene and cell type list. First column includes genes and second column includes cell types in the last level.

genelist <- matrix(c("genea","geneb","Tr1","Microglia"),nrow = 2,byrow = F)
colnames(genelist) <- c("Gene","cell Type")
CellType_new <- addGene(geneCellTypeList = genelist,celltype = CellType)

```

__Hierarchical database visualization__

```R
# Note: it will generate a Database.html file within your current path.

CellType <- get_celltype(species="Hs",db="TcellAI")
databaseVisual(celltype = CellType)
```

__Seurat input__

```R
library(scMRMA)
load(system.file("data", "colon1.Rdata", package = "scMRMA"))

# Creat Seurat object
library(Seurat)
colon1 <- CreateSeuratObject(colon1)
colon1 <- NormalizeData(colon1, verbose=FALSE)
colon1 <- FindVariableFeatures(colon1, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
colon1 <- ScaleData(colon1,features =VariableFeatures(colon1), verbose=FALSE)
colon1 <- RunPCA(colon1,features = VariableFeatures(colon1), npcs = 50,verbose=FALSE)
colon1 <- RunUMAP(colon1, reduction = "pca", dims = 1:50, verbose=FALSE)

# scMRMA annotation
result <-scMRMA(input=colon1, species="Mm")

# UMAP plot
colon1[["scMRMA"]] <- result$multiR$annotationResult[colnames(colon1),ncol(result$multiR$annotationResult)]
DimPlot(colon1,reduction = "umap",group.by = "scMRMA",label = TRUE,repel = TRUE)
colon1[["UniformR"]] <- result$uniformR$annotationResult[colnames(colon1),1]
DimPlot(colon1,reduction = "umap",group.by = "UniformR",label = TRUE,repel = TRUE)
```

<p align="center">
  <img width="900"  src="https://github.com/JiaLiVUMC/scMRMA/blob/main/scMRMA_panglaodb.png">
</p>

```R
# Use user-provided cluster information
# Note: cluster information should be provided as factor

colon1 <- FindNeighbors(colon1,verbose = F)
colon1 <- FindClusters(colon1,resolution = 0.5,verbose=F)
result <- scMRMA(input=colon1, species="Mm",selfClusters=Idents(colon1)
```


