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

## Example

After installing scMRMA, use following codes to run [example](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3580745):

```R
# Note: example will take two minutes.

library(scMRMA)
load(system.file("data", "MouseBrain.Rdata", package = "scMRMA"))
result <- scMRMA(input = Brain_GSM3580745,
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

## Other functions

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

<p align="center">
  <img width="450"  src="https://github.com/JiaLiVUMC/scMRMA/blob/main/Database_TcellAI.png">
</p>

__Incorporate with Seurat__

```R
library(scMRMA)
load(system.file("data", "Brain_GSM3580745.Rdata", package = "scMRMA"))

# Create Seurat object
library(Seurat)
brain <- CreateSeuratObject(Brain_GSM3580745)
brain <- NormalizeData(brain, verbose=FALSE)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
brain <- ScaleData(brain,features =VariableFeatures(brain), verbose=FALSE)
brain <- RunPCA(brain,features = VariableFeatures(brain), npcs = 50,verbose=FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:50, verbose=FALSE)

# scMRMA annotation
result <-scMRMA(input=brain, species="Mm")

# UMAP plot
brain[["scMRMA"]] <- result$multiR$annotationResult[colnames(colon1),ncol(result$multiR$annotationResult)]
DimPlot(brain,reduction = "umap",group.by = "scMRMA",label = TRUE,repel = TRUE)
```

<p align="center">
  <img width="900"  src="https://github.com/JiaLiVUMC/scMRMA/blob/main/scMRMA_panglaodb.png">
</p>

```R
# Use user-provided cluster information
# Note: cluster information should be provided as factor

brain <- FindNeighbors(brain,verbose = F)
brain <- FindClusters(brain,resolution = 0.5,verbose=F)
result <- scMRMA(input=brain, species="Mm",selfClusters=Idents(brain)
```


