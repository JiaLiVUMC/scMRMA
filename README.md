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
