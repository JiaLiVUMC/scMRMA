# scMRMA: single cell Multi-Resolution Marker-based Annotatation Algorithm

<p align="center">
  <img width="900"  src="https://github.com/JiaLiVUMC/scMRMA/blob/main/overview_scMRMA.png">
</p>

## Installation

scMRMA R package can be easily installed from Github using devtools:  

```
devtools::install_github("JiaLiVUMC/scMRMA")
```
some users might have issues when installing scMRMA pacakge due to the version of C++, please check possible solution through [this website](https://teuder.github.io/rcpp4everyone_en/020_install.html)

## Usage

After installing scMRMA, use following codes to run examples

```R
library(scMRMA)
load(system.file("data", "colon1.Rdata", package = "scMRMA"))
result <-scMRMA(input=colon1, species="Mm")
head(result$multiR$annotationResult)
```
