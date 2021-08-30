##This code includes all five methods (scMRMA, scCATCH, Garnett, SingleR, CellAssign) for humanPBMC dataset.
#It also includes code for UMAP (Supplementary Fig.S3) and barplot figure (Fig.3)

##PBMC data
#GSM2486333
library(Seurat)
library(scMRMA)
library(scCATCH)
library(garnett)
library(SingleR)
library(org.Hs.eg.db)
library(celldex)
library(cellassign)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

pbmc <- readRDS("PBMC.rds")

#Ref
pbmc$Ref <- pbmc$ct.ori
pbmc$Ref <- gsub("DC","Myeloid",pbmc$Ref)
DimPlot(pbmc,reduction = "umap",group.by = "Ref",repel = T,label = T)

##scMRMA
#running time: ~2mins
pbmc.scMRMA <- scMRMA(pbmc,species = "Hs")
pbmc[["scMRMA"]] <- pbmc.scMRMA$multiR$annotationResult$Level4
DimPlot(pbmc,reduction = "umap",group.by = "scMRMA",repel = T,label = T)

##scCATCH
#running time: ~2mins
pbmc <- FindNeighbors(pbmc,reduction = "pca", dims = 1:30, verbose=FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5)#8 clusters
pbmc.marker <- findmarkergenes(pbmc, species ="Human",match_CellMatch=T,
                                           tissue = c('Blood','Peripheral blood', 'Plasma', 'Serum', 'Umbilical cord blood', 'Venous blood'))
pbmc.ann <- scCATCH(pbmc.marker$clu_markers, 'Human', 
                                tissue = c('Blood','Peripheral blood', 'Plasma', 'Serum', 'Umbilical cord blood', 'Venous blood'))

a <- data.frame(scCATCH=plyr::mapvalues(Idents(pbmc),0:(dim(pbmc.ann)[1]-1),pbmc.ann$cell_type))

pbmc[["scCATCH"]] <- a$scCATCH
DimPlot(pbmc,reduction = "umap",group.by = "scCATCH",repel = T,label = T)

##Garnett
#use pre-trained human PBMC classifier: https://cole-trapnell-lab.github.io/garnett/classifiers/
#running time: ~8s
mat <- GetAssayData(pbmc,slot = "counts")
fdata <- as.data.frame(matrix(1,nrow = nrow(mat),ncol = 1))
rownames(fdata)<- rownames(mat)
pdata <- as.data.frame( Idents(pbmc))
pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
cds <- newCellDataSet(as(mat, "dgCMatrix"), phenoData = pd,featureData = fd)
cds <- estimateSizeFactors(cds)

set.seed(100)
classifier <- readRDS("hsPBMC_20191017.RDS")
system.time(cds <- classify_cells(cds, classifier, db = org.Hs.eg.db,cluster_extend = TRUE,cds_gene_id_type = "SYMBOL"))#8s
pbmc[["Garnett"]] <- pData(cds)[colnames(pbmc),]$cluster_ext_type
DimPlot(pbmc, reduction = "umap",group.by = "Garnett",label = TRUE,repel = TRUE)

##SingleR
#reference: HPCA main label
#running time: ~35s
rawsample <- GetAssayData(pbmc,slot = "counts")
hpca.se <- HumanPrimaryCellAtlasData()
common.hpca <- intersect(rownames(rawsample), rownames(hpca.se))
trained.hpca <- trainSingleR(hpca.se[common.hpca,], labels=hpca.se$label.main)
hpca <- classifySingleR(rawsample[common.hpca,], trained.hpca)
pbmc[["SingleR"]] <- hpca$pruned.labels
DimPlot(pbmc, reduction = "umap",group.by = "SingleR",label = TRUE,repel = TRUE)

##CellAssign
#Use garnett trained classifier file as marker file
#running time: ~15mins
#pbmc.sce <- as.SingleCellExperiment(pbmc)
pbmc_diet <- DietSeurat(pbmc, graphs = "pca")
pbmc.sce <- as.SingleCellExperiment(pbmc_diet)
s <- scran::calculateSumFactors(pbmc.sce)#33s

marker.garnett <- read.table("Garnett_PBMC_marker.txt",header = T,sep = "\t",stringsAsFactors = F)
garnett.celltype <- tapply(marker.garnett$Gene,marker.garnett$GarnettCelltype,list)
marker_mat <- marker_list_to_mat(garnett.celltype,include_other = FALSE)
sce_marker <- pbmc.sce[intersect(rownames(marker_mat), rownames(pbmc.sce)),]
marker_mat_sce <- marker_mat[rownames(sce_marker),]
fit <- cellassign(exprs_obj = sce_marker,marker_gene_info = marker_mat_sce,
                              s=s,learning_rate = 1e-2,shrinkage = T,verbose = F)
pbmc[["CellAssign"]] <- celltypes(fit)
DimPlot(pbmc, reduction = "umap",group.by = "CellAssign",label = TRUE,repel = TRUE)

##figure3 for PBMC dataset
data <- read.table("PBMC_annotation.txt",header = T,sep = "\t",stringsAsFactors = F)
data2 <- melt(data[,c(1,3:7)])
colnames(data2) <- c("Celltype","Algorithm","Number")
data2$log <- log10(data2$Number)
data2[data2$Number==0,]$log<-0
data2$label <- data2$Number
data2$label[data2$label==0] <- ""
data2$hline <- log10(rep(data[,2],(ncol(data)-2)))
data2$Celltype <- factor(data2$Celltype,levels = data[,1])

color <- brewer.pal(8, "Set1")
ggplot(data2,aes(x=Celltype,y=log,fill=factor(Algorithm)))+
  geom_bar(position = "dodge",stat = "identity")+
  ylab("log10(Cell Number)")+xlab("")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 20),panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black"),axis.text = element_text(colour = "black",size = 15),
        axis.text.x = element_text(angle = 45,hjust=1),axis.title.y = element_text(size = 15),
        legend.position = "none",
        #legend.position = "bottom",legend.title  = element_blank(),legend.text = element_text(size = 20),
  )+
  geom_errorbar(aes(ymax=hline, ymin=hline), linetype="dashed",color="red")+
  scale_y_continuous(limits = c(0,max(data2$log)+1),expand = c(0,0))+
  guides(fill=guide_legend(nrow = 1,byrow = T))+
  #scale_fill_brewer(palette = 'Set1')+
  scale_fill_manual(values = color[2:6])+
  geom_text(mapping = aes(label = Number),position = position_dodge(0.9),vjust=0.5,size=4,angle=90,hjust=-0.1)

