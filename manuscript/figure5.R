#Use to show the impact of clustering resolution on scMRMA and Uniform (Fig.5).


library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scMRMA)

uniformR_anno <- function(CellType,scdata,logcounts,p,selfClusters,nhvg,res){
  marker <- CellType
  
  if(! is.null(selfClusters)){
    clusters <- selfClusters
  }else{
    hvg <- get_variableGene(scdata)
    hvg.gene <- rownames(hvg)[order(hvg$variance.standardized,decreasing = T)][1:nhvg]
    data.hvg <- logcounts[hvg.gene,]
    scaledata <- t(scale(t(as.matrix(data.hvg))))
    
    set.seed(42)
    pca <- irlba::irlba(t(scaledata),nv=50)
    cell.pca <- pca$u %*% diag(pca$d)
    rownames(cell.pca) <- colnames(scaledata)
    colnames(cell.pca) <- paste0("PC_", 1:50)
    
    knnMatrix <- RANN::nn2(cell.pca,cell.pca,k=20,searchtype = "standard")[[1]]
    snn <- getSNN(knnMatrix,1/15)
    rownames(snn) <- rownames(cell.pca)
    colnames(snn) <- rownames(cell.pca)
    
    clusters <- RunModularityClusteringCpp(
      SNN =snn ,modularityFunction = 1,resolution = res,algorithm = 1,nRandomStarts = 10,
      nIterations = 10,randomSeed = 0,printOutput = F
    )
    names(clusters) <- rownames(cell.pca)
    clusters <- remove_singleton(clusters)
    clusters <- as.factor(clusters)
  }
  #print(table(clusters))
  newCluster <- ori_annotation(marker,clusters,i=ncol(marker),sub = logcounts,p)
  newClusterIDs <- newCluster$newClusterIDs
  layer <- as.data.frame(plyr::mapvalues(clusters,0:(length(newClusterIDs)-1),newClusterIDs))
  colnames(layer) <- "UniformR"
  
  uniformR <- list()
  uniformR$layer <- layer
  uniformR$annoDetails <- newCluster$annoDetails_tem
  uniformR$annoDetails$clusters <- paste("Root",uniformR$annoDetails$clusters,sep = "_")
  return(uniformR)
}

##Mouse brain dataset GSM3580745
brain <- readRDS("mouseBrain.rds")
anno.brain <- scMRMA(brain,species = "Mm")

#Uniform with low resolution (0.6)
CellType <- get_celltype("Mm","panglaodb")
anno.low <- uniformR_anno(CellType,
                          scdata = GetAssayData(brain,slot = "count"),
                          logcounts  = GetAssayData(brain,slot = "data"),
                          p=0.05, selfClusters=NULL, nhvg = 2000, res = 0.6)
brain$uni0.6 <- as.character(anno.low$layer$UniformR)
brain$uni0.6[ brain$uni0.6 %in% c("Smooth muscle cells")] <- "SMCs"
brain$uni0.6[ brain$uni0.6 %in% c("Neural stem/precursor cells")] <- "NPCs"
brain$uni0.6[ brain$uni0.6 %in% c("Oligodendrocyte progenitor cells")] <- "OPCs"
brain$uni0.6[ brain$uni0.6 %in% c("Erythroid-like and erythroid precursor cells")] <- "Precursor cells"

#Uniform with high resolution (1.2)
anno.high <- uniformR_anno(CellType,
                          scdata = GetAssayData(brain,slot = "count"),
                          logcounts  = GetAssayData(brain,slot = "data"),
                          p=0.05, selfClusters=NULL, nhvg = 2000, res = 1.2)
brain$uni1.2 <- as.character(anno.high$layer$UniformR)
brain$uni1.2[ brain$uni1.2 %in% c("Smooth muscle cells")] <- "SMCs"
brain$uni1.2[ brain$uni1.2 %in% c("Neural stem/precursor cells")] <- "NPCs"
brain$uni1.2[ brain$uni1.2 %in% c("Oligodendrocyte progenitor cells")] <- "OPCs"
brain$uni1.2[ brain$uni1.2 %in% c("Erythroid-like and erythroid precursor cells")] <- "Precursor cells"

##scMRMA
brain$scMRMA <- as.character(anno.brain$multiR$annotationResult$Level4)
brain$scMRMA[ brain$scMRMA %in% c("Smooth muscle cells")] <- "SMCs"
brain$scMRMA[ brain$scMRMA %in% c("Neural stem/precursor cells")] <- "NPCs"
brain$scMRMA[ brain$scMRMA %in% c("Oligodendrocyte progenitor cells")] <- "OPCs"
brain$scMRMA[ brain$scMRMA %in% c("Erythroid-like and erythroid precursor cells")] <- "Precursor cells"

color <- c(brewer.pal(12, "Set3"),brewer.pal(7, "Set1"))
p1 <- DimPlot(brain,reduction = "umap",group.by = "uni0.6",label = TRUE,repel = TRUE,
              label.size = 4.5,cols = color[c(1,19,3,4,15,7:11,13,12,14)])+ggtitle("Uniform Low Resolution")
p2 <- DimPlot(brain,reduction = "umap",group.by = "uni1.2",label = TRUE,repel = TRUE,
              label.size = 4.5,cols = color[c(1,19,3,4,6,15,7:11,13,12,14)])+ggtitle("Uniform High Resolution")
p3 <- DimPlot(brain,reduction = "umap",group.by = "scMRMA",label = TRUE,repel = TRUE,
              label.size = 4.5,cols = color[c(1,19,3:11,13,12,14)])

##Number of cells identified as Macrophages, Microglia, Pericytes and SMCs
ct <- c("Macrophages","Microglia","Pericytes","SMCs")
cellNum <- data.frame(CellType = ct,
                      scMRMA = as.numeric(table(brain$scMRMA)[ct]),
                      Low = as.numeric(table(brain$uni0.6)[ct]),
                      High = as.numeric(table(brain$uni1.2)[ct]))
cellNum[is.na(cellNum)] <- 0
colnames(cellNum)[3:4] <- c("Uniform Low Resolution","Uniform High Resolution")
cellNum.plot <- melt(cellNum)

p4 <-ggplot(cellNum.plot,aes(x=CellType,y=value,fill=factor(variable)))+
  geom_bar(position = "dodge",stat = "identity")+
  ylab("Cell Number")+xlab("")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 20),panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black"),axis.text = element_text(colour = "black",size = 15),
        axis.text.x = element_text(hjust=0.5),axis.title.y = element_text(size = 20),
        legend.position = "top",legend.title  = element_blank(),legend.text = element_text(size = 12),
  )+
  scale_y_continuous(limits = c(0,max(cellNum.plot$value)+100),expand = c(0,0))+
  guides(fill=guide_legend(nrow = 1,byrow = T))+
  geom_text(mapping = aes(label = value),position = position_dodge(0.9),vjust=-0.3,size=5,hjust=0.5)

Rmisc::multiplot(p1,p2,p3,p4,layout=matrix(c(1,2,3,4),nrow=2,byrow = T))
