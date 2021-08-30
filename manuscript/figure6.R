#Use to show the advantage of iterative clustering of scMRMA (Fig.6).

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scMRMA)

##Human lung dataset GSM3489182
#Fig6A, Fig6B
lung <- readRDS(file = "Fig6_humanLung.rds")
anno.lung <- scMRMA(lung,species = "Hs")

color <- brewer.pal(9, "Set1")
lung$scMRMA <- as.character(anno.lung$multiR$annotationResult$Level4)
lung$scMRMA[! lung$scMRMA %in% c("Pulmonary alveolar type I cells","Pulmonary alveolar type II cells")] <- NA
p1<-DimPlot(lung,reduction = "umap",group.by = "scMRMA",label = TRUE,repel = TRUE,
            label.size = 5,cols = color[c(1,2,9)])+NoLegend()

lung$uniform <- as.character(anno.lung$uniformR$annotationResult$UniformR)
lung$uniform[! lung$uniform %in% c("Pulmonary alveolar type I cells","Pulmonary alveolar type II cells")] <- NA
p2<-DimPlot(lung,reduction = "umap",group.by = "uniform",label = TRUE,repel = TRUE,
            label.size = 5,cols = color[c(2,9)])+NoLegend()+ggtitle("Uniform")

##Fig6C: scMRMA level1
cluster3 <- data.frame(clu=anno.lung$multiR$meta$Level1_cluster)
rownames(cluster3) <- colnames(lung)
cluster3 <- tidyr::separate(cluster3,col="clu",into = c("Celltype","Clusters"),sep = "_")
lung$clusters <- factor(cluster3[colnames(lung),2],levels = c(1:14))
color <- color <- c(brewer.pal(8, "Set2"),brewer.pal(12,"Paired"))
p3<-DimPlot(lung,reduction = "umap",group.by = "clusters",label = TRUE,repel = TRUE,
            label.size = 5,cols = color[c(1:14)])+
  ggtitle("scMRMA Level1 Clustering")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 15),
        legend.text = element_text(size = 15),
  )

##Fig6D: scMRMA level3
cluster3 <- data.frame(clu=anno.lung$multiR$meta$Level3_cluster)
rownames(cluster3) <- colnames(lung)
cluster3 <- tidyr::separate(cluster3,col="clu",into = c("Celltype","Clusters"),sep = "_")
cluster3[! cluster3$Celltype %in% "Lung epithelial cell",2] <- " "
lung$clusters <- factor(cluster3[colnames(lung),2],levels = c(" ",1:9))
color <- c(brewer.pal(8, "Set2"),brewer.pal(9,"Set1"))
p4<-DimPlot(lung,reduction = "umap",group.by = "clusters",label = TRUE,repel = TRUE,
            label.size = 5,cols = color[c(17,1:7,11:12)])+
  ggtitle("scMRMA Level3 Clustering")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 15),
        legend.text = element_text(size = 15),
  )

##Fig6E, Fig6F: Markers
CellType <- get_celltype(species = "Hs",db = "panglaodb")
marker.ct <- unique(CellType[CellType$Level4 %in% c("Pulmonary alveolar type I cells","Pulmonary alveolar type II cells"),][,1])
marker.ct <- intersect(marker.ct,rownames(lung))#64
lung <- ScaleData(lung,features = marker.ct)

Idents(lung) <- lung$clusters
lung.heatmap <- subset(lung,idents=c(7,9))
lung.heatmap$clusters <- factor(lung.heatmap$clusters,levels = c(7,9))
marker <- FindAllMarkers(lung.heatmap,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker <- marker[intersect(rownames(marker),marker.ct),]
marker <- marker[marker$p_val_adj <0.05,]
marker$cluster <- as.character(marker$cluster)
marker <- marker[order(marker$cluster,-marker$avg_log2FC),]
marker$gene <- factor(marker$gene,levels = marker$gene)
marker.all <- rownames(marker)
color <- c(brewer.pal(8, "Set2"),brewer.pal(9,"Set1"))
p5 <- DoHeatmap(lung.heatmap,features=marker.all,angle=0,hjust = 0.5,slot="data",size = 4,
                group.colors = color[c(7,12)],group.bar = T)+
  scale_y_discrete(labels = rev(levels(marker$gene)) )+
  theme(axis.text.y = element_text(colour = "black"),axis.title.y = element_text(angle = 0))+
  ylab("\n\n\n\n\n\n\nType II\n\n\n\n\n\n\n\nType I")


p6 <- FeaturePlot(lung,features = c("AGER","RTKN2"),slot = "data")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(angle = 0,hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 10),
  )

p7 <- FeaturePlot(lung,features = c("SFTPC","SFTPA1"),slot = "data")+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(angle = 0,hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 10),
  )
Rmisc::multiplot(p1,p2,p3,p4,p5,p6,p7,layout=matrix(c(1,2,1,2,3,4,3,4,5,6,5,7),nrow=6,byrow = T))
