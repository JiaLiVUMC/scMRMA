#Use to show the advantage of scMRMA annotation method from broad to specific (Fig.4).

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(scMRMA)

##Mouse brain dataset GSM3580745
brain <- readRDS("mouseBrain.rds")
anno.brain <- scMRMA(brain,species = "Mm")

##Fig4A
brain$Ref <- brain$ct.ori
brain$Ref[! brain$Ref %in% c("Pericytes","Vascular SMCs")] <- NA
brain$Ref[ brain$Ref %in% c("Vascular SMCs")] <- "SMCs"
brain$scMRMA <- as.character(anno.brain$multiR$annotationResult$Level4)
brain$scMRMA[brain$scMRMA=="Smooth muscle cells"] <- "SMCs"
brain$scMRMA[! brain$scMRMA %in% c("Pericytes","SMCs")] <- NA
brain$Uniform <- as.character(anno.brain$uniformR$annotationResult$UniformR)
brain$Uniform[brain$Uniform=="Smooth muscle cells"] <- "SMCs"
brain$Uniform[! brain$Uniform %in% c("Pericytes","SMCs")] <- NA

color <- brewer.pal(9, "Set1")
p1 <- DimPlot(brain,reduction = "umap",group.by = "Ref",label = TRUE,repel = TRUE,
              label.size = 5,cols = color[c(1,2,9)])+NoLegend()
p2 <- DimPlot(brain,reduction = "umap",group.by = "scMRMA",label = TRUE,repel = TRUE,
              label.size = 5,cols = color[c(1,2,9)])+NoLegend()
p3 <- DimPlot(brain,reduction = "umap",group.by = "Uniform",label = TRUE,repel = TRUE,
              label.size = 5,cols = color[c(1,2,9)])+NoLegend()
Rmisc::multiplot(p1,p2,p3,layout=matrix(c(1,2,3),nrow=1,byrow = T))

##Fig4B
cellNum <- data.frame(CellType = c("Pericytes","SMCs"),
                      Ref=data.frame(table(brain$Ref))[,2],
                      scMRMA=data.frame(table(brain$scMRMA))[,2],
                      Uniform=data.frame(table(brain$Uniform))[,2])
cellNum.plot <- melt(cellNum)
ggplot(cellNum.plot,aes(x=variable,y=value,fill=factor(CellType)))+
  geom_bar(position = "dodge",stat = "identity")+
  ylab("Cell Number")+xlab("")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 20),panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black"),axis.text = element_text(colour = "black",size = 15),
        axis.text.x = element_text(hjust=0.5),axis.title.y = element_text(size = 20),
        legend.position = "top",legend.title  = element_blank(),legend.text = element_text(size = 20),
  )+
  scale_y_continuous(limits = c(0,max(cellNum.plot$value)+100),expand = c(0,0))+
  guides(fill=guide_legend(nrow = 1,byrow = T))+
  scale_fill_manual(values =color[c(1,2)])+
  geom_text(mapping = aes(label = value),position = position_dodge(0.9),vjust=-0.3,size=6,hjust=0.5)

##Fig4C
brain$level1 <- as.character(anno.brain$multiR$annotationResult$Level1)
brain$level2 <- as.character(anno.brain$multiR$annotationResult$Level2)
brain$level3 <-  as.character(anno.brain$multiR$annotationResult$Level3)
brain$level1[! brain$level1 %in% c("Connective tissue cell","Muscle cell")] <- NA
brain$level1[ brain$level1 %in% c("Connective tissue cell")] <- "Connective"
brain$level2[! brain$level2 %in% c("Fibroblast","Pericytes","Smooth muscle cell")] <- NA
brain$level2[ brain$level2 %in% c("Smooth muscle cell")] <- "SMC"
brain$level3[! brain$level3 %in% c("Fibroblasts","Pericytes","Smooth muscle cells")] <- NA
brain$level3[brain$level3 %in% c("Smooth muscle cells")] <- "SMCs"

p5 <- DimPlot(brain,reduction = "umap",group.by = "level1",label = TRUE,repel = TRUE,
              label.size = 5,cols = color[c(4,5)])+ggtitle("Level 1")
p6 <- DimPlot(brain,reduction = "umap",group.by = "level2",label = TRUE,repel = TRUE,
              label.size = 5,cols = color[c(8,1,6)])+ggtitle("Level 2")
p7 <- DimPlot(brain,reduction = "umap",group.by = "level3",label = TRUE,repel = TRUE,
              label.size = 5,cols = color[c(3,1,2)])+ggtitle("Level 3")
Rmisc::multiplot(p5,p6,p7,layout=matrix(c(1,2,3),nrow=1,byrow = T))

##Fig4D
CellType <- get_celltype("Mm","panglaodb")
marker.ct <- unique(CellType[CellType$Level4 %in% c("Smooth muscle cells","Pericytes"),][,1])
marker.ct <- intersect(marker.ct,rownames(brain))
brain <- ScaleData(brain,features = marker.ct)

#scMRMA markers
brain.heatmap <- brain[,colnames(brain) %in% names(brain$scMRMA)[!is.na(brain$scMRMA)]]
Idents(brain.heatmap) <- factor(brain.heatmap$scMRMA,levels = c("Pericytes","SMCs"))
marker <- FindAllMarkers(brain.heatmap,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker <- marker[intersect(rownames(marker),marker.ct),]
marker$cluster <- as.character(marker$cluster)
marker <- marker[order(marker$cluster,-marker$avg_log2FC),]
marker$gene <- factor(marker$gene,levels = marker$gene)
marker.all <- rownames(marker)

p8 <- DoHeatmap(brain.heatmap,features=marker.all,angle=0,hjust = 0.5,slot="scale.data",size = 6,
                group.colors = color[c(1,2)],group.bar = T)+NoLegend()+
  ylab("SMCs                                   Pericytes          ")+ggtitle("scMRMA")+
  scale_y_discrete(labels = rev(levels(marker$gene)) )+
  theme(axis.text.y = element_text(colour = "black"),axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20,hjust = 0))
p9 <- ggplot(marker,aes(x=avg_log2FC,y=gene))+
  geom_bar(position = "dodge",stat = "identity")+
  ylab("")+xlab("log2FC")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text = element_text(colour = "black"),
        axis.text.x = element_text(hjust=1),axis.title.x = element_text(size = 20,vjust = 1),
        panel.border = element_blank(),axis.ticks.y=element_blank(),
  )+
  scale_x_continuous(limits = c(0,max(marker$avg_log2FC)),expand = c(0,0),position = "top")+
  scale_y_discrete(limits = rev(levels(marker$gene)))
Rmisc::multiplot(p8,p9,layout=matrix(c(1,1,1,2),nrow=1,byrow = T))

#Uniform markers
brain.heatmap <- brain[,colnames(brain) %in% names(brain$Uniform)[!is.na(brain$Uniform)]]
Idents(brain.heatmap) <- factor(brain.heatmap$Uniform,levels = c("Pericytes","SMCs"))
marker <- FindAllMarkers(brain.heatmap,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker <- marker[intersect(rownames(marker),marker.ct),]
marker$cluster <- as.character(marker$cluster)
marker <- marker[order(marker$cluster,-marker$avg_log2FC),]
marker$gene <- factor(marker$gene,levels = marker$gene)
marker.all <- rownames(marker)

p10 <- DoHeatmap(brain.heatmap,features=marker.all,angle=0,hjust = 0.5,slot="scale.data",size = 7,
                 group.colors = color[c(1,2)],group.bar = T)+NoLegend()+
  ylab("SMCs                                        Pericytes          ")+ggtitle("Uniform")+
  scale_y_discrete(labels = rev(levels(marker$gene)) )+
  theme(axis.text.y = element_text(colour = "black"),axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 20,hjust = 0))
p11 <- ggplot(marker,aes(x=avg_log2FC,y=gene))+
  geom_bar(position = "dodge",stat = "identity")+
  ylab("")+xlab("log2FC")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text = element_text(colour = "black"),
        axis.text.x = element_text(hjust=1),axis.title.x = element_text(size = 20,vjust = 1),
        panel.border = element_blank(),axis.ticks.y=element_blank(),
  )+
  scale_x_continuous(limits = c(0,max(marker$avg_log2FC)),expand = c(0,0),position = "top")+
  scale_y_discrete(limits = rev(levels(marker$gene)))

Rmisc::multiplot(p10,p11,layout=matrix(c(1,1,1,2),nrow=1,byrow = T))

##Fig4E
#scMRMA
brain.heatmap <- brain[,colnames(brain) %in% names(brain$scMRMA)[!is.na(brain$scMRMA)]]
Idents(brain.heatmap) <- factor(brain.heatmap$scMRMA,levels = c("Pericytes","SMCs"))
p12 <- VlnPlot(brain.heatmap, features = c("Kcnj8"),cols = color[c(1,2)],y.max = 5)+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(angle = 0,hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),
        axis.text =element_text(colour = "black",size = 10),
  )+
  annotate("text",label="LogFC = 2.404",x=2,y=5,size=3.5,color="black")+
  stat_compare_means(comparisons = list(c("Pericytes","SMCs")), label = "p.signif",label.y = 4.3)+
  ylab("Normalized Expression Level")+xlab("")+ggtitle("Kcnj8")
p13 <- VlnPlot(brain.heatmap, features = c("Myh11"),cols = color[c(1,2)],y.max = 5)+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(angle = 0,hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),
        axis.text =element_text(colour = "black",size = 10),
  )+
  annotate("text",label="LogFC = 1.994",x=2,y=5,size=3.5,color="black")+
  stat_compare_means(comparisons = list(c("Pericytes","SMCs")), label = "p.signif",label.y = 4.2)+
  ylab("Normalized Expression Level")+xlab("")+ggtitle("Myh11")
p14 <- VlnPlot(brain.heatmap, features = c("Acta2"),cols = color[c(1,2)],y.max = 7)+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(angle = 0,hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 10),
  )+
  annotate("text",label="LogFC = 2.217",x=2,y=7,size=3.5,color="black")+
  stat_compare_means(comparisons = list(c("Pericytes","SMCs")), label = "p.signif",label.y = 5.6)+
  ylab("Normalized Expression Level")+xlab("")+ggtitle("Acta2")

#Uniform
brain.heatmap <- brain[,colnames(brain) %in% names(brain$Uniform)[!is.na(brain$Uniform)]]
Idents(brain.heatmap) <- factor(brain.heatmap$Uniform,levels = c("Pericytes","SMCs"))
p15 <- VlnPlot(brain.heatmap, features = c("Kcnj8"),cols = color[c(1,2)],y.max = 5)+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(angle = 0,hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),
        axis.text =element_text(colour = "black",size = 10),
  )+
  annotate("text",label="LogFC = 2.226",x=2,y=5,size=3.5,color="black")+
  stat_compare_means(comparisons = list(c("Pericytes","SMCs")), label = "p.signif",label.y = 4.3)+
  ylab("Normalized Expression Level")+xlab("")+ggtitle("Kcnj8")
p16 <- VlnPlot(brain.heatmap, features = c("Myh11"),cols = color[c(1,2)],y.max = 5)+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(angle = 0,hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),
        axis.text =element_text(colour = "black",size = 10),
  )+
  annotate("text",label="LogFC = 0.946",x=2,y=5,size=3.5,color="black")+
  stat_compare_means(comparisons = list(c("Pericytes","SMCs")), label = "p.signif",label.y = 4.2)+
  ylab("Normalized Expression Level")+xlab("")+ggtitle("Myh11")
p17 <- VlnPlot(brain.heatmap, features = c("Acta2"),cols = color[c(1,2)],y.max = 7)+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text.x = element_text(angle = 0,hjust = 0.5,size = 15),
        axis.title = element_text(colour = "black",size = 15),
        axis.text =element_text(colour = "black",size = 10),
  )+
  annotate("text",label="LogFC = 0.789",x=2,y=7,size=3.5,color="black")+
  stat_compare_means(comparisons = list(c("Pericytes","SMCs")), label = "p.signif",label.y = 5.6)+
  ylab("Normalized Expression Level")+xlab("")+ggtitle("Acta2")
Rmisc::multiplot(p12,p13,p14,p15,p16,p17,layout=matrix(c(1,2,3,4,5,6),nrow=1,byrow = T))
