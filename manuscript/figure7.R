#Use to show the advantage of iterative clustering of scMRMA (Fig.6).

library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scMRMA)
library(ggpubr)

##Fig7A TcellAI hierarchical structure
#Generate a Database.html file under current folder
CellType <- get_celltype(species="Hs",db="TcellAI")
databaseVisual(celltype = CellType)

##Fig7B
scc <- readRDS("SCC.rds")
anno.scc <- scMRMA(scc,species = "Hs",db = "TcellAI")
scc$scMRMA <- as.character(anno.scc$multiR$annotationResult$Level2)
scc$Treatment <-  Hmisc::capitalize(meta.scc[colnames(scc),2])
scc$Treatment <- factor(scc$Treatment,levels = c("Pre","Post"))

color <- brewer.pal(8, "Accent")
p1 <- DimPlot(scc,reduction = "umap",group.by = "scMRMA",label = TRUE,repel = TRUE,
              label.size = 5,cols = color[c(1:3,5:7)])
color <- brewer.pal(8, "Set1")
p2 <- DimPlot(scc,reduction = "umap",group.by = "Treatment",label = TRUE,repel = TRUE,
              label.size = 5,cols = color[c(1:2)])

Rmisc::multiplot(p1,p2,layout=matrix(c(1,2),nrow=1,byrow = T))

##Fig7C
#Celltype proportion change pre and post treatment.
scc.anno.all <- data.frame(patient=scc$sample,treatment=scc$Treatment,cluster=scc$ct.ori,
                           scMRMA=scc$scMRMA)

scc.compare <- matrix(0,nrow = 2*length(unique(scc.anno.all$patient)),ncol = length(unique(scc.anno.all$scMRMA)))
rownames(scc.compare) <- paste(rep(unique(scc.anno.all$patient),each=2),c("Post","Pre"),sep="_")
colnames(scc.compare) <- sort(unique(scc.anno.all$scMRMA))
m=1
for (i in unique(scc.anno.all$patient)) {
  print(i)
  submeta <- subset(scc.anno.all,patient==i)
  print(table(submeta[,c(2,4)]))
  scc.compare[m:(m+1),colnames(table(submeta[,c(2,4)]))]<- as.matrix(table(submeta[,c(2,4)]))
  m <- m+2
}
scc.compare <- scc.compare/rowSums(scc.compare)

scc.plot <- data.frame("Per"=as.numeric(scc.compare),"Patient"=rep(rep(unique(scc.anno.all$patient),each=2),ncol(scc.compare)),
                       "Treatment"=rep(c("Post","Pre"),ncol(scc.compare)),
                       "Cluster"=rep(colnames(scc.compare),each=nrow(scc.compare)))
scc.plot$newCluster <- paste(scc.plot$Cluster,scc.plot$Treatment,sep = "_")
scc.plot$Treatment <- factor(scc.plot$Treatment,levels = c("Pre","Post"))

ggplot(scc.plot,aes(x=Treatment,y=Per,color=Patient,shape=Treatment))+
  geom_point(size=2)+geom_line( aes(group=Patient) )+ facet_grid(.~Cluster)+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 15),
        legend.text = element_text(colour = "black",size = 15),legend.title =element_text(colour = "black",size = 15),
        strip.text.x = element_text(colour = "black",size = 12),
        legend.position = "right",legend.box="vertical")+
  labs(x = "Treatment", y = "Cell Percentage")

