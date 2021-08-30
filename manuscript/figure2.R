##figure2A panglaodb hierarchical structure
#Generate a Database.html file under current folder
library(scMRMA)
CellType <- get_celltype(species="Hs",db="panglaodb")
databaseVisual(celltype = CellType)


##figure2B toy data
#Input toy.exp a seurat object with cluster information
#Annotation were manully added to each level
library(Seurat)
library(ggplot2)
library(RColorBrewer)
#setwd("/Users/lavender/vandy/sample/annotation/manusript/NAR/revision/code")
toy.exp <- readRDS("Fig2_toy.rds")

toy <- as.data.frame(as.numeric(t(t(Idents(toy.exp)))))
rownames(toy) <- colnames(toy.exp)
colnames(toy) <- "Level1"
toy$Level2<- toy$Level1;toy$Level3 <- toy$Level1; toy$Level4 <- toy$Level1
toy$Level1 <- rep("Hematopoietic cell",ncol(toy.exp))
toy$Level2[toy$Level2==5] <- "NK cells"
toy$Level2[toy$Level2==3] <- "B"
toy$Level2[toy$Level2 %in% c(1,2,4)] <- "T"
table(toy$Level2)

toy$Level3[toy$Level3==5] <- "NK cells"
toy$Level3[toy$Level3==3] <- "B cells"
toy$Level3[toy$Level3==4] <- "T cytotoxic cells"
toy$Level3[toy$Level3 %in% c(1,2)] <- "T helper cell"
table(toy$Level3)

toy$Level4[toy$Level4==5] <- "NK cells"
toy$Level4[toy$Level4==3] <- "B cells"
toy$Level4[toy$Level4==4] <- "T cytotoxic cells"
toy$Level4[toy$Level4==1] <- "T helper cells"
toy$Level4[toy$Level4==2] <- "T follicular helper cells"
table(toy$Level4)


color <- brewer.pal(10, "Paired")
toy.exp[["Level1"]] <- toy[colnames(toy.exp),1]
p1<- DimPlot(toy.exp, reduction = "umap",group.by = "Level1",label = TRUE,repel = TRUE,
             cols = color[1])+ggtitle("Level 1")+
  theme(plot.title = element_text(hjust = 0,size = 20),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 15),
        legend.text = element_text(size = 15),
  )
toy.exp[["Level2"]] <- toy[colnames(toy.exp),2]
p2<- DimPlot(toy.exp, reduction = "umap",group.by = "Level2",label = TRUE,repel = TRUE,
             cols = color[c(4,2,6)])+ggtitle("Level 2")+
  theme(plot.title = element_text(hjust = 0,size = 20),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 15),
        legend.text = element_text(size = 15),
  )
toy.exp[["Level3"]] <- toy[colnames(toy.exp),3]
p3<- DimPlot(toy.exp, reduction = "umap",group.by = "Level3",label = TRUE,repel = TRUE,
             cols=color[c(3,2,8,5)])+ggtitle("Level 3")+
  theme(plot.title = element_text(hjust = 0,size = 20),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 15),
        legend.text = element_text(size = 15),
  )
toy.exp[["Level4"]] <- toy[colnames(toy.exp),4]
p4<- DimPlot(toy.exp, reduction = "umap",group.by = "Level4",label = TRUE,repel = TRUE,
             cols = color[c(3,2,8,9,7)])+ggtitle("Level 4")+
  theme(plot.title = element_text(hjust = 0,size = 20),
        axis.title = element_text(colour = "black",size = 15),axis.text =element_text(colour = "black",size = 15),
        legend.text = element_text(size = 15),
  )

Rmisc::multiplot(p1,p2,p3,p4,layout=matrix(c(1,2,3,4),nrow=2,byrow = T))
