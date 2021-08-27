#' scMRMA
#'
#' @description scMRMA single cell annotation

#' @param input
#' @param species
#' @param db
#' @param p
#' @param normalizedData
#' @param selfClusters
#'
#' @return Result

#' @examples
#' library(scMRMA)
#' load(system.file("data", "MouseBrain.Rdata", package = "scMRMA"))
#' result<-scMRMA(input= Brain_GSM3580745, species="Mm")
#' result
#' @import RANN irlba plyr tidyr data.tree networkD3 stats
#' @importFrom devtools session_info
#' @export
#'
#' @description scMRMA single cell annotation
#'
#' @export
scMRMA <- function(input,species,db="panglaodb",p=0.05,
                           normalizedData=F,selfDB=NULL,selfClusters=NULL,k=20){
  
  if(is(input,"SingleCellExperiment")){
    scdata <- counts(input)
  }else if(is(input,"Seurat")){
    if (!require("Seurat",character.only = TRUE)) { stop("Please install Seurat")}
    scdata <- Seurat::GetAssayData(input,slot="counts")
  }else if(is(input,"dgCMatrix")){
    scdata <- input
  }else if(is(input, "matrix")){
    scdata <- Matrix::Matrix(input,sparse=TRUE)
  }else if(is(input,"data.frame")){
    scdata <- Matrix::Matrix(as.matrix(input),sparse=TRUE)
  }else{
    stop("Please provide correct format of input. SingleCellExperiment, Seurat, dgCMatrix, matrix and data.frame are accepted.")
  }
  
  if(nrow(scdata) < 200){stop("Less than 200 genes.")}
  if(nrow(scdata) > 2000){
    nhvg <- 2000
  }else{
    warning("Less than 2000 genes.")
    nhvg <- nrow(scdata)
  }
  
  if(!is.null(selfClusters)){
    if(! length(selfClusters) == ncol(scdata)){
      stop("Provided clusters should have the same cell number as input matrix.")
    }else if(! identical(sort(names(selfClusters)),sort(colnames(scdata)))){
      stop("Cell name for provided clusters and input matrix are not match.")
    }
  }
  
  if(!is.null(selfDB)){
    cat("User-provided cell type database will be used.\n")
    CellType <- selfDB
  }else{
    cat("Pre-defined cell type database",db,"will be used.\n")
    CellType <- get_celltype(species,db)
    #print(unique(CellType[,2]))
  }
  
  if(length(intersect(CellType[,1],rownames(scdata))) ==0){
    stop("No overlapped genes for input matrix and reference database. Please check the structure of database or species.")
  }
  
  if(normalizedData){
    if(is(normalizedData,"SingleCellExperiment")){
      logcounts <- counts(normalizedData)
    }else if(is(normalizedData,"Seurat")){
      if (!require("Seurat",character.only = TRUE)) { stop("Please install Seurat")}
      logcounts <- Seurat::GetAssayData(normalizedData,slot="counts")
    }else if(is(normalizedData,"dgCMatrix")){
      logcounts <- normalizedData
    }else if(is(normalizedData, "matrix")){
      logcounts <- Matrix::Matrix(normalizedData,sparse=TRUE)
    }else if(is(normalizedData,"data.frame")){
      logcounts <- Matrix::Matrix(as.matrix(normalizedData),sparse=TRUE)
    }
  }else{
    logcounts <- LogNorm(scdata)
    rownames(logcounts) <- rownames(scdata)
    colnames(logcounts) <- colnames(scdata)
  }
  
  annoResult <- list()
  anno <- matrix("unknown",nrow = dim(scdata)[2],ncol = (ncol(CellType)))
  rownames(anno)<- colnames(scdata)
  colnames(anno) <- c(colnames(CellType)[2:ncol(CellType)],"extend")
  
  annoDetails <- matrix(NA,nrow = dim(scdata)[2],ncol = (ncol(CellType)-1)*3)
  rownames(annoDetails)<- colnames(scdata)
  colnames(annoDetails)<- paste(rep(colnames(CellType)[2:ncol(CellType)],each=3),
                                rep(c("cluster","celltypeActivity","pValue"),ncol(CellType)-1),sep = "_")
  annoDetails <- as.data.frame(annoDetails,stringsAsFactors = FALSE)
  
  cat("Multi Resolution Annotation Started.","\n")
  for (i in 2:ncol(CellType)){
    cat("Level",i-1,"annotation started.","\n")
    if(i==2){
      marker <- CellType
      if(!is.null(selfClusters)){
        clusters <- selfClusters
      }else{
        clusters <- rough_cluster(scdata = scdata,logcounts = logcounts,nhvg = nhvg,k=k)
      }
      #print(table(clusters))
      newCluster <- ori_annotation(marker,clusters,i,sub = logcounts,p)
      newClusterIDs <- newCluster$newClusterIDs
      #print(newClusterIDs)
      layer <- as.data.frame(plyr::mapvalues(clusters,0:(length(newClusterIDs)-1),newClusterIDs))
      anno[rownames(layer),1] <- as.character(layer[,1])
      anno[,i]<-anno[,i-1]
      annoDetails[rownames(newCluster$annoDetails_tem),1:3] <- newCluster$annoDetails_tem[,1:3]
      #print(head(annoDetails))
    }else{
      for (m in unique(anno[,i-2])){
        if (m != "Unassigned"){
          marker <- subset(CellType,CellType[,i-1]==m)
          if(length(unique(marker[,i])) > 1 ){##stop if there is no subtype
            name <- rownames(anno)[anno[,i-2]==m]
            sub <- scdata[,name]
            subLog <- logcounts[,name]
            if(ncol(sub)>20){
              npcs=50;
              if(ncol(sub)<50){npcs=ncol(sub)}
              if(! is.null(selfClusters)){
                clusters <- selfClusters[name]
                clusters <- factor(clusters)
                clusters <- plyr::mapvalues(clusters,levels(clusters),0:(length(levels(clusters))-1))
              }else{
                clusters <- precise_cluster(rawdata = sub,logdata = subLog,npcs,nhvg,k=k)
              }
              #print(table(clusters))
              #clusters <- precise_cluster(rawdata = sub,logdata = subLog,npcs,nhvg)
              newCluster <- ori_annotation(marker,clusters,i,sub=subLog,p)
              newClusterIDs <- newCluster$newClusterIDs
              layer <- as.data.frame(plyr::mapvalues(clusters,0:(length(newClusterIDs)-1),newClusterIDs))
              anno[rownames(layer),i-1] <- as.character(layer[,1])
              #print(head(newCluster$annoDetails_tem))
              #print(class(annoDetails))
              annoDetails[rownames(newCluster$annoDetails_tem),c(3*i-5,3*i-4,3*i-3)] <- newCluster$annoDetails_tem[,1:3]
            }else{
              cat("single cluster")
              clusters <- t(t(rep(0,ncol(sub))))
              rownames(clusters) <- colnames(sub)
              clusters <- as.factor(clusters[,1])
              newCluster  <- ori_annotation(marker,clusters,i,sub,p)
              newClusterIDs <- newCluster$newClusterIDs
              anno[names(clusters),i-1] <- rep(newClusterIDs,length(clusters))
              annoDetails[rownames(newCluster$annoDetails_tem),c(3*i-5,3*i-4,3*i-3)] <- newCluster$annoDetails_tem[,1:3]
            }
          }
        }
      }
      anno[,i]<-anno[,i-1]
    }
  }
  for (i in 1:nrow(anno)) {
    if(anno[i,ncol(anno)] == "Unassigned"){
      if(which(anno[i,]=="Unassigned")[1]>1){anno[i,which(anno[i,]=="Unassigned")[1]:ncol(anno)]=anno[i,which(anno[i,]=="Unassigned")[1]-1]}
    }
  }
  
  anno <- as.data.frame(anno)
  for (i in seq(from=1, to=ncol(annoDetails), by=3)) {
    if(i == 1){
      annoDetails[,i] <- paste("Root",annoDetails[,i],sep = "_")
    }else{
      annoDetails[,i] <- paste(anno[,(i-1)/3],annoDetails[,i],sep="_")
    }
  }
  annoDetails[annoDetails=="Unassigned_NA"] <- NA
  
  annoResult$multiR$annotationResult <- anno[,1:ncol(anno)-1]
  annoResult$multiR$meta <- annoDetails
  cat("Uniform Resolution Annotation Started.","\n")
  uniformR <- uniformR_anno(CellType,scdata,logcounts,p,selfClusters,nhvg,k)
  annoResult$uniformR$annotationResult <- uniformR$layer
  annoResult$uniformR$meta <- uniformR$annoDetails
  return(annoResult)
}

#' uniform resolution annotation
#'
#' @description uniform resolution annotation
#' 
#' @param CellType reference cell type
#' @param scdata single cell count data
#' @param logcounts single cell normalized data
#' @param p p value cutoff
#' @param selfClusters user-provided cluster information
#'
#' @return uniformR
#' @import RANN irlba plyr tidyr
#'
#' @export
uniformR_anno <- function(CellType,scdata,logcounts,p,selfClusters,nhvg,k){
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
    
    knnMatrix <- RANN::nn2(cell.pca,cell.pca,k=k,searchtype = "standard")[[1]]
    snn <- getSNN(knnMatrix,1/15)
    rownames(snn) <- rownames(cell.pca)
    colnames(snn) <- rownames(cell.pca)
    
    clusters <- RunModularityClusteringCpp(
      SNN =snn ,modularityFunction = 1,resolution = 1,algorithm = 1,nRandomStarts = 10,
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

#' get hierarchical database
#'
#' @description get hierarchical database
#' 
#' @param species choose species
#' @param db hierarchical database

#' @return CellType
#'
#' @export
get_celltype <- function(species,db){
  Species = "Hs"
  if(db == "panglaodb"){
    load(system.file("data", "Human_PanglaoDB.Rdata", package = "scMRMA"))
    if(toupper(species) == "MM"){
      load(system.file("data", "Mouse_PanglaoDB.Rdata", package = "scMRMA"))
    }
    return(CellType)
  }
  if(db == "TcellAI"){
    load(system.file("data", "Human_TcellAI.Rdata", package = "scMRMA"))
    if(toupper(species) == "MM"){
      stop("No TcellAI for mouse.")
    }
  }
  return(CellType)
}

#' get_variableGene
#'
#' @description find variable genes based on normalized data
#' 
#' @param data normalized data
#'
#' @return hvf.info
#' @import Matrix
#'
#' @export
get_variableGene <- function(data,chunk=1000){
  gene_mean <- Matrix::rowMeans(data)
  gene_var <- SparseRowVar2(data,gene_mean,display_progress = F)
  
  hvf.info <- data.frame(mean = gene_mean,variance = gene_var)
  hvf.info$variance.expected <- 0
  hvf.info$variance.standardized <- 0
  fit <- loess(formula = log10(variance) ~ log10(mean),data = hvf.info[hvf.info$variance > 0, ],span = 0.3)
  hvf.info$variance.expected[hvf.info$variance > 0] <- 10 ^ fit$fitted
  hvf.info$variance.standardized <- SparseRowVarStd(data,hvf.info$mean,sd=sqrt(hvf.info$variance.expected),
                                                    vmax = sqrt(ncol(data)),display_progress = F)
  rownames(hvf.info) <- rownames(data)
  return(hvf.info)
}

#' remove singleton
#'
#' @description remove clusters with only one cell

#' @param clusters clusters
#'
#' @return clusters
#'
#' @export
remove_singleton <- function(clusters){
  singleton <- names(which(table(clusters)==1))
  singleton <- intersect(singleton,clusters)
  leftCluster <- setdiff(unique(clusters),singleton)
  
  for (i in singleton){
    set.seed(1)
    closestCluster <- sample(leftCluster,1)
    clusters[names(which(clusters==i))] <- closestCluster
  }
  return(clusters)
}

#' rough cluster
#'
#' @description clustering process for level 1
#' 
#' @param scdata single cell count data
#' @param logcounts normalized data
#'
#' @return clusters
#' @import irlba RANN
#'
#' @export
rough_cluster <- function(scdata,logcounts,nhvg,k){
  hvg <- get_variableGene(scdata)
  hvg.gene <- rownames(hvg)[order(hvg$variance.standardized,decreasing = T)][1:nhvg]
  data.hvg <- logcounts[hvg.gene,]
  #scaledata <- Seurat::ScaleData(data.hvg,verbose=F)
  scaledata <- t(scale(t(as.matrix(data.hvg))))
  scaledata[scaledata>10] <- 10
  
  set.seed(42)
  pca <- irlba::irlba(t(scaledata),nv=50)
  cell.pca <- pca$u %*% diag(pca$d)
  rownames(cell.pca) <- colnames(scaledata)
  colnames(cell.pca) <- paste0("PC_", 1:50)
  
  knnMatrix <- RANN::nn2(cell.pca,cell.pca,k=k,searchtype = "standard")[[1]]
  snn <- getSNN(knnMatrix,1/15)
  rownames(snn) <- rownames(cell.pca)
  colnames(snn) <- rownames(cell.pca)
  
  clusters <- RunModularityClusteringCpp(
    SNN =snn ,modularityFunction = 1,resolution = 0.8,algorithm = 1,nRandomStarts = 10,
    nIterations = 10,randomSeed = 0,printOutput = F
  )
  names(clusters) <- rownames(cell.pca)
  clusters <- remove_singleton(clusters)
  clusters <- as.factor(clusters)
  return(clusters)
}

#' precise cluster
#'
#' @description clustering process after level 1
#' 
#' @param scdata single cell count data
#' @param logcounts normalized data
#' @param npcs number of princple components
#'
#' @return clusters
#' @import irlba RANN stats
#'
#' @export
precise_cluster <- function(rawdata,logdata,npcs,nhvg,k){
  hvg <- get_variableGene(rawdata)
  hvg.top <- hvg[order(hvg$variance.standardized,decreasing = T)[1:nhvg],]
  hvg.gene <- rownames(hvg.top)[hvg.top$variance>0]
  data.hvg <- logdata[hvg.gene,]
  scaledata <- t(scale(t(as.matrix(data.hvg))))
  scaledata[scaledata>10] <- 10
  
  if( ncol(scaledata)> 5000){
    set.seed(42)
    pca <- irlba::irlba(t(scaledata),nv=50)
    cell.pca <- pca$u %*% diag(pca$d)
    rownames(cell.pca) <- colnames(scaledata)
    colnames(cell.pca) <- paste0("PC_", 1:50)
  }else{
    set.seed(42)
    pca <- stats::prcomp(t(scaledata),rank. = npcs)
    cell.pca <- pca$x %*% diag(pca$sdev[1:npcs]^2)
    colnames(cell.pca) <- paste0("PC_", 1:npcs)
  }
  
  knnMatrix <- RANN::nn2(cell.pca,cell.pca,k=k,searchtype = "standard")[[1]]
  snn <- getSNN(knnMatrix,1/15)
  rownames(snn) <- rownames(cell.pca)
  colnames(snn) <- rownames(cell.pca)
  
  clusters <- RunModularityClusteringCpp(
    SNN =snn ,modularityFunction = 1,resolution = 1,algorithm = 1,nRandomStarts = 10,
    nIterations = 10,randomSeed = 0,printOutput = F
  )
  names(clusters) <- rownames(cell.pca)
  clusters <- remove_singleton(clusters)
  clusters <- as.factor(clusters)
  return(clusters)
}

#' ori_annotation
#'
#' @description  get annotation process for each cluster
#' 
#' @param marker reference cell type
#' @param clusters cluster information
#' @param i the number of level + 1
#' @param sub normalized data
#' @param p p value cut off
#'
#' @return newcluters cluster with annotaiton information
#' @import Matrix
#'
#' @export
ori_annotation <- function(marker,clusters,i,sub,p){
  celltype <- tapply(marker[,1],marker[,i],list)
  celltype <- lapply(celltype, unique)
  freq<-sort((table(unlist(celltype)))/length(celltype))
  if(max(freq)-min(freq)==0){
    weight <- rep(1,length(freq))
    names(weight) <- names(freq)
  }else{
    weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))
  }
  
  #sub <- sub[Matrix::rowSums(sub)>0,]
  meanexp <- matrix(NA,nrow = nrow(sub),ncol = length(unique(clusters)))
  for (j in 0:(ncol(meanexp)-1)) {
    groupCell <- names(clusters)[clusters==j]
    logcluster <- sub[,groupCell]
    meanexp[,j+1] <- Matrix::rowMeans(logcluster)
  }
  rownames(meanexp) <- rownames(sub)
  colnames(meanexp) <- 0:(ncol(meanexp)-1)
  
  if(i==2){
    predict_celltype<-ORA_celltype(weight,meanexp,cellType = celltype,rough = T)
  }else{
    predict_celltype<-ORA_celltype(weight,meanexp,cellType = celltype,rough = F)
  }
  newClusterIDs<-names(predict_celltype$max_cta)
  annoP <- predict_celltype$ora[unique(newClusterIDs),]
  
  annoDetails_tem <- as.data.frame(clusters)
  annoDetails_tem$clusters <- as.numeric(annoDetails_tem$clusters)
  annoDetails_tem$celltypeActivity <- predict_celltype$max_cta[annoDetails_tem$clusters]
  if(is.matrix(annoP)){
    annoDetails_tem$pValue <- apply(annoDetails_tem,1,function(x){annoP[names(predict_celltype$max_cta[x[1]]),x[1]]})
  }else{
    annoDetails_tem$pValue <- annoP[annoDetails_tem$clusters]
  }
  
  newCluster <- list()
  newCluster$annoDetails_tem <- annoDetails_tem
  
  if(nrow(meanexp) < 10000){
    newCluster$newClusterIDs <- newClusterIDs
    return(newCluster)
  }else{
    if(is.matrix(annoP)){
      for (j in 1:length(newClusterIDs)){if(annoP[newClusterIDs[j],j] >= p){newClusterIDs[j]="Unassigned"}}
    }else{
      for (j in 1:length(newClusterIDs)){if(annoP[j] >= p){newClusterIDs[j]="Unassigned"}}
    }
    newCluster$newClusterIDs <- newClusterIDs
    return(newCluster)
  }
}

#' ORA celltype
#'
#' @description detailed annotation algorithm
#' 
#' @param weight weight for each gene
#' @param meanexp mean normalized data for each cluster
#' @param cellType reference cell type
#' @param rough rough=T when level 1
#'
#' @return list
#'
#' @export
ORA_celltype<-function(weight,meanexp,cellType,rough){
  ORA_result<-matrix(NA, nrow=length(cellType),ncol=dim(meanexp)[2])
  CTA_result<-matrix(0,nrow=length(cellType),ncol=dim(meanexp)[2])
  exp_z<-scale(meanexp)
  genenames<-rownames(meanexp)
  n=0.3
  total_length <- dim(meanexp)[1]
  if(rough){
    exp_z <- meanexp
    genenames <- rownames(meanexp)
    n=0.6
  }
  for (j in 1: dim(meanexp)[2]){
    clusterexp<-meanexp[,j]
    clusterexp_z<-exp_z[,j]
    for (i in 1:length(cellType)){
      ct_exp<-length(intersect(genenames[clusterexp>0],cellType[[i]]))
      ct_not_exp<-length(cellType[[i]])-ct_exp
      exp_not_ct<-sum(clusterexp>0)-ct_exp
      #not_exp_not_ct<-length(clusterexp)-ct_not_exp
      not_exp_not_ct<-total_length-ct_not_exp
      cont.table<-matrix(c(ct_exp,ct_not_exp,exp_not_ct,not_exp_not_ct),nrow=2)
      ORA_result[i,j]<-fisher.test(cont.table,alternative="greater")$p.value
      
      weight_ss<-weight[names(weight)%in%cellType[[i]]]
      ind<-match(names(weight_ss),genenames)
      exp_ss<-clusterexp_z[ind[!is.na(ind)]]
      weight_ss<-weight_ss[!is.na(ind)]
      #print(exp_ss)
      CTA_result[i,j]<-sum(exp_ss*weight_ss)/(length(exp_ss)^n)
    }
  }
  rownames(ORA_result)<-rownames(CTA_result)<-names(cellType)
  minp_ora_ind<- apply(ORA_result,2,function(x){which.min(x)})
  minp_ora<-apply(ORA_result,2,min)
  names(minp_ora)<-rownames(ORA_result)[minp_ora_ind]
  
  max_cta_ind<- apply(CTA_result,2,function(x){which.max(x)})
  max_cta<-apply(CTA_result,2,max,na.rm=T)
  names(max_cta)<-rownames(CTA_result)[max_cta_ind]
  return(list(ora=ORA_result,cta=CTA_result,min_ora=minp_ora,max_cta=max_cta))
}

#' self defined database
#'
#' @description user-provided hierarchical cell type
#' 
#' @param file user-provided reference database
#'
#' @return celltype
#' @import tidyr
#'
#' @export
selfDefinedDatabase <- function(file,sep=","){
  file.ori <- scan(file,what = "",sep="\n")
  if(! (length(file.ori) %% 2) ==0){
    stop("Please provide correct format of self-defined database.")
  }
  
  marker.ori <- matrix("undefined",nrow=length(file.ori)/2,ncol = 2)
  for (i in 1:length(file.ori)){if(i%%2==1){marker.ori[(i+1)/2,1]=substr(file.ori[i],2,nchar(file.ori[i]))}else{marker.ori[i/2,2]=file.ori[i]}}
  ref <- as.data.frame(marker.ori[,1])
  ref <- tidyr::separate(ref,col=colnames(ref),into = c("celltype","subtypeOf"),sep = sep)
  rownames(ref) <- ref$celltype
  celltype.ori <- sapply(tapply(marker.ori[,2],marker.ori[,1],list), function(x) strsplit(x,","))
  
  celltype.tem <- as.data.frame(matrix(unlist(celltype.ori)))
  celltype.tem[,2] <- rep(names(celltype.ori),lengths(celltype.ori))
  celltype.tem <- tidyr::separate(celltype.tem,col=V2,into = c("celltype","subtypeOf"),sep = ",")
  
  ct <- as.data.frame(matrix("Undefined",nrow = nrow(celltype.tem),ncol = length(unique(celltype.tem$subtypeOf))+2),stringsAsFactors=FALSE)
  ct[,1:2] <- celltype.tem[,1:2]
  for (i in 3:ncol(ct)){
    for (j in 1:nrow(ct)) {if(ct[j,i-1] %in% rownames(ref)){ct[j,i]<-ref[ct[j,i-1],2]}else{ct[j,i]<-ct[j,i-1]}}
  }
  
  layer <- ncol(ct)
  if(! identical(ct[,ncol(ct)],ct[,ncol(ct)-1])){layer <- ncol(ct)}else{
    for (i in 1:ncol(ct)) {if(identical(ct[,i],ct[,i+1])){layer <- i;break}}
  }
  
  celltype.trim <- ct[,1:layer]
  for (i in 1:nrow(celltype.trim)){
    tag <- which(celltype.trim[i,]==celltype.trim[i,ncol(celltype.trim)])
    if(length(tag)>1){
      if(ncol(celltype.trim)+1-length(tag)==3){suf <- ct[i,3]}
      else{suf <- ct[i,3:(ncol(celltype.trim)+1-length(tag))]}
      pre <- rep(celltype.trim[i,2],length(tag)-1)
      celltype.trim[i,] <- c(celltype.trim[i,][1:2],pre,suf)
    }
  }
  
  for (i in 1:ncol(celltype.trim)){if(i == 1){colnames(celltype.trim)[i]="gene"}else{colnames(celltype.trim)[i]=paste0("Level",(ncol(celltype.trim)+1-i),"")}}
  
  ct.final <- celltype.trim
  for (i in 1:ncol(ct.final)){
    if(i >1){ct.final[,i]=celltype.trim[,ncol(celltype.trim)+2-i];colnames(ct.final)[i] <- colnames(celltype.trim)[ncol(celltype.trim)+2-i]}
  }
  #ct.final[,1] <- toupper(ct.final[,1])
  #ct.final <- ct.final[,apply(ct.final, 2, function(x) {length(unique(x))})>1]
  ct.final[] <- lapply(ct.final, as.character)
  
  return(ct.final)
}

#' generate a html file to visualize hierarchical cell type structure
#' 
#' @description generate a html file to visualize hierarchical cell type structure
#' 
#' @param celltype reference cell type for visualization
#'
#' @return html file
#' @import data.tree networkD3
#'
#' @export
databaseVisual <- function(celltype,repeated=F){
  celltype.tree <- celltype[,2:ncol(celltype)]
  celltype.tree <- celltype.tree[complete.cases(celltype.tree),]
  
  celltype.path <- as.character()
  for (i in 1:nrow(celltype.tree)) {
    if(repeated){
      celltype.path[i] <- paste(as.character(celltype.tree[i,]),collapse = "/")
    }else{
      celltype.path[i] <- paste(unique(as.character(celltype.tree[i,])),collapse = "/")
    }
    if(length(unique(celltype.tree[,1])) > 1){
      celltype.path[i] <- paste("Database",celltype.path[i],sep = "/")
    }
  }
  celltype.tree$pathString <- celltype.path
  population <- data.tree::as.Node(celltype.tree)
  
  useRtreeList <- data.tree::ToListExplicit(population, unname = TRUE)
  Fsize <- ifelse(length(unique(celltype.tree[,ncol(celltype.tree)]))>=50,10,15)
  networkD3::saveNetwork(networkD3::radialNetwork(useRtreeList,nodeStroke = "orange",
                                                  fontSize = Fsize),file = "./Database.html")
}

#' add genes
#'
#' @description add genes to existed cell type
#' 
#' @param geneCellTypeList user-provided genes and cell types added to existed reference
#' @param celltype existed reference cell types
#'
#' @return cell type with added genes
#'
#' @export
addGene <- function(geneCellTypeList,celltype){
  if(!is(class(geneCellTypeList),"data.frame")){
    geneCellTypeList <- as.data.frame(geneCellTypeList)
  }
  
  if(ncol(geneCellTypeList) !=2 ){
    #cat("Please provide correct format of gene to celltype list.\n")
    stop("Please provide correct format of gene to celltype list.\n")
  }else{
    ct <- unique(geneCellTypeList[,2])
    ct_db <- unique(as.character(celltype[,ncol(celltype)]))
    for (i in ct) {
      if(i %in% ct_db){cat(i,"in the databse.\n")}else{cat(i,"not in the databse.\n")}
    }
    ct_list <- geneCellTypeList[geneCellTypeList[,2] %in% ct_db,]
    if(ncol(ct_list)==0){stop("No provided celltype in the databse.\n")}
    db_list <- unique(celltype[,2:ncol(celltype)])
    db_list <- db_list[complete.cases(db_list),]
    rownames(db_list) <- db_list[,ncol(db_list)]
    for (i in 1:nrow(ct_list)) {
      celltype <- rbind(celltype,c(ct_list[i,1],as.character(db_list[ct_list[i,2],])))
    }
    return(celltype)
  }
}
