setwd("~/Projects/Frankfurt/StellateCellsACLF/AnalysisSpatTrans/")

library(Seurat) 
library(ggplot2)
#library(SeuratData)
#library(ggplot2)
#library(patchwork)
library(dplyr)
#for faster computation of Wilcoxon Rank Sum tests
#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')
#library("presto")



seu.obj <- readRDS("liver1_seurat_object.Rds")
seu.obj2 <- readRDS("liver2_seurat_object.Rds")

#test frmom vignette
head(seu.obj@meta.data)

#names of all measured gene markers:
Foxrow=which(rownames(seu.obj@assays$Nanostring@counts) == "Foxo1")
Actarow=which(rownames(seu.obj2@assays$Nanostring@counts) == "Acta2")
Col1a1row=which(rownames(seu.obj@assays$Nanostring@counts) == "Col1a1")
Pdgfrbrow=which(rownames(seu.obj@assays$Nanostring@counts) == "Pdgfrb")
#ausserdem Perycite marker PECAM1, RGS5, DES, VTN, (according to PanglaoDB)

#create expression data frame
  

#Create expression dataframe for plotting
#first tissue 1
data=as.data.frame(t(as.matrix(seu.obj@assays$Nanostring@counts)))
data$cellType=seu.obj@meta.data$nb_clus
data$fov=as.numeric(seu.obj@meta.data$fov)
data$leiden_clus=seu.obj@meta.data$leiden_clus
data$tissue=seu.obj@meta.data$tissue

data$condition=data$fov
#add UMAP coordinate of cells
data$umap1=seu.obj@reductions$umap@cell.embeddings[,1]
data$umap2=seu.obj@reductions$umap@cell.embeddings[,2]
#annotate fovs
healthy=c(1,2,15,16,17,18,19,20,21,22,23,45)
cirr=c(5,6,24,25,26,30,31,34,35,40,41,44)
steat=c(3,4,9,10,13,14,27,28,29,38,39)
aclf=c(7,8,11,12,32,33,36,37,42,43)
data$condition[data$condition %in% healthy] = "healthy"
data$condition[data$condition %in% cirr] = "cirrhosis"
data$condition[data$condition %in% steat] = "steatosis"
data$condition[data$condition %in% aclf] = "aclf"


#now tissue 2
data2=as.data.frame(t(as.matrix(seu.obj2@assays$Nanostring@counts)))
data2$cellType=seu.obj2@meta.data$nb_clus
data2$fov=as.numeric(seu.obj2@meta.data$fov)
data2$leiden_clus=seu.obj2@meta.data$leiden_clus
data2$tissue=seu.obj2@meta.data$tissue

data2$condition=data2$fov
#add UMAP coordinate of cells
data2$umap1=seu.obj2@reductions$umap@cell.embeddings[,1]
data2$umap2=seu.obj2@reductions$umap@cell.embeddings[,2]
#annotate fovs
cirr=c(1,2,3,4,5,6,7,8,13,14,17,18,19,20,21,22,23,24,25,36,37,44,45)
aclf=c(9,10,11,12,15,16,26,27,28,29,30,31,32,33,34,35,38,39,40,41,42,43)
data2$condition[data2$condition %in% cirr] = "cirrhosis"
data2$condition[data2$condition %in% aclf] = "aclf"
#merge both dataframes
data=as.data.frame(rbind(data,data2))

data$cellType=as.factor(data$cellType)



#plot per cell cluster
#increase default font size
theme_set(theme_gray(base_size = 18))

p <- ggplot(data, aes(x=cellType, y=Foxo1)) + 
  geom_violin()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p <- ggplot(data, aes(x=cellType, y=Acta2)) + 
  geom_violin()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,7)

p <- ggplot(data, aes(x=cellType, y=Rgs5)) + 
  geom_violin()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,7)




p <- ggplot(data, aes(x=cellType, y=Col1A1)) + 
  geom_violin()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,7)

p <- ggplot(data, aes(x=cellType, y=Pdgfrb)) + 
  geom_violin()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,7)

#plot per condition
p <- ggplot(data, aes(x=condition, y=Acta2)) + 
  geom_violin() + ylim(0,10)
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#plot the clustering with condition labels
ggplot(data,aes(x=umap1,y=umap2,col=condition))+
  geom_point()
ggsave("UMAP-Condition.pdf")

ggplot(data,aes(x=umap1,y=umap2,col=cellType,shape=condition))+
  geom_point()
ggsave("UMAP_cluster_cellType.pdf")

ggplot(subset(data,condition=="aclf"),aes(x=umap1,y=umap2,shape=cellType,col=cellType))+
  geom_point() #+ scale_color_gradientn(colours = rainbow(4))
ggsave("UMAP_aclf.pdf")
ggplot(subset(data,condition=="healthy"),aes(x=umap1,y=umap2,shape=cellType,col=cellType))+
  geom_point() #+ scale_color_gradientn(colours = rainbow(4))
ggsave("UMAP_healthy.pdf")
ggplot(subset(data,condition=="cirrhosis"),aes(x=umap1,y=umap2,shape=cellType,col=cellType))+
  geom_point() #+ scale_color_gradientn(colours = rainbow(4))
ggsave("UMAP_cirrhosis.pdf")
ggplot(subset(data,condition=="steatosis"),aes(x=umap1,y=umap2,shape=cellType,col=cellType))+
  geom_point() #+ scale_color_gradientn(colours = rainbow(4))
ggsave("UMAP_steatodsis.pdf")

ggplot(data,aes(x=umap1,y=umap2,shape=condition,col=Casp3))+
  geom_point() + scale_color_gradientn(colours = rainbow(4))
  

###compute conditions per cluster
CondPerCluster <- function(){
  res=c()
  for (i in unique(data$cellType)){
    print(i)
 
    sub=subset(data, cellType == i)
    tab=table(sub$condition)
    for (j in 1:length(tab)){
      print(paste(c(i,names(tab)[j],tab[j]),sep="-"))
      res=rbind(res,c(i,names(tab)[j],tab[j]))
    }
  }
  return(res)
} 

result=as.data.frame(CondPerCluster())
names(result)=c("cellType","Condition","CellCount")
result$CellCount=as.numeric(result$CellCount)


ggplot(result,aes(fill=Condition,x=cellType,y=CellCount))+
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#percentage of cell type per Condition
ggplot(result,aes(fill=cellType,x=Condition,y=CellCount))+
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#cells per condition
ggplot(result,aes(x=Condition,y=CellCount))+
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("CellsPerCondition.pdf")


#estimate average expression and per cell expression per category
#######################
AveragePerCategory <- function(){
  res=c()
  for (i in unique(data$cellType)){
    for (j in unique(data$condition)){
      sub=subset(data, cellType == i & condition == j)
      print(paste(c("LookingAt:",i,j),sep="-"))
      res=rbind(res,c(as.numeric(apply(sub[,1:970],2,mean)),i,j))
    }
  }
  res=as.data.frame(res)
  names(res)[1:970]=names(data)[1:970]
  names(res)[971:972]=c("cellType","condition")
  return(res)
} 

#compute aveargae for all genes per category
avg=AveragePerCategory()
avg$Rgs5=as.numeric(avg$Rgs5)
avg$Pdgfrb=as.numeric(avg$Pdgfrb)
avg$Acta2=as.numeric(avg$Acta2)
avg$Pecam1=as.numeric(avg$Pecam1)
avg$Csf1r=as.numeric(avg$Csf1r)
#Apoptosis Add1, Apc , Atm, Bcl2, BDnf, Casp3 , Cd14, 
avg$Casp3=as.numeric(avg$Casp3)
avg$Apc=as.numeric(avg$Apc)

p <- ggplot(avg, aes(x=cellType, y=condition,color=Rgs5,size=Rgs5)) + 
  geom_point()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("RGS5ExpressionPerCategory.pdf")
p <- ggplot(avg, aes(x=cellType, y=condition,color=Pdgfrb,size=Pdgfrb)) + 
  geom_point()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("PdgfrbExpressionPerCategory.pdf")
p <- ggplot(avg, aes(x=cellType, y=condition,color=Acta2,size=Acta2)) + 
  geom_point()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Acta2ExpressionPerCategory.pdf")

p <- ggplot(avg, aes(x=cellType, y=condition,color=Pecam1,size=Pecam1)) + 
  geom_point()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Pecam1ExpressionPerCategory.pdf")

p <- ggplot(avg, aes(x=cellType, y=condition,color=Csf1r,size=Csf1r)) + 
  geom_point()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Csf1rExpressionPerCategory.pdf")
p <- ggplot(avg, aes(x=cellType, y=condition,color=Casp3,size=Casp3)) + 
  geom_point()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Casp3ExpressionPerCategory.pdf")
p <- ggplot(avg, aes(x=cellType, y=condition,color=Apc,size=Apc)) + 
  geom_point()
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Casp3ExpressionPerCategory.pdf")


#estimate overall counts of all genes
mean1=apply(as.matrix(seu.obj@assays$Nanostring@counts),1,mean)
mean2=apply(as.matrix(seu.obj2@assays$Nanostring@counts),1,mean)
meanAll=mean1+mean2/2
meanAll=sort(meanAll,decreasing=T)
TotalCounts=data.frame(names=names(meanAll),meanAll)
write.table(TotalCounts,"MeanALlGenes.csv",sep=";",quote=F,col.names=T,row.names=F)
# Image names. Each slide is stored as a separate image within the object.
Images(seu.obj)
# Get name of the first image
image1 <- Images(seu.obj)[1]



# Plot all cells.
# We recommend setting the border color to 'NA' as the default 'white' often masks all cells when zoomed out, leading to a fully white plot.
ImageDimPlot(seu.obj, fov = image1, axes = TRUE, border.color = NA)

use_fov=20
use_slide_image <- Images(seu.obj)[1] # Slide desired, as named in images 

use_slide_metadata <- seu.obj@meta.data$Run_Tissue_name[1] # Slide desired, as named in the metadata column ‘Run_Tissue_name’ 

# First get cells in your FOV 

cells_of_interest <- seu.obj$id[(seu.obj$fov == use_fov) & ( seu.obj$Run_Tissue_name == use_slide_metadata)] 

# Then find spatial boundaries of the rectangle containing the centroids of these cells 

centroid_data <- seu.obj@images[[use_slide_image]]$centroids 

zoom_fov <- apply(centroid_data@coords[centroid_data@cells %in% cells_of_interest,], 2, range)

ImageDimPlot(seu.obj,
            fov = use_slide_image,
            border.color = "black") + xlim(c(-20000,-10000)) + ylim(c(-170000,-145000))

ImageFeaturePlot(seu.obj,fov = use_slide_image, border.color = NA,features = "Rgs5")+ xlim(c(-20000,-10000)) + ylim(c(-170000,-155000))

ImageDimPlot(seu.obj,
             fov = use_slide_image,
             border.color = "black",
             alpha = 0.5, # Reduce alpha of cell fills to better visualize the overlaying molcules
             molecules = c("Rgs5", "Pdgfrb", "Acta2","Pdgfra"),
             mols.size = 0.8,
             nmols = 100000, # Set the total number of molecules to visualize, noting that this is across the whole slide
             axes = FALSE)+ xlim(c(-20000,-10000)) + ylim(c(-168000,-158000))

#
ggsave("ExampleImage.pdf")

FeaturePlot(seu.obj, 
            features = "Rgs5",
            order = TRUE)

seu.obj.markers <- FindAllMarkers(seu.obj, only.pos = TRUE)
seu.obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

write.csv2(seu.obj.markers,"ClusterMarkerGenes.csv",sep = ";",quote=F,row.names=F)


