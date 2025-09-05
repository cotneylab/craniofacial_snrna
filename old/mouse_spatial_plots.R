library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
#E15.5 spatial plots
#loading mouse facial transcriptomics and E15 spatial data from Pina et al 
cds <- readRDS("face_mouse_final.rds")
cds1 <- readRDS("data/combined_HQ.rds")
cds2 <- subset(cds,cells = rownames(cds@meta.data[cds$stage == "E15",]))
Idents(cds2) <- cds2$annotation
#project the scRNA from our data onto the E15.5 head section spatial transcriptomics 
anchors <- FindTransferAnchors(reference = cds, query = cds1,dims = 1:30,  features =intersect(VariableFeatures(cds),rownames(cds1)),reduction = "cca")
predicted.labels <- TransferData(anchorset = anchors, refdata = cds$annotation , weight.reduction = cds2[['pca']],dims = 1:30 )
cds1$annotation <- predicted.labels$predicted.id
predicted.labels <- TransferData(anchorset = anchors, refdata = cds$cell_type1 , dims = 1:30 )
cds1$cell_type <- predicted.labels$predicted.id
cds1$cell_type <- factor(cds1$cell_type,levels = levels(cds$cell_type1))
#Plot the transfered annotations
Idents(cds1) <- "cell_type"
SpatialPlot(cds1,images = c("E15_1"),pt.size.factor = 10)
Idents(cds1) <- "annotation"
SpatialPlot(cds1,images = c("E15_1"),pt.size.factor = 10)

dir.create("figures/inspection")
for (n in levels(cds1)){
  ids = WhichCells(cds1, idents = n)
  p1 = SpatialPlot(cds1, cells.highlight = ids,images = c("E15_1"), pt.size.factor = 10, label = F)+NoAxes()+NoLegend()
  png(paste0("figures/inspection/cluster_",n,".png"), w=3000, h=3000)
  plot(p1)
  dev.off()
}


#E11.5 Spatial plots
#convert Mosta files to seurat readable format
Convert("E11.5_E1S1.MOSTA.h5ad", dest = "h5seurat", overwrite = FALSE)
Convert("E11.5_E1S3.MOSTA.h5ad", dest = "h5seurat", overwrite = FALSE)
Convert("E11.5_E1S4.MOSTA.h5ad", dest = "h5seurat", overwrite = FALSE)
#load spatial data
seu1 <-  LoadH5Seurat("/Users/naghamkhourifarah/Downloads/E11.5_E1S1.MOSTA.h5seurat",assays ="RNA")
p <- DimPlot(seu1, reduction = "spatial", group.by = "annotation")+NoAxes()+NoLegend()
Embedding <- as.data.frame(seu1@reductions$spatial@cell.embeddings)
seu1$x <- as.numeric(Embedding$spatial_1)+150
seu1$y <- -as.numeric(Embedding$spatial_2)+400
seu1$section <- "S01"

seu2 <-  LoadH5Seurat("/Users/naghamkhourifarah/Downloads/E11.5_E1S3.MOSTA.h5seurat",assays ="RNA")
p <- DimPlot(seu2, reduction = "spatial", group.by = "annotation")+NoAxes()+NoLegend()
Embedding <- as.data.frame(seu2@reductions$spatial@cell.embeddings)
seu2$x <- as.numeric(Embedding$spatial_1)-150
seu2$y <- -as.numeric(Embedding$spatial_2)-50
seu2$section <- "S02"

seu3 <-  LoadH5Seurat("/Users/naghamkhourifarah/Downloads/E11.5_E1S4.MOSTA.h5seurat",assays ="RNA")
p <- DimPlot(seu3, reduction = "spatial", group.by = "annotation")+NoAxes()+NoLegend()
Embedding <- as.data.frame(seu3@reductions$spatial@cell.embeddings)
seu3$x <- as.numeric(Embedding$spatial_1)+100
seu3$y <- -as.numeric(Embedding$spatial_2)+300
seu3$section <- "S03"



obj.list <- list(seu1,seu2,seu3)
#create an image for each section
for (j in 1:3) {
  
positions = obj.list[[j]]@meta.data[, c("x","y")]
obj.list[[j]][["image"]] <- new(
  Class = 'SlideSeq',
  assay = 'RNA',
  coordinates = positions
)
}
seu <- cds

for (i in levels(seu)) {
  sig <- FindMarkers(seu,ident.1 = i,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                     min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                     only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                     latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                     pseudocount.use = 1, return.thresh = 0.01)
  #top  <- sig %>% top_n(-100, pct.2);top
  top  <- sig %>% top_n(100, avg_log2FC);top
    obj.list[[1]] <- AddModuleScore(object = obj.list[[1]],features = list(rownames(top)),name = c(paste0(i,"_")))
  

}

library(patchwork)
p <- list()



for (i in levels(seu)) {
  p[[i]]<- FeaturePlot(obj.list[[1]],paste0(i,"_1"), reduction = "spatial",pt.size = 0.2)+scale_color_viridis_b()+NoAxes()+NoLegend()+ scale_y_reverse() 
  
}

pdf("E11_all_top100_section1.pdf",w=12,h=24)
p
dev.off()
for (i in levels(seu)) {
  p[[i]]<- FeaturePlot(obj.list[[1]],paste0(i,"_1"), reduction = "spatial",cols = c("grey90","navy"),pt.size = 0.2)+NoAxes()+NoLegend()+ scale_y_reverse() 
  
}

for (i in levels(seu)) {
  p[[i]]<- SpatialFeaturePlot(obj.list[[1]],paste0(i,"_1"),pt.size.factor = 1.4)+NoAxes()+ scale_y_reverse() 
  
}
pdf("E11_allModuleScore_top100_pct20_SpatialFeaturePlot_sec3.pdf",w=5,h=5)
p
dev.off()
pdf("E11_mcncc_top100_singleSection_SpatialFeaturePlot.pdf",w=12,h=24)
wrap_plots(p,ncol = 2)
dev.off()

p <- list()

figFeat <- c("Nr2f2","Pax3","Tfap2a","Tfap2b","Alx4","Crabp1","Hand2","Nkx2-1")
for (i in figFeat) {
  p[[i]] <- SpatialFeaturePlot(obj.list[[1]],paste0(i),pt.size.factor = 1.4)+NoAxes()+ scale_y_reverse()
  }
pdf("E11_figure4_markers.pdf",w=12,h=24)
wrap_plots(p,ncol = 2)
dev.off()
p <- list()

for (i in figFeat) {
  p[[i]] <- FeaturePlot(obj.list[[1]],paste0(i), reduction = "spatial",pt.size = 0.2,)+NoAxes()+NoLegend()+ scale_y_reverse()+scale_color_viridis_b()
}
pdf("E11_figure4_markers_image.pdf",w=12,h=24)
wrap_plots(p,ncol = 2)
dev.off()

#calculate module scores for all subtypes
for (i in levels(seu)) {
  for (j in 1:3) {
    sig <- FindMarkers(seu,ident.1 = i,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                       min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                       only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                       latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                       pseudocount.use = 1, return.thresh = 0.01)
    #top  <- sig %>% top_n(-100, pct.2);top
    top  <- sig %>% top_n(100, avg_log2FC);top
    obj.list[[j]] <- AddModuleScore(object = obj.list[[j]],features = list(rownames(top)),name = c(paste0(i,"_")))
    
  }
}

p <- list()

for (i in levels(seu)) {
  p1 <- FeaturePlot(obj.list[[1]],paste0(i,"_1"), reduction = "spatial",pt.size = 0.2)+scale_color_viridis_c()+NoAxes()+NoLegend()+ scale_y_reverse() + ggtitle(NULL)
  p2 <- FeaturePlot(obj.list[[2]],paste0(i,"_1"), reduction = "spatial", pt.size = 0.2)+scale_color_viridis_c()+NoAxes()+NoLegend()+ scale_y_reverse() + ggtitle(NULL)
  p3 <- FeaturePlot(obj.list[[3]],paste0(i,"_1"), reduction = "spatial", pt.size = 0.2)+scale_color_viridis_c()+NoAxes()+NoLegend()+ scale_y_reverse() + ggtitle(NULL)
  p[[i]]<- wrap_plots(p1,p2,p3,ncol = 3)+ plot_annotation(i) & theme(plot.tag = element_text(size = 8,face = "bold"))
}


pdf("E11_all_top100_3sections.pdf",w=18,h=9)
p
dev.off()




