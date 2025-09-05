library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
#loading E15 mouse and spatial data
cds <- readRDS("datacds.rds")
cds1 <- readRDS("/Users/naghamkhourifarah/Downloads/GSE245469_RAW/data/combined_HQ.rds")
cds2 <- subset(cds,cells = rownames(cds@meta.data[cds$stage == "E15",]))
Idents(cds2) <- cds2$ectmes

anchors <- FindTransferAnchors(reference = cds2, query = cds1,dims = 1:30,  features =intersect(rownames(cds2),rownames(cds1)),reduction = "cca")
predicted.labels <- TransferData(anchorset = anchors, refdata = cds2$ectmes , dims = 1:30,weight.reduction =cds1[['harmony']] )
cds1$ectmes <- predicted.labels$predicted.id

SpatialPlot(cds1,ncol = 2,images = c("E14_1","E15_1"),pt.size.factor = 10,group.by = "ectmes")
Idents(cds1) <- "ectmes"
dir.create("figures/inspection_ectmes")
for (n in levels(cds1)){
  ids = WhichCells(cds1, idents = n)
  p1 = SpatialPlot(cds1, cells.highlight = ids,images = c("E15_1"), pt.size.factor = 10, label = F)+NoAxes()+NoLegend()
  png(paste0("figures/inspection_ectmes/cluster_",n,".png"), w=3000, h=10000)
  plot(p1)
  dev.off()
}
