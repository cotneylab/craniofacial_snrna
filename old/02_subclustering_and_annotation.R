library(Seurat)
library(tidyverse)
library(harmony)
library(SeuratWrappers)
library(patchwork)
#load full processed human data
cds2<- readRDS("data/cds2.rds")
#create an ensembl gene annotation table linking mouse and human genes
require(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
annot_table <- getLDS(
  mart = human,
  attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name','chromosome_name','gene_biotype','description'),
  martL = mouse,
  attributesL = c('ensembl_gene_id','mgi_symbol','external_gene_name','chromosome_name'))

#subset and cluster ectoderm from human full dataset
ect <- subset(cds2, idents = "ectoderm")
ect <-   NormalizeData(ect) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE) %>%
  RunHarmony(group.by.vars = "orig.ident")
ElbowPlot(ect,ndims = 50)
ect <- RunUMAP(ect,reduction = "harmony", dims = 1:30,min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1,2,3), verbose = FALSE)
#plot cluster tree at various resolutions 
library(clustree)
p4 <- clustree(ect, prefix = "RNA_snn_res.")
p4
Idents(ect) <- "RNA_snn_res.0.8"
p1 <- DimPlot(ect, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(ect, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(ect,group.by = "stage")+NoAxes()
pdf("figures/ect_sub_all_final.pdf",width = 16,height = 8)
p1+p2
p4
dev.off()

#load mouse ectoderm subtypes
mect <- readRDS("/Volumes/Extreme SSD/Mouse_cellranger2/data/ect_all.rds")
#generate a table of mouse/human orthologs of the features in the object  
orth <- as.data.frame(unique(annot_table$HGNC.symbol[match(rownames(mect),annot_table$MGI.symbol)]))
colnames(orth) <- "hGenes"
orth$mGenes <- annot_table$MGI.symbol[match(orth$hGenes,annot_table$HGNC.symbol)]
orth1 <- orth[!duplicated(orth$mGenes), ]
common <- intersect(orth1$mGenes,rownames(mect))
length(common)
#  15907
#recreate an orthologous seurat object 
metaData <- mect@meta.data
count <- GetAssayData(mect, assay = "RNA",slot = "counts")
count <- count[common,]
rownames(count) <- orth1$hGenes[match(rownames(count),orth1$mGenes)]
mect1 <- CreateSeuratObject(counts = count,meta.data = metaData)
#process and re-annotate the resulting dataset
mect1 <- NormalizeData(mect1) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE)
mect1 <- RunHarmony(mect1,group.by.vars = "orig.ident")
ElbowPlot(mect1,ndims = 50)
mect1 <- RunUMAP(mect1,reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors( dims = 1:30, reduction = "harmony",verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1,1.5,2,3), verbose = FALSE)
#refine clusters annotation based on original metadata 
Idents(mect1) <- "RNA_snn_res.1"
mect1 <- RenameIdents(mect1, c("0"="surface2","1"="ect.EBF","2"="palate1","3"="lOE","4"="palate.surface","5"="dental2","6"="surface3","7"="dental1","8"="palate2",
                               "9"="surface1","10"="periderm","11"="NaP","12"="surface4","13"="eOE.PAX7","14"="cilia.OE",
                               "15"="iOE","16"="fusion zone","17"="auditory","18"="pituitary","19"="OE3","20"="myeloid","21"="ect.EBF"))

#transfer the annotations fromm mouse (with human features) to human ectoderm subset
anchors <- FindTransferAnchors(reference = mect1, query = ect, reduction = 'cca', normalization.method = "LogNormalize", dims = 1:30,features = rownames(mect1))
predicted.id <- TransferData(anchorset = anchors, refdata = Idents(mect1), weight.reduction = ect[['harmony']],dims = 1:30)
ect$predicted <- predicted.id$predicted.id
p1 <- DimPlot(ect, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(ect,group.by = "predicted",label = T)+NoAxes()+NoLegend()
Idents(ect) <- "RNA_snn_res.0.8"
#Plot unique features and get specific markers (and ontology) for more accurate annotations
FeaturePlot(ect, "SHH",label = T)

#find marker genes for ectoderm subtypes
sig_genes <- FindAllMarkers(ect, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
top <- sig_genes %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
p <- DotPlot(object = ect, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
df <- data.frame(Gene = p$data$features.plot, avg.exp = p$data$avg.exp.scaled, pct.exp = p$data$pct.exp, cluster = p$data$id)
p2 <- df %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p2

#annotate ectoderm sub-clusters
ect <- RenameIdents(ect, c("0"="dental","1"="surface2","2"="ect.GDNF","3"="palate","4"="auditory1", "5"="iOE1","6"="ect.EBF",
                           "7"="auditory2", "8"="NaP.PT", "9"="palate.surface","10"="fusion.zone","11"="surface3",
                           "12"="periderm","13"="pituitary","14"="iOE2","15"="myeloid","16"="OE3","17"="auditory3",
                           "18"="eOE","19"="thyroid","20"="e.dental", "21"="eye.ect"))
ect$anno_res0.8 <- Idents(ect)
#combine similar annotations
ect <- RenameIdents(ect, c("auditory1"="auditory","auditory2"="auditory","auditory3"="auditory","myeloid"="ect.EBF",
                           "surface2"="surface","surface3"="surface","eOE"="surface", "palate.surface"="surface","e.dental"="dental"))
ect$anno_reduced <- Idents(ect)
p3 <- DimPlot(ect,group.by = "anno_res0.8",label = T)+NoAxes()+NoLegend()
p4 <- DimPlot(ect,group.by = "anno_reduced",label = T)+NoAxes()+NoLegend()
pdf("figures/ect_sub_all_annotated_final.pdf",width = 16,height = 8)
p1+p2
p3+p4
dev.off()

saveRDS(anchors,"data/ect_allMtoH_anchors.rds")
saveRDS(mect,"data/mect_final.rds")
saveRDS(mect1,"data/mect1_final.rds")
saveRDS(ect,"data/ect_final.rds")


#subset and cluster mesenchyme from human full dataset
mes <- subset(cds2, idents = "mesenchyme")
Idents(mes) <- "stage"
mes <- subset(mes, idents=c("CS20"))
mes <- CellCycleScoring(mes, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
mes <-   NormalizeData(mes) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE)
mes <- harmony::RunHarmony(mes,group.by.vars = "orig.ident")
ElbowPlot(mes,ndims = 50)
mes <- RunUMAP(mes,reduction = "harmony", dims = 1:30,min.dist = 0.3,) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1,2,3), verbose = FALSE)
#plot cluster tree at various resolutions 
library(clustree)
p4 <- clustree(mes, prefix = "RNA_snn_res.")
p4
Idents(mes) <- "RNA_snn_res.0.8"
p1 <- DimPlot(mes, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(mes, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(mes,split.by = "stage",label = T)+NoAxes()+NoLegend()
pdf("figures/mes_sub_all_final.pdf",width = 16,height = 8)
p1+p2
p4
dev.off()

#load mouse mesenchyme subtypes
mmes <- readRDS("/Volumes/Extreme SSD/Mouse_cellranger2/data/mes_all.rds")
#generate a table of mouse/human orthologs of the features in the object  
orth <- as.data.frame(unique(annot_table$HGNC.symbol[match(rownames(mmes),annot_table$MGI.symbol)]))
colnames(orth) <- "hGenes"
orth$mGenes <- annot_table$MGI.symbol[match(orth$hGenes,annot_table$HGNC.symbol)]
orth1 <- orth[!duplicated(orth$mGenes), ]
common <- intersect(orth1$mGenes,rownames(mmes))
length(common)
#  15907
#recreate an orthologous seurat object 
metaData <- mmes@meta.data
count <- GetAssayData(mmes, assay = "RNA",slot = "counts")
count <- count[common,]
rownames(count) <- orth1$hGenes[match(rownames(count),orth1$mGenes)]
mmes1 <- CreateSeuratObject(counts = count,meta.data = metaData)
#process and re-annotate the resulting dataset
mmes1 <- NormalizeData(mmes1) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE)
mmes1 <- RunHarmony(mmes1,group.by.vars = "orig.ident")
ElbowPlot(mmes1,ndims = 50)
mmes1 <- RunUMAP(mmes1,reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors( dims = 1:30, reduction = "harmony",verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1,1.5,2,3), verbose = FALSE)
DimPlot(mmes1, label = T)+NoAxes()+NoLegend()
p4 <- clustree(mmes1, prefix = "RNA_snn_res.")
p4
Idents(mmes1) <- "RNA_snn_res.0.8"

#transfer the annotations from mouse mesenchyme (with human features) to human mesenchyme subset
anchors <- FindTransferAnchors(reference = mmes1, query = mes, reduction = 'cca', normalization.method = "LogNormalize", dims = 1:30,features = VariableFeatures(mmes1))
predicted.id <- TransferData(anchorset = anchors, refdata = mmes1$cellType4, weight.reduction = mes[['harmony']],dims = 1:30)
mes$predicted <- predicted.id$predicted.id
p1 <- DimPlot(mes, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(mes,group.by = "predicted",label = T)+NoAxes()+NoLegend()
Idents(mes) <- "RNA_snn_res.0.8"
#Plot unique features and get specific markers (and ontology) for more accurate annotations
FeaturePlot(mes, "TFAP2A",label = T)

#find marker genes for mesenchymal subtypes
sig_genes <- FindAllMarkers(mes, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
#remove "LIN" genes from markers lists
LINFeatures <- grep("LIN",rownames(mes))
mes1 <- rownames(mes)[-LINFeatures]
sig_genes2 <- filter(sig_genes, gene %in% mes1)
top <- sig_genes2 %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
p <- DotPlot(object = mes, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
df <- data.frame(Gene = p$data$features.plot, avg.exp = p$data$avg.exp.scaled, pct.exp = p$data$pct.exp, cluster = p$data$id)
p2 <- df %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p2

#annotate mesenchyme sub-clusters
mes <- RenameIdents(mes, c("0"="mandibular arch1","1"="MxP2","2"="palatal shelf2.1","3"="LNP1","4"="mandibular arch2", "5"="palatal shelf2.2","6"="LNP2",
                           "7"="mandibular arch3", "8"="palatal.fusion1", "9"="pLNP2","10"="fusion.mes1","11"="palatal.fusion2",
                           "12"="cartilage1","13"="MxP.aLNP","14"="palatalShelf1.1","15"="pLNP.fusion","16"="palatalShelf1.2","17"="MxP.surface",
                           "18"="mandibular arch3","19"="l.osteoblat","20"="fusion.mes2", "21"="l.osteoblast","22"="cartilage2"))
mes$anno_res0.8 <- Idents(mes)

#annotate mesenchyme sub-clusters
mes <- RenameIdents(mes, c("0"="cartilage1","1"="e.osteoblat","2"="palatal shelf2.1","3"="LNP1","4"="mandibular arch2", "5"="palatal shelf2.2","6"="LNP2",
                           "7"="mandibular arch3", "8"="palatal.fusion1", "9"="pLNP2","10"="fusion.mes1","11"="palatal.fusion2",
                           "12"="cartilage1","13"="MxP.aLNP","14"="palatalShelf1.1","15"="pLNP.fusion","16"="palatalShelf1.2","17"="MxP.CALB2",
                           "18"="mandibular arch3","19"="e.osteoblast","20"="fusion.mes2", "21"="l.osteoblast","23"="cartilage2"))
mes$anno_res0.8 <- Idents(mes)
#combine similar annotations
p3 <- DimPlot(mes,group.by = "anno_res0.8",label = T)+NoAxes()+NoLegend()
p4 <- DimPlot(mes,group.by = "anno_reduced",label = T)+NoAxes()+NoLegend()
pdf("figures/mes_sub_all_annotated_final.pdf",width = 16,height = 8)
p1+p2
p3+p4
dev.off()

saveRDS(anchors,"data/mes_allMtoH_anchors.rds")
saveRDS(mmes,"data/mmes_final.rds")
saveRDS(mmes1,"data/mmes1_final.rds")
saveRDS(mes,"data/mes_final.rds")



#subset and cluster ectoderm from human full dataset
cncc <- subset(cds2, idents = "CNCC")
cncc <-   NormalizeData(cncc) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE) %>%
  RunHarmony(group.by.vars = "orig.ident")
ElbowPlot(cncc,ndims = 50)
cncc <- RunUMAP(cncc,reduction = "harmony", dims = 1:30,min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1,2,3), verbose = FALSE)
#plot cluster tree at various resolutions 
library(clustree)
p4 <- clustree(cncc, prefix = "RNA_snn_res.")
p4
Idents(cncc) <- "RNA_snn_res.0.4"
p1 <- DimPlot(cncc, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(cncc, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(cncc,group.by = "stage")+NoAxes()
pdf("figures/cncc_sub_all_final.pdf",width = 16,height = 8)
p1+p2
p4
dev.off()

#find marker genes for major cell types
sig_genes <- FindAllMarkers(cncc, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)

# #remove "LIN" genes from markers lists
# LINFeatures <- grep("LIN",rownames(cncc))
# cncc1 <- rownames(cncc)[-LINFeatures]
# sig_genes2 <- filter(sig_genes, gene %in% cncc1)
# 
# top <- sig_genes2 %>% group_by(cluster) %>% top_n(20, avg_log2FC);top
# p <- DotPlot(object = cncc, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
# df <- data.frame(Gene = p$data$features.plot, avg.exp = p$data$avg.exp.scaled, pct.exp = p$data$pct.exp, cluster = p$data$id)
# p2 <- df %>% 
#   filter(avg.exp > 0, pct.exp > 1) %>% 
#   ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
#   geom_point() +
#   scale_color_viridis_c() + 
#   cowplot::theme_cowplot() + 
#   theme(axis.line  = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   ylab('') + xlab("")+
#   theme(axis.ticks = element_blank()) 
# p2

# #annotate cncc sub-clusters res0.8
Idents(cncc) <- "RNA_snn_res.0.8"
cncc <- RenameIdents(cncc, c("0"="eCNCC","1"="iCNCC","2"="ilCNCC1","3"="lCNCC1","4"="ilCNCC1", "5"="ilCNCC2","6"="lCNCC2",
                           "7"="ilCNCC3", "8"="eCNCC", "9"="optic","10"="lCNCC3","11"="ncALX4",
                           "12"="pigment"))
cncc$anno_res0.8 <- Idents(cncc)

#annotate cncc sub-clusters res0.4
Idents(cncc) <- "RNA_snn_res.0.4"
cncc <- RenameIdents(cncc, c("0"="iCNCC1","1"="eCNCC1","2"="lCNCC1","3"="iCNCC2","4"="iCNCC3", "5"="lCNCC2","6"="iCNCC4",
                             "7"="cnl1", "8"="cnl2", "9"="cnl3","10"="cnl4"))
cncc$anno_res0.4 <- Idents(cncc)
saveRDS(cncc,"data/cncc_final.rds")

DimPlot(mes)
cncc$tempanno <- Idents(cncc)
ect$tempanno <- ect$anno_res1
mes$tempanno <- paste0("mes",Idents(mes))

merged <- merge(mes, c(ect,cncc))
Idents(merged) <- "tempanno"
for (i in levels(merged)) {
  id <- rownames(merged@meta.data[merged@meta.data$tempanno == i,])
  cds2 <- SetIdent(cds2,id,value = i)
}        
DimPlot(cds2)
cds2$tempanno <- Idents(cds2)
sce <- as.SingleCellExperiment(cds2)
saveRDS(sce, "data/tempanno_sce.rds")
sce

ect@meta.data
met <- ect@meta.data[,c(1:26,28:29)]

ect1 <- ect

ect@meta.data <- met
saveRDS(ect1, "data/ect_final.rds")

for (i in levels(mes)) {
  id <- rownames(mes@meta.data[Idents(mes) == i,])
  cds2 <- SetIdent(cds2,id,value = i)
}        

cds2$tempanno <- Idents(cds2)
