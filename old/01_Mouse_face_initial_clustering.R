library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
#loading inHouse data
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E9_1_matrix.h5")
E9_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E9_1")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E9_2_matrix.h5")
E9_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E9_2")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E9_3_matrix.h5")
E9_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E9_3")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E10_1_matrix.h5")
E10_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E10_1")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E10_2_matrix.h5")
E10_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E10_2")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E10_3_matrix.h5")
E10_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E10_3")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E11_1_matrix.h5")
E11_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E11_1")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E11_2_matrix.h5")
E11_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E11_2")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E11_3_matrix.h5")
E11_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E11_3")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E12_1_matrix.h5")
E12_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E12_1")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E12_2_matrix.h5")
E12_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E12_2")
counts <- Read10X_h5("/Volumes/Extreme SSD/Mouse_cellranger/E12_3_matrix.h5")
E12_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "E12_3")
cds <- merge(E10_1,c(E10_2,E10_3,E11_1,E11_2,E11_3,E12_1,E12_2,E12_3))
cds$data <- "inHouse"
cds@meta.data <- separate(cds@meta.data ,col = "orig.ident" ,sep = "_",into = c("stage","rep"),remove = F)

# Combine inHouse data with later stages from DSouza
e13 <- Read10X("/Users/naghamkhourifarah/Downloads/GSE245469_RAW/scRNA/E13")
e13 <- CreateSeuratObject(e13,assay = "RNA",project = "E13_1")
e13$data <- "DSouza"
e13$stage <- "E13"
e14 <- Read10X("/Users/naghamkhourifarah/Downloads/GSE245469_RAW/scRNA/E14")
e14 <- CreateSeuratObject(e14,assay = "RNA",project = "E14_1")
e14$data <- "DSouza"
e14$stage <- "E14"
e15 <- Read10X("/Users/naghamkhourifarah/Downloads/GSE245469_RAW/scRNA/E15")
e15 <- CreateSeuratObject(e15,assay = "RNA",project = "E15_1")
e15$data <- "DSouza"
e15$stage <- "E15"

e13m <- Read10X("/Users/naghamkhourifarah/Downloads/GSE245469_RAW/Multiome/E13_wt")
e13m <- CreateSeuratObject(e13m$`Gene Expression`,assay = "RNA",project = "E13_2")
e13m$data <- "DSouza"
e13m$stage <- "E13"

e15m <- Read10X("/Users/naghamkhourifarah/Downloads/GSE245469_RAW/Multiome/E15_wt")
e15m <- CreateSeuratObject(e15m$`Gene Expression`,assay = "RNA",project = "E15_2")
e15m$data <- "DSouza"
e15m$stage <- "E15"

cds2 <- merge(cds,c(e13,e13m,e14,e15,e15m))
#Filtering and QC
Idents(cds2) <- "orig.ident"
cds2[["percent.mt"]] <- PercentageFeatureSet(cds2, pattern = "^mt-")
p1 <- VlnPlot(object = cds2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
cds.filt <- subset(cds2, subset = nFeature_RNA < 7000 & nCount_RNA < 20000 & percent.mt < 10 & nCount_RNA > 500 & nFeature_RNA > 500)
p2 <- VlnPlot(object = cds.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p3 <- FeatureScatter(object = cds2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p4 <- FeatureScatter(object = cds2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p5 <- FeatureScatter(object = cds.filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
p6 <- FeatureScatter(object = cds.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("figures/QC_E10toE15.pdf",w=12,height = 8)
p1+plot_annotation("Before Filtration")
p2+plot_annotation("After Filtration")
(p3+p4)+plot_annotation("Before Filtration")
(p5+p6)+plot_annotation("After Filtration")
dev.off()

cds2$stage <- factor(cds2$stage, levels = c("E10","E11","E12","E13","E14","E15"))
cds2 <- JoinLayers(cds2)
require(biomaRt)
mart_m <- useMart("ENSEMBL_MART_ENSEMBL", host = "useast.ensembl.org")
mart_m <- useDataset("mmusculus_gene_ensembl", mart_m)
annotLookup2 <- getBM(
  mart = mart_m,
  attributes = c(
    "mgi_symbol","description",
    "entrezgene_id",
    "ensembl_gene_id",
    "gene_biotype"),
  filter = "mgi_symbol",
  values = rownames(cds2),
  uniqueRows=TRUE)
dim(annotLookup2)

annotLookup <- annotLookup2 %>% drop_na()
dim(annotLookup)
MitoFeatures <- grep("mt-",rownames(cds2))
count <- GetAssayData(cds2,assay = "RNA",slot = "count")
count <- count[-MitoFeatures,]
cds2 <- subset(cds2, features = intersect(rownames(count),annotLookup$mgi_symbol))
#processing and clustering of combined data
cds2 <- NormalizeData(cds2) %>% 
  FindVariableFeatures() 
cds2 <- CellCycleScoring(cds2, s.features = str_to_sentence(cc.genes$s.genes), g2m.features = str_to_sentence(cc.genes$g2m.genes), set.ident = F)
cds2 <- ScaleData(cds2,vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")

cds2 <- RunUMAP(cds2, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1), verbose = FALSE)

p1 <- DimPlot(cds2, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(cds2, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(cds2,split.by = "stage",label = T)+NoAxes()+NoLegend()
p4 <- DimPlot(cds2,group.by = "Phase",label = T)+NoAxes()+NoLegend()
p1+p2
p3

pdf("figures/UMAP_harmony_E10_E15.pdf",width = 16,height = 8)
p1+p2
dev.off()
pdf("figures/UMAP_byStage_E10_E15.pdf",width = 20,height = 4)
p3
dev.off()

dir.create("data")
saveRDS(cds2,"data/E10_E15.rds")


#Filter neural cluster
Idents(cds2) <-  "RNA_snn_res.0.1"
DimPlot(cds2, label = T)+NoAxes()+NoLegend()
FeaturePlot(cds2, "Tubb3")
cds <- subset(cds2, idents = "2", invert=T)
cds <- NormalizeData(cds) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")
ElbowPlot(cds,ndims = 50)
cds <- RunUMAP(cds, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1), verbose = FALSE)
p1 <- DimPlot(cds, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(cds, group.by = "stage")+NoAxes()
p3 <- DimPlot(cds,split.by = "stage",label = T)+NoAxes()+NoLegend()
p4 <- DimPlot(cds, label = T, group.by = "RNA_snn_res.0.1")+NoAxes()+NoLegend()
p5 <- DimPlot(cds, group.by = "stage")+NoAxes()
p1+p2
p4+p5
p3
pdf("figures/UMAP_harmony_E10_E15_NoNeur.pdf",width = 16,height = 8)
p1+p2
p4+p5
dev.off()
pdf("figures/UMAP_byStage_E10_E15_NoNeur.pdf",width = 20,height = 4)
p3
dev.off()
library(clustree)
clustree(cds, prefix = "RNA_snn_res.")

#find marker genes for major cell types
Idents(cds) <-  "RNA_snn_res.0.2"
DimPlot(cds, label = T)+NoAxes()+NoLegend()
mGenes <- readRDS(file="/Users/naghamkhourifarah/Documents/Annotation1.rds")
sig_genes <- FindAllMarkers(cds, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
sig_genes$TF <- "no"
id <- intersect(mGenes$mgi_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$TF[sig_genes$gene %in% id] <- mGenes$Type[match(x$gene,mGenes$mgi_symbol)]
sig_genes$Entrez <- mGenes$entrezgene[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes$description <- mGenes$description[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes <- sig_genes[, c("gene","TF","Entrez","description","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
openxlsx::write.xlsx(sig_genes,"tables/SigGenes_E10_E15_res0.1.xlsx")

#remove "Rik" genes from markers lists
# RikFeatures <- grep("Rik",rownames(count))
# count1 <- count[-RikFeatures,]
# sig_genes2 <- filter(sig_genes, gene %in% row.names(count1))
top <- sig_genes %>% group_by(cluster) %>% top_n(50, avg_log2FC) %>% top_n(5, pct.diff);top
p <- DotPlot(object = cds, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
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

#Annotate major cell types
DimPlot(cds, label = T)+NoAxes()+NoLegend()
Idents(cds) <-  "RNA_snn_res.0.2"
DimPlot(cds, label = T, group.by = "cellType")+NoAxes()+NoLegend()
cds <- RenameIdents(cds, c("0"="mesenchyme","1"="mesenchyme","2"="mesenchyme","4"="mesenchyme","9"="mesenchyme","11"="mesenchyme",
                           "7"="muscle","5"="red blood cells", "8"="unknown","3"="ectoderm","12"="ectoderm","6"="endothelium",
                           "10"="other blood cells","13"="other blood cells"))
levels(cds) <- c("mesenchyme","ectoderm","muscle","red blood cells","endothelium","unknown","other blood cells")
cds$cell_type1 <- Idents(cds)

# CNCC module score doi: 10.1242/dev.105445
# https://www.sigmaaldrich.com/US/en/product/mm/scc049
# https://www.liebertpub.com/doi/pdf/10.1089/scd.2012.0155
cnccMarkers <- c("Sox10","Foxd3","Lmo4","Sox5","Sox6","Tfap2a","Ebf1","Pax3", "Rxrg","Twist1", "Myc", "Ets1", "Cald1","Rac1","Cfl1", "Crabp1", "Epha2","Itgb1", "Cd44",
                 "Dlx1", "Dlx2","Dlx5", "Dlx6","Sox9","Id1","Id2","Id3","Id4","Myc","Myb","Mef2c")
FeaturePlot(cds, c("Sox10","Foxd3","Lmo4","Sox5","Sox6","Tfap2a","Ebf1","Pax3", "Rxrg","Twist1", "Myc", "Ets1", "Cald1","Rac1","Cfl1", "Crabp1", "Epha2","Itgb1", "Cd44",
                   "Dlx1", "Dlx2","Dlx5", "Dlx6","Sox9","Id1","Id2","Id3","Id4","Myc","Myb","Mef2c"))
cds1 <- cds
cds1 <- AddModuleScore(cds1, cnccMarkers,name = "CNCC_ModuleScore" )
cds$CNCC_ModuleScore <- cds1$CNCC_ModuleScore1
cds <- RunTSNE(cds, dims = 1:30)
p1 <- FeaturePlot(cds, "CNCC_ModuleScore",cols = c("grey90","red"))+NoAxes()
p2 <- FeaturePlot(cds, "CNCC_ModuleScore",cols = c("grey90","red"),reduction = "tsne")+NoAxes()
pdf("figures/CNCC_ModuleScore_E10_E15.pdf",w=6,h=6)
p1
p2
dev.off()
cds <- RenameIdents(cds, c("unknown"="CNCC"))
levels(cds) <- c("mesenchyme","ectoderm","muscle","red blood cells","endothelium","CNCC","other blood cells")
cds$cell_type <- Idents(cds)
p2 <- DimPlot(cds, label = T, reduction="tsne")+NoAxes()+NoLegend()
p1 <- DimPlot(cds, label = T)+NoAxes()+NoLegend()
pdf("figures/UMAP_TSNE_MajorCellTypes_E10_E15.pdf",w=6,h=6)
p1
p2
dev.off()

saveRDS(cds,"data/E10_E15_filtered.rds")
#Find markers for the annotated major cell types
sig_genes <- FindAllMarkers(cds, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
sig_genes$TF <- "no"
id <- intersect(mGenes$mgi_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$TF[sig_genes$gene %in% id] <- mGenes$Type[match(x$gene,mGenes$mgi_symbol)]
sig_genes$Entrez <- mGenes$entrezgene[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes$description <- mGenes$description[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes <- sig_genes[, c("gene","TF","Entrez","description","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
openxlsx::write.xlsx(sig_genes,"tables/SigGenes_E10_E15_MajorCellTypes.xlsx")
#remove "Rik" genes from markers lists for plotting
RikFeatures <- grep("Rik",rownames(count))
count1 <- count[-RikFeatures,]
sig_genes2 <- filter(sig_genes, gene %in% row.names(count1))
top <- sig_genes2 %>% group_by(cluster) %>% top_n(10, avg_log2FC)# %>% top_n(10, pct.diff);top
p <- DotPlot(object = cds, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
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

markers <- c("Prrx1","Alx1","Epcam","Ttn","Kel","Hemgn","Cdh5","Sox10","Fcer1g")
p <- FeaturePlot(cds, markers,combine = F,reduction = "tsne",cols = c("grey90","red"))
for (i in 1:length(markers)) {
  p[[i]] <- p[[i]]+NoAxes()
}
pdf("figures/markers_MajorCellTypes_E10_E15_tsne.pdf",w=16,h=16)
patchwork::wrap_plots(p)
dev.off()
saveRDS(cds,"data/cds.rds")

