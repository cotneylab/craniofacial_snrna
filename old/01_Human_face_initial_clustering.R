library(Seurat)
library(tidyverse)
library(harmony)
#loading inHouse data
counts <- Read10X_h5("CS12_1_matrix.h5")
CS12_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS12_1")
counts <- Read10X_h5("CS12_2_matrix.h5")
CS12_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS12_2")
counts <- Read10X_h5("CS12_3_matrix.h5")
CS12_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS12_3")
counts <- Read10X_h5("CS12_4_matrix.h5")
CS12_4 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS12_4")
counts <- Read10X_h5("CS13_1_matrix.h5")
CS13_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS13_1")
counts <- Read10X_h5("CS13_2_matrix.h5")
CS13_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS13_2")
counts <- Read10X_h5("CS13_3_matrix.h5")
CS13_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS13_3")
counts <- Read10X_h5("CS13_4_matrix.h5")
CS13_4 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS13_4")
counts <- Read10X_h5("CS13_5_matrix.h5")
CS13_5 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS13_5")
counts <- Read10X_h5("CS13_6_matrix.h5")
CS13_6 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS13_6")
counts <- Read10X_h5("CS14_1_matrix.h5")
CS14_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS14_1")
counts <- Read10X_h5("CS14_2_matrix.h5")
CS14_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS14_2")
counts <- Read10X_h5("CS14_3_matrix.h5")
CS14_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS14_3")
counts <- Read10X_h5("CS16_1_matrix.h5")
CS16_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS16_1")
counts <- Read10X_h5("CS16_2_matrix.h5")
CS16_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS16_2")
counts <- Read10X_h5("CS16_3_matrix.h5")
CS16_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS16_3")
counts <- Read10X_h5("CS17_1_matrix.h5")
CS17_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS17_1")
counts <- Read10X_h5("CS17_2_matrix.h5")
CS17_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS17_2")
counts <- Read10X_h5("CS17_3_matrix.h5")
CS17_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS17_3")
counts <- Read10X_h5("CS20_1_matrix.h5")
CS20_1 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS20_1")
counts <- Read10X_h5("CS20_2_matrix.h5")
CS20_2 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS20_2")
counts <- Read10X_h5("CS20_3_matrix.h5")
CS20_3 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS20_3")
counts <- Read10X_h5("CS20_4_matrix.h5")
CS20_4 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS20_4")
counts <- Read10X_h5("CS20_5_matrix.h5")
CS20_5 <- CreateSeuratObject(counts$`Gene Expression`,assay = "RNA",project = "CS20_5")
cds <- merge(CS12_1,c(CS12_2,CS12_3,CS12_4,CS13_1,CS13_2,CS13_3,CS13_4,CS13_5,CS13_6,CS14_1,CS14_2,CS14_3,CS16_1,CS16_2,CS16_3,CS17_1,CS17_2,CS17_3,CS20_1,CS20_2,CS20_3,CS20_4,CS20_5))
cds$data <- "inHouse"
#Filtration and QC
cds[["percent.mt"]] <- PercentageFeatureSet(cds, pattern = "^MT-")
p1 <- VlnPlot(object = cds, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
cds.filt <- subset(cds, subset = nFeature_RNA < 7000 & nCount_RNA < 25000 & percent.mt < 10 & nCount_RNA > 500 & nFeature_RNA > 500)
p2 <- VlnPlot(object = cds.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p3 <- FeatureScatter(object = cds, feature1 = "nCount_RNA", feature2 = "percent.mt")
p4 <- FeatureScatter(object = cds, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p5 <- FeatureScatter(object = cds.filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
p6 <- FeatureScatter(object = cds.filt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dir.create("figures")
pdf("figures/QC.pdf",w=18,height = 6)
p1+plot_annotation("Before Filtration")
p2+plot_annotation("After Filtration")
(p3+p4)+plot_annotation("Before Filtration")
(p5+p6)+plot_annotation("After Filtration")
dev.off()
dim(cds)
#[1] 36601 86359
dim(cds.filt)
#[1] 36601 77653
cds <- cds.filt
#inHouse data processing and clustering
cds <- JoinLayers(cds)

require(biomaRt)
mart <-  useMart("ENSEMBL_MART_ENSEMBL", host = "https://useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = rownames(cds),
  uniqueRows=TRUE)

annotLookup1 <- annotLookup %>% drop_na()
dim(annotLookup)
#23954     4
dim(annotLookup1)
#21811     4
features <- unique(intersect(rownames(cds),annotLookup1$hgnc_symbol))
to.filter <- c(grep("MT-",rownames(cds),value = T),grep("-AS",rownames(cds),value = T))
features1 <- features[!(features %in% to.filter)]
cds2 <- subset(cds, features = unique(features1))
cds2 <- NormalizeData(cds2)
cds2 <- CellCycleScoring(cds2, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
cds2 <- FindVariableFeatures(cds2) %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")
cds2 <- RunUMAP(cds2, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 1, verbose = FALSE)
cds2@meta.data <- separate(cds2@meta.data ,col = "orig.ident" ,sep = "_",into = c("stage","rep"),remove = F)
cds2$stage <- factor(cds2$stage, levels = c("CS12","CS13","CS14","CS16","CS17","CS20"))
p1 <- DimPlot(cds2, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(cds2, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(cds2,split.by = "stage",label = T)+NoAxes()+NoLegend()
p1+p2
p3
pdf("figures/UMAP_harmony.pdf",width = 16,height = 8)
p1+p2
dev.off()

pdf("figures/UMAP_byStage.pdf",width = 16,height = 4)
p3
dev.off()


cds2 <- FindClusters(cds2,resolution = c(0.1,0.2,0.4,0.6,0.8,1), verbose = FALSE)
Idents(cds2) <- cds2$RNA_snn_res.0.1
p <- DimPlot(cds2,group.by = "orig.ident")
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]] # Pull the data used for the plot
pdata$sample <- cds2$stage
pdata <-  pdata[order(pdata$sample), ] # Order the plot data by group
ucols <- unique(pdata$colour) # Get a vector of unique colors
naect(ucols) <- unique(pdata$sample) #
barplot(prop.table(table( cds2$orig.ident,Idents(cds2)), margin=2), col = ucols, legend.text = T)
FeaturePlot(cds2, "TUBB3")
DimPlot(cds2, label = T)
# Filter out Neural clusters
cds2 <- subset(cds2, idents = c(1,2,4,5,8), invert=T)
cds2 <- NormalizeData(cds2) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")

cds2 <- RunUMAP(cds2, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 1, verbose = FALSE)
p1 <- DimPlot(cds2, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(cds2, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(cds2,split.by = "stage",label = T)+NoAxes()+NoLegend()
p1+p2
p3

pdf("figures/UMAP_harmony_NoNeur.pdf",width = 16,height = 8)
p1+p2
dev.off()
pdf("figures/UMAP_byStage_NoNeur.pdf",width = 20,height = 4)
p3
dev.off()
cds2 <- FindClusters(cds2,resolution = c(0.1,0.2,0.4,0.6,0.8,1), verbose = FALSE)

# lower resolution clustering to label 
Idents(cds2) <- cds2$RNA_snn_res.0.2
DimPlot(cds2, label = T)+NoAxes()+NoLegend()
markers <- str_to_upper(c("Prrx1","Alx1","Epcam","Ttn","Kel","Hemgn","Cdh5","Sox10","Fcer1g"))
cds2 <- RunTSNE(cds2,reduction = "harmony",dims = 1:30)
p <- FeaturePlot(cds2, markers,combine = F,reduction = "tsne",cols = c("grey90","red"))
for (i in 1:length(markers)) {
  p[[i]] <- p[[i]]+NoAxes()
}
pdf("figures/markers_MajorCellTypes_tsne.pdf",w=16,h=16)
patchwork::wrap_plots(p)
dev.off()
saveRDS(cds2,"data/cds2.rds")
cds2 <- RenameIdents(cds2, c("0"="mesenchyme","1"="mesenchyme","2"="ectoderm","3"="mesenchyme","4"="muscle","5"="unknown","6"="endothelium",
                           "7"="ectoderm","8"="red blood cells","9"="mesenchyme","10"="other blood cells","11"="mesenchyme"))
levels(cds2) <- c("mesenchyme","ectoderm","muscle","red blood cells","endothelium","unknown","other blood cells")
cds2 <- RenameIdents(cds2,c("unknown"="CNCC"))
levels(cds2) <- c("mesenchyme","ectoderm","muscle","red blood cells","endothelium","CNCC","other blood cells")

cds2$cell_type1 <- Idents(cds2)
#find marker genes for major cell types
sig_genes <- FindAllMarkers(cds2, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
sig_genes$TF <- "no"
id <- intersect(annotLookup$hgnc_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$Entrez <- annotLookup$entrezgene_id[match(as.character(sig_genes$gene),annotLookup$hgnc_symbol)]
sig_genes$description <- annotLookup$description[match(as.character(sig_genes$gene),annotLookup$hgnc_symbol)]
sig_genes <- sig_genes[, c("gene","TF","Entrez","description","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
openxlsx::write.xlsx(sig_genes,"tables/SigGenes_MajorCellTypes.xlsx")

top <- sig_genes %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
p <- DotPlot(object = cds2, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
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
dir.create("tables/majorcellTypes_stages")
for (i in levels(cds2)) {
  seu_temp <- subset(cds2, idents = i)
  Idents(seu_temp) <- "stage"
  sig_genes <- FindAllMarkers(seu_temp, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                              min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                              only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                              latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                              pseudocount.use = 1, return.thresh = 0.01)
  sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
  openxlsx::write.xlsx(sig_genes,paste0("tables/majorcellTypes_stages/SigGenes_",i,".xlsx"))
  
}


# CNCC module score doi: 10.1242/dev.105445
# https://www.sigmaaldrich.com/US/en/product/mm/scc049
# https://www.liebertpub.com/doi/pdf/10.1089/scd.2012.0155
#cnccMarkers <- str_to_upper(c("Sox10","Foxd3","Lmo4","Sox5","Sox6","Tfap2a","Ebf1","Pax3", "Rxrg","Twist1", "Myc", "Ets1", "Cald1","Rac1","Cfl1", "Crabp1", "Epha2","Itgb1", "Cd44","Dlx1", "Dlx2","Dlx5", "Dlx6","Sox9","Id1","Id2","Id3","Id4","Myc","Myb","Mef2c"))
cnccMarkers <- str_to_upper(c("Ets1", "Foxd3", "Sox10", "Tfap2a", "NR2F1", "NR2F2"))

cds1 <- cds2
cds1 <- AddModuleScore(cds1, list(cnccMarkers),name = "CNCC_ModuleScore" )
cds2$CNCC_ModuleScore <- cds1$CNCC_ModuleScore1
p1 <- FeaturePlot(cds2, "CNCC_ModuleScore",cols = c("grey90","red"))+NoAxes()
p2 <- FeaturePlot(cds2, "CNCC_ModuleScore",cols = c("grey90","red"),reduction = "tsne")+NoAxes()
pdf("figures/CNCC_ModuleScore_UMAP_TSNE.pdf",w=6,h=6)
p1
p2
dev.off()
