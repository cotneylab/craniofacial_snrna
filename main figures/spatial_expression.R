library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(DropletUtils)
library(tidyverse)
library(harmony)
s1_filter_matrix <- Read10X("./s1_filter_mx/")
s2_filter_matrix <- Read10X("./s2_filter_mx/")


write10xCounts("./s1_filter_mx/filtered_feature_bc_matrix.h5", s1_filter_matrix, type = "HDF5",
               genome = "hg38", version = "3", overwrite = TRUE,
               gene.id = rownames(s1_filter_matrix),
               gene.symbol = rownames(s1_filter_matrix))

write10xCounts("./s2_filter_mx/filtered_feature_bc_matrix.h5", s2_filter_matrix, type = "HDF5",
               genome = "hg38", version = "3", overwrite = TRUE,
               gene.id = rownames(s2_filter_matrix),
               gene.symbol = rownames(s2_filter_matrix))

s1_counts <- Seurat::Read10X("s1_filter_mx")  
s2_counts <- Seurat::Read10X("s2_filter_mx") 


s1_data = Seurat::CreateSeuratObject(
  counts = s1_counts , 
  project = 's1_spatial', 
  assay = 'Spatial')

s1_data$slice = 1 
s1_data$region = 'cs13'
s1_imgpath = "s1_filter_mx/spatial"
s1_img = Seurat::Read10X_Image(image.dir = s1_imgpath, image.name = "tissue_lowres_image.png", assay = "Spatial", slice = "slice1")  
Seurat::DefaultAssay(object = s1_img) <- 'Spatial'
s1_img = s1_img[colnames(x = s1_data)]  
s1_data[['image1']] = s1_img 
s1_data[["percent.mt"]] <- PercentageFeatureSet(s1_data, pattern = "^MT-")

s2_data = Seurat::CreateSeuratObject(
  counts = s2_counts , 
  project = 's2_spatial', 
  assay = 'Spatial')

s2_data$slice = 1 
s2_data$region = 'cs13'

s2_imgpath = "s2_filter_mx/spatial"
s2_img = Seurat::Read10X_Image(image.dir = s2_imgpath, image.name = "tissue_lowres_image.png", assay = "Spatial", slice = "slice2")  
Seurat::DefaultAssay(object = s2_img) <- 'Spatial'

s2_img = s2_img[colnames(x = s2_data)]  
s2_data[['image2']] = s2_img  
s2_data[["percent.mt"]] <- PercentageFeatureSet(s2_data, pattern = "^MT-")

visium_list <- list(s1_data,s2_data)
#need to increase size for SCTransform
options(future.globals.maxSize = 8000 * 1024^2)
visium_list <- lapply(X = visium_list , FUN = SCTransform, assay = "Spatial", vars.to.regress = "percent.mt", method = "glmGamPoi", return.only.var.genes =FALSE, verbose = TRUE)
integration_features <- SelectIntegrationFeatures(visium_list)

s1_data <- SCTransform(s1_data, assay = "Spatial", vars.to.regress = "percent.mt", method = "glmGamPoi", return.only.var.genes =FALSE, verbose = TRUE)
s2_data <- SCTransform(s2_data, assay = "Spatial", vars.to.regress = "percent.mt", method = "glmGamPoi", return.only.var.genes =FALSE, verbose = TRUE)

#merge_data
emb.merge <- merge(s1_data, s2_data)
VariableFeatures(emb.merge) <- integration_features
plot1 <- VlnPlot(emb.merge, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(emb.merge, features = "nCount_Spatial", image.alpha = 0.3) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
pdf("initial_qc.pdf")
plot1
plot2
dev.off()
emb.merge  <- PrepSCTFindMarkers(emb.merge)
emb.merge <- SCTransform(emb.merge, assay = "Spatial", vars.to.regress = "percent.mt", method = "glmGamPoi", return.only.var.genes =FALSE, verbose = TRUE)
emb.merge <- RunPCA(emb.merge, assay = "SCT")
print(emb.merge[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(emb.merge, reduction = "pca") + NoLegend()
DimHeatmap(emb.merge, dims = 1, cells = 500, balanced = TRUE)
emb.merge <- RunUMAP(emb.merge, dims = 1:20)
p1 <- DimPlot(emb.merge, group.by = "orig.ident", label = TRUE)
p2 <- SpatialDimPlot(emb.merge, label = TRUE, label.size = 3,  image.alpha = 0.3)
p1 + p2
#harmonize
emb.merge <- RunHarmony(emb.merge, assay.use="SCT", group.by.vars = "orig.ident")
emb.merge <- RunUMAP(emb.merge, dims = 1:20, reduction = "harmony")
emb.merge <- FindNeighbors(emb.merge, reduction = "harmony", dims = 1:20)
emb.merge <- FindClusters(emb.merge, resolution = 0.8)

p3 <- DimPlot(emb.merge, group.by = "orig.ident", label = TRUE)
p4<- DimPlot(emb.merge, group.by = "ident", split.by = 'orig.ident')
p5 <- DimPlot(emb.merge, reduction = "umap", label = TRUE)
p6 <- SpatialDimPlot(emb.merge, label = TRUE, label.size = 3,  image.alpha = 0.3)
p7 <- SpatialDimPlot(emb.merge, label = FALSE, label.size = 3,  image.alpha = 0.3)
p8 <- VlnPlot(emb.merge, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p9 <- SpatialFeaturePlot(emb.merge, features = "nCount_Spatial", image.alpha = 0.3) + theme(legend.position = "right")

pdf(file = "harmonize_cluster_qc.pdf")
p1
p2
p3
p4
p5
p6
p7
p8
p9
dev.off()

sig_genes <- FindAllMarkers(emb.merge, assay = "SCT" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)

sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)

require(biomaRt)
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = rownames(emb.merge),
  uniqueRows=TRUE)

id <- intersect(annotLookup$hgnc_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$Entrez <- annotLookup$entrezgene_id[match(as.character(sig_genes$gene),annotLookup$hgnc_symbol)]
sig_genes$description <- annotLookup$description[match(as.character(sig_genes$gene),annotLookup$hgnc_symbol)]
sig_genes$biotype <- annotLookup$gene_biotype[match(as.character(sig_genes$gene),annotLookup$hgnc_symbol)]

sig_genes <- sig_genes[, c("gene","Entrez","description","biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
openxlsx::write.xlsx(sig_genes,"SigGenes_Main_Clusters.xlsx")

sig_genes_protein <- subset(sig_genes, biotype == "protein_coding")

top <- sig_genes_protein %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top100 <- sig_genes_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

library("clusterProfiler")
library(enrichplot)
library("org.Hs.eg.db")
library("AnnotationHub")
library(ReactomePA)
library(DOSE)

#gene ontology enrichment across clusters
#need to use Entrez id in clusterprofiler
df1 <- top100pval[,c(2,10)]
df1sample <- split(df1$Entrez,df1$cluster)
length(df1sample)
main_types <- names(df1sample)
genelist <- as.list(df1sample)

BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)


DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)


options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 10, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 10, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 10, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 10, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 10, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 10, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))


cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "main_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
dev.off()

go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "gene_onology_bp_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "gene_onology_cc_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "gene_onology_mf_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

#rename based on functional enrichments, position, and original labels
#use code to isolate specific clusters.
#SpatialDimPlot(emb.merge, cells.highlight = CellsByIdentities(object = emb.merge, idents = c(5)), facet.highlight = TRUE, ncol = 3, image.alpha = 0.3)


emb.merge <- RenameIdents(emb.merge, c("0"="posterior neural crest","1"="trunk neural crest","2"="somatic LPM","3"="somites","4"="IM","5"="neuron","6"="heart",
                     "7"="head mesoderm","8"="posterior IM","9"="trunk IM","10"="somites","11"="blood", "12" = "lower digestive", "13" = "phyarngeal arches",
                     "14" = "blood", "15" = "hindbrain", "16" = "forebrain","17" = "endothelium", "18" = "liver", "19" = "craniofacial_ectoderm", "20" = "upper digestive", "21" = "epidermis"))


pdf("cluster_naming_check.pdf", h = 8.5, w = 11)
SpatialDimPlot(emb.merge, label = TRUE, label.size = 3,  image.alpha = 0.3) + NoLegend()
dev.off()

saveRDS(emb.merge, file = "emb.merge.rds")

#regenerate markers based on new labels
sig_genes <- FindAllMarkers(emb.merge, assay = "SCT" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)

sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$hgnc_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$Entrez <- annotLookup$entrezgene_id[match(as.character(sig_genes$gene),annotLookup$hgnc_symbol)]
sig_genes$description <- annotLookup$description[match(as.character(sig_genes$gene),annotLookup$hgnc_symbol)]
sig_genes$biotype <- annotLookup$gene_biotype[match(as.character(sig_genes$gene),annotLookup$hgnc_symbol)]

sig_genes <- sig_genes[, c("gene","Entrez","description","biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
openxlsx::write.xlsx(sig_genes,"SigGenes_Main_Clusters.xlsx")

sig_genes_protein <- subset(sig_genes, biotype == "protein_coding")

top <- sig_genes_protein %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top100 <- sig_genes_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)
df1 <- top100pval[,c(2,10)]
df1sample <- split(df1$Entrez,df1$cluster)
length(df1sample)
main_types <- names(df1sample)
genelist <- as.list(df1sample)

go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "gene_onology_bp_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "gene_onology_cc_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "gene_onology_mf_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)


cnccMarkers_fulllist <- str_to_upper(c("Arpc1b", "Barx1", "Bmi1", "Cald1", "Cd44", "Cd47", "Cd9", "Cfl1", "Col1a1", "Crabp1", "Cyba", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Ebf1", "Epha2", "Ets1", "Foxd3", "Gli3", "Hk2", "Id1", "Id2", "Id3", "Id4", "Itgb1", "Klf4", "Lmo4", "Mef2c", "Msh6", "Msx1", "Myb", "Myc", "Myo10", "Pax3", "Pax7", "Prrx1", "Rac1", "Rarg", "Runx2", "Rxrg", "Snai1", "Snai2", "Sox10", "Sox5", "Sox6", "Sox9", "Spp1", "Tfap2a", "Twist1"))
#smaller list https://www.nature.com/articles/s42003-021-01970-0
cnccMarkers <- str_to_upper(c("Ets1", "Foxd3", "Tfap2a", "Tfap2b", "NR2F1", "NR2F2", "Tcf3", "Tcf12"))
migratory_cnccMarkers <- str_to_upper(c("Pax7", "Tfap2a"))
premigratory_cnccMarkers <- str_to_upper(c("Sox9", "Sox10", "Pax3"))

hescvcnccdegenes <- read.table(file ="../res_h9ncc.txt", header = TRUE, row.names = 1, sep = "\t")
yankeecnccMarkers <- subset(hescvcnccdegenes, padj < 0.05)
yankeecnccMarkers <- subset(yankeecnccMarkers, log2FoldChange > 4)
yankeecnccMarkers <- subset(yankeecnccMarkers, Symbol != "NA")
yankeeCNCC <- yankeecnccMarkers$Symbol
alx <- c("ALX1", "ALX3", "ALX4")
treachercollinsgenes <- c("TCOF1", "POLR1C", "POLR1D")

human_main_markers <- readxl::read_xlsx("tables/SigGenes_Main_CellTypes.xlsx", sheet = 1)
human_subtype_markers <- readxl::read_xlsx("tables/SigGenes_CellSubTypes.xlsx", sheet = 1)
human_mes_subtype_markers <- readxl::read_xlsx("tables/SigGenes_Cell_mes_SubTypes.xlsx", sheet = 1)
human_cncc_subtype_markers <- readxl::read_xlsx("tables/SigGenes_Cell_cncc_SubTypes.xlsx", sheet = 1)
human_ect_subtype_markers <- readxl::read_xlsx("tables/SigGenes_Cell_ect_SubTypes.xlsx", sheet = 1)
human_subtype_markers <- subset(human_subtype_markers, biotype == "protein_coding")
top <- human_subtype_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top100 <- human_subtype_markers %>% group_by(cluster) %>% top_n(100, avg_log2FC);top
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

craniofacial_subtypes <- unique(top100pval$cluster)

extract_subtype_genelists <- function(x) {
  subtype_genelists[[x]] <- subset(top, cluster == x)$gene
}

subtype_genelists <- list()

subtype_genelists <- lapply(craniofacial_subtypes, extract_subtype_genelists)
names(subtype_genelists) <- craniofacial_subtypes

top <- human_mes_subtype_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top50 <- human_mes_subtype_markers %>% group_by(cluster) %>% top_n(50, avg_log2FC);top
top50pval <- subset(top50, rowSums(top50[5] < 0.05) > 0)
mes_subtype_marker_genes <- unique(top50pval$gene)
top <- human_ect_subtype_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top50 <- human_ect_subtype_markers %>% group_by(cluster) %>% top_n(50, avg_log2FC);top
top50pval <- subset(top50, rowSums(top50[5] < 0.05) > 0)
ect_subtype_marker_genes <- unique(top50pval$gene)
top <- human_cncc_subtype_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top50 <- human_cncc_subtype_markers %>% group_by(cluster) %>% top_n(50, avg_log2FC);top
top50pval <- subset(top50, rowSums(top50[5] < 0.05) > 0)
cncc_subtype_marker_genes <- unique(top50pval$gene)

human_mesenchyme <- subset(human_main_markers, cluster == "mesenchyme")
human_cncc <- subset(human_main_markers, cluster == "CNCC")
human_ectoderm <- subset(human_main_markers, cluster == "ectoderm")
human_muscle <- subset(human_main_markers, cluster == "muscle")
human_endothelium <- subset(human_main_markers, cluster == "endothelium")
human_blood <- subset(human_main_markers, cluster == "red blood cells")
human_immune <- subset(human_main_markers, cluster == "other blood cells")




emb.merge <- AddModuleScore(emb.merge, list(cnccMarkers), name = "CNCC_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(cnccMarkers_fulllist), name = "CNCC_Full_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(yankeeCNCC), name = "YankeeCNCC_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(treachercollinsgenes), name = "TreacherCollins_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(alx), name = "ALX_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(human_mesenchyme$gene), name = "Craniofacial_Mesenchyme_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(human_cncc$gene), name = "Craniofacial_CNCC_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(human_ectoderm$gene), name = "Craniofacial_Ectoderm_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(human_endothelium$gene), name = "Craniofacial_Endothelium_ModuleScore")
emb.merge <- AddModuleScore(emb.merge, list(human_muscle$gene), name = "Craniofacial_Muscle_ModuleScore")




for (n in craniofacial_subtypes){
emb.merge <- AddModuleScore(emb.merge, list(subtype_genelists[[n]]), name = gsub(" ",".",paste("Craniofacial",n,"ModuleScore", sep = "_")))
}

all_nsclp_genes <- readxl::read_xlsx("../modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 1)
nsclp_genes <- read.table("../modules/all_nsclp_genes.txt", header=TRUE, sep ="\t")
clp_genes <- read.table("../modules/clp_genes.txt", header=TRUE, sep ="\t")
all_clp_genes <- readxl::read_xlsx("../modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 2)
cp_genes <- read.table("../modules/cp_genes.txt", header=TRUE, sep ="\t")
all_cp_genes <- readxl::read_xlsx("../modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 3)
cleftgenedb <- read.table("../modules/cleftgenedb.txtr", header = TRUE)
cleftgenedb$Gene = toupper(cleftgenedb$Gene)
cleftgenedb <- cleftgenedb$Gene
novel_genes <- readxl::read_xlsx("../modules/Supplementary_Table_S7_update.xlsx", sheet = 2)
novel_genes <- novel_genes$SYMID
cfse_genes <- read.table("../modules/all_cfse_genes.txt", header = FALSE)
gnomad <- read.table("../gnomad.v4.1.constraint_metrics.tsv", header = TRUE)
gnomad_dec1 <- subset(gnomad, lof.oe_ci.upper_bin_decile == 1)
gnomad_dec9 <- subset(gnomad, lof.oe_ci.upper_bin_decile == 9)
gnomad_pli <- subset(gnomad, lof.pLI == 1)
gnomad_dec1_genes <- unique(gnomad_dec1$gene)
gnomad_dec9_genes <- unique(gnomad_dec9$gene)
gnomad_pli_genes <- unique(gnomad_pli$gene)

other_gene_lists <- list(cleftgenedb, novel_genes, gnomad_dec1_genes, gnomad_dec9_genes, gnomad_pli_genes, neanderthal_genes, all_har_genes)
names(other_gene_lists) <- c("cleftgenedb", "prioritized_genes", "gnomad_dec1", "gnomad_dec9", "gnomad_pli", "neanderthal", "har")

for (n in names(other_gene_lists))
  emb.merge <- AddModuleScore(emb.merge, list(other_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))


VlnPlot(emb.merge, "CNCC_ModuleScore1", sort = T)
VlnPlot(emb.merge, "CNCC_Full_ModuleScore1", sort = T)
VlnPlot(emb.merge, "YankeeCNCC_ModuleScore1", sort = T)
VlnPlot(emb.merge, "TreacherCollins_ModuleScore1", sort = T)
VlnPlot(emb.merge, "ALX_ModuleScore1", sort = T)
SpatialFeaturePlot(emb.merge, features = "ALX_ModuleScore1", image.alpha = 0.3) + theme(legend.position = "right")
SpatialFeaturePlot(emb.merge, features = "CNCC_ModuleScore1", image.alpha = 0.3) + theme(legend.position = "right")
SpatialFeaturePlot(emb.merge, features = "CNCC_Full_ModuleScore1", image.alpha = 0.3) + theme(legend.position = "right")
SpatialFeaturePlot(emb.merge, features = "Craniofacial_Mesenchyme_ModuleScore1", image.alpha = 0.3) + theme(legend.position = "right")
SpatialFeaturePlot(emb.merge, features = "Craniofacial_CNCC_ModuleScore1", image.alpha = 0.3) + theme(legend.position = "right")
SpatialFeaturePlot(emb.merge, features = "Craniofacial_Ectoderm_ModuleScore1", image.alpha = 0.3) + theme(legend.position = "right")
SpatialFeaturePlot(emb.merge, features = "Craniofacial_Muscle_ModuleScore1", image.alpha = 0.3) + theme(legend.position = "right")
SpatialFeaturePlot(emb.merge, features = "Craniofacial_Endothelium_ModuleScore1", image.alpha = 0.3) + theme(legend.position = "right")

pdf(file="whole_embryo_subtype_modulescores.pdf")
for (n in craniofacial_subtypes) {
  print(SpatialFeaturePlot(emb.merge, features = gsub(" ",".",paste("Craniofacial",n,"ModuleScore1", sep = "_")), image.alpha = 0.3) + theme(legend.position = "right"))
}
dev.off()


pdf(file="whole_embryo_other_lists_modulescores.pdf")
for (n in names(other_gene_lists)) {
  print(SpatialFeaturePlot(emb.merge, features = gsub(" ",".",paste(n,"ModuleScore1", sep = "_")), image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

pdf(file="whole_embryo_other_lists_violin_modulescores.pdf")
for (n in names(other_gene_lists)) {
  print(VlnPlot(emb.merge, gsub(" ",".",paste(n,"ModuleScore1", sep = "_"))) + theme(legend.position = "right")) 
}
dev.off()

pdf(file="whole_embryo_cleftgenedb.pdf")
for (n in sort(intersect(rownames(emb.merge),cleftgenedb))) {
  print(SpatialFeaturePlot(emb.merge, features = n, image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

pdf(file="plots/whole_embryo_mes_subtype_marker_genes.pdf")
for (n in sort(intersect(rownames(emb.merge),mes_subtype_marker_genes))) {
  print(SpatialFeaturePlot(emb.merge, features = n, image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

pdf(file="plots/whole_embryo_ect_subtype_marker_genes.pdf")
for (n in sort(intersect(rownames(emb.merge),ect_subtype_marker_genes))) {
  print(SpatialFeaturePlot(emb.merge, features = n, image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

pdf(file="plots/whole_embryo_cncc_subtype_marker_genes.pdf")
for (n in sort(intersect(rownames(emb.merge),cncc_subtype_marker_genes))) {
  print(SpatialFeaturePlot(emb.merge, features = n, image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

pdf(file="whole_embryo_prioritized_genes.pdf")
for (n in sort(intersect(rownames(emb.merge),novel_genes))) {
  print(SpatialFeaturePlot(emb.merge, features = n, image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

#subclustering of phyarngeal arches
emb.merge <- FindSubCluster(
      emb.merge,
      "SCT_snn",
      cluster = "phyarngeal arches",
      subcluster.name = "pa.sub.cluster",
      resolution = 1,
      algorithm = 1
  )
#assign labels to full data set
Idents(emb.merge) <- "pa.sub.cluster"

#subset only craniofacial spots and renormalize
craniofacial <- subset(emb.merge, idents = c("phyarngeal arches_0", "phyarngeal arches_1", "phyarngeal arches_2", "phyarngeal arches_3","craniofacial_ectoderm", "head mesoderm", "epidermis"))
craniofacial <- SCTransform(craniofacial, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
DimPlot(craniofacial, group.by = "ident", label = TRUE)


for (n in craniofacial_subtypes)
  craniofacial <- AddModuleScore(craniofacial, list(subtype_genelists[[n]]), name = gsub(" ",".",paste("Craniofacial",n,"ModuleScore", sep = "_")))

pdf(file="craniofacial_subtype_modulescores.pdf")
for (n in craniofacial_subtypes) {
  print(SpatialFeaturePlot(craniofacial, features = gsub(" ",".",paste("Craniofacial",n,"ModuleScore1", sep = "_")), image.alpha = 0.3) + theme(legend.position = "right"))
}
dev.off()


#integrate with single cell data
cds1 <- readRDS("../cds_face_human_temp.rds")
Idents(cds1) <- "tempanno"
DimPlot(cds1, group.by = "ident", label = TRUE)
cds1 <- SCTransform(cds1, verbose = TRUE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
DimPlot(cds1, group.by = "ident", label = TRUE)
#work in progress
#craniofacial integration
anchors <- FindTransferAnchors(reference = cds1, query = craniofacial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = cds1$tempanno, prediction.assay = TRUE,
                                  weight.reduction = craniofacial[["pca"]], dims = 1:30)
craniofacial[["predictions"]] <- predictions.assay
celltype_predictions <- rownames(craniofacial[["predictions"]]$data)


DefaultAssay(craniofacial) <- "predictions"

plot_spatial_predictions_craniofacial <- function(x){
  SpatialFeaturePlot(craniofacial, features = x, pt.size.factor = 1.6, ncol = 3, crop = TRUE)
}

pdf("craniofacial_spatial_predictions.pdf", w = 11, h = 8.5)
for (n in celltype_predictions){
  print(SpatialFeaturePlot(craniofacial, features = n, pt.size.factor = 1.6, crop = TRUE))
}
dev.off()

#full embryo integration
anchors <- FindTransferAnchors(reference = cds1, query = emb.merge, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = cds1$tempanno, prediction.assay = TRUE,
                                  weight.reduction = emb.merge[["pca"]], dims = 1:30)
emb.merge[["predictions"]] <- predictions.assay
celltype_predictions <- rownames(emb.merge[["predictions"]]$data)

plot_spatial_predictions_full_embryo <- function(x){
  SpatialFeaturePlot(emb.merge, features = x, pt.size.factor = 1.6, ncol = 3, crop = TRUE)
}


pdf("full_embryo_spatial_predictions.pdf", w = 11, h = 8.5)
lapply(celltype_predictions, plot_spatial_predictions_full_embryo)
dev.off()
