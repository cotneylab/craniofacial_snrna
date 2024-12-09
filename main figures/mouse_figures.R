library(ggplot2)
library(dplyr)
library(magrittr)
library(Seurat)
library(scCustomize)
library(tidyverse)
library(harmony)

mouse_cds1 <- readRDS("cds_face_mouse.rds")
mouse_cds2 <- mouse_cds1
mouse_cds_ect <- readRDS("ect_all.rds")
mouse_cds_mes <- readRDS("mes_all.rds")
Idents(mouse_cds1) <- "cell_type1"
mouse_cds_cncc <- readRDS("mcncc.rds")
Idents(mouse_cds2) <- "ectmes"
Idents(mouse_cds_ect) <- "cellType4"
Idents(mouse_cds_mes) <- "cellType4"
Idents(mouse_cds_cncc) <- "tempanno"

cluster_stats_main <- Cluster_Stats_All_Samples(seurat_object = mouse_cds1)
cluster_stats_subtypes <- Cluster_Stats_All_Samples(seurat_object = mouse_cds2)

require(biomaRt)
#need to use version 105 to get linked human gene symbols 
mart <-  useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 105)
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
annotLookup <- getLDS(
  mart = mart,
  attributes = c(
    "mgi_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filters = "mgi_symbol",
  attributesL = c("hgnc_symbol", "entrezgene_id"),
  martL = human,
  values = rownames(mouse_cds1),
  uniqueRows=TRUE)


features <- str_to_title(c("Prrx1","Alx1","Epcam","Ttn","Kel","Hemgn","Cdh5","Sox10","Fcer1g", "Irf6", "Pax7", "Sp7", "TCOF1", "TRP63", "TFAP2A"))
hox_genes <- str_to_title(c("HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10", "HOXA11", "HOXA13", "HOXB1", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXB9", "HOXB13", "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXC9", "HOXC10", "HOXC11", "HOXC12", "HOXC13", "HOXD1", "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HOXD10", "HOXD11", "HOXD12", "HOXD13"))
anterior_hox <- str_to_title(c("HOXA1", "HOXA2", "HOXA3", "HOXB1", "HOXB2", "HOXB3", "HOXD1", "HOXD3"))
central_hox <- str_to_title(c("HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXD4", "HOXD8"))
posterior_hox <-str_to_title(c("HOXA9", "HOXA10", "HOXA11", "HOXA13", "HOXB9", "HOXB13", "HOXC9", "HOXC10", "HOXC11", "HOXC12", "HOXC13", "HOXD9", "HOXD10", "HOXD11", "HOXD12", "HOXD13"))
alx <-str_to_title(c("ALX1", "ALX3", "ALX4"))


VlnPlot(mouse_cds1, features = features)
FeaturePlot(mouse_cds1, features = features)
DotPlot(mouse_cds1, features = features) + RotatedAxis()



p1 <- DimPlot(mouse_cds1, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(mouse_cds1, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(mouse_cds1,split.by = "stage",label = T)+NoAxes()+NoLegend()
p1+p2
p3

VlnPlot(mouse_cds1, features = "Xist", split.by = "orig.ident", group.by = "orig.ident") + RotatedAxis() + NoLegend()

sig_genes <- FindAllMarkers(mouse_cds1, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)

sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$MGI.symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID[match(as.character(sig_genes$gene),annotLookup$MGI.symbol)]
sig_genes$description <- annotLookup$Gene.description[match(as.character(sig_genes$gene),annotLookup$MGI.symbol)]
sig_genes$biotype <- annotLookup$Gene.type[match(as.character(sig_genes$gene),annotLookup$MGI.symbol)]
sig_genes$hgnc_symbol <- annotLookup$HGNC.symbol[match(as.character(sig_genes$gene),annotLookup$MGI.symbol)]
sig_genes$Human_Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID.1[match(as.character(sig_genes$gene),annotLookup$MGI.symbol)]

sig_genes <- sig_genes[, c("gene","Entrez","description","hgnc_symbol", "Human_Entrez", "biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
openxlsx::write.xlsx(sig_genes,"Mouse_SigGenes_Main_CellTypes.xlsx")

sig_genes_protein <- subset(sig_genes, biotype == "protein_coding")
sig_genes_protein <- subset(sig_genes_protein, Human_Entrez != "NA")
top <- sig_genes_protein %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top100 <- sig_genes_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top
top100pval <- subset(top100, rowSums(top100[7] < 0.05) > 0)

sig_genes_subtype <- FindAllMarkers(mouse_cds2, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                    min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                    only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                    latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                    pseudocount.use = 1, return.thresh = 0.01)

sig_genes_subtype <- mutate(sig_genes_subtype,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$MGI.symbol,sig_genes_subtype$gene)
x <- sig_genes_subtype[sig_genes_subtype$gene %in% id,]
sig_genes_subtype$Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID[match(as.character(sig_genes_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_subtype$description <- annotLookup$Gene.description[match(as.character(sig_genes_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_subtype$biotype <- annotLookup$Gene.type[match(as.character(sig_genes_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_subtype$hgnc_symbol <- annotLookup$HGNC.symbol[match(as.character(sig_genes_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_subtype$Human_Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID.1[match(as.character(sig_genes_subtype$gene),annotLookup$MGI.symbol)]

sig_genes_subtype <- sig_genes_subtype[, c("gene","Entrez","description","hgnc_symbol", "Human_Entrez", "biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_subtype)
openxlsx::write.xlsx(sig_genes_subtype,"Mouse_SigGenes_CellSubTypes.xlsx")

sig_genes_subtype_protein <- subset(sig_genes_subtype, biotype == "protein_coding")
sig_genes_subtype_protein <- subset(sig_genes_subtype_protein, Human_Entrez != "NA")

sig_genes_subtype_lncRNA <- subset(sig_genes_subtype, biotype == "lncRNA")



top_subtype <- sig_genes_subtype_protein %>% group_by(cluster) %>% top_n(5, avg_log2FC);top_subtype
top100_subtype <- sig_genes_subtype_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top_subtype
top100pval_subtype <- subset(top100_subtype, rowSums(top100_subtype[7] < 0.05) > 0)



sig_genes_cncc_subtype <- FindAllMarkers(mouse_cds_cncc, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                         min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                         only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                         latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                         pseudocount.use = 1, return.thresh = 0.01)

sig_genes_cncc_subtype <- mutate(sig_genes_cncc_subtype,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$MGI.symbol,sig_genes_cncc_subtype$gene)
x <- sig_genes_cncc_subtype[sig_genes_cncc_subtype$gene %in% id,]
sig_genes_cncc_subtype$Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID[match(as.character(sig_genes_cncc_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_cncc_subtype$description <- annotLookup$Gene.description[match(as.character(sig_genes_cncc_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_cncc_subtype$biotype <- annotLookup$Gene.type[match(as.character(sig_genes_cncc_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_cncc_subtype$hgnc_symbol <- annotLookup$HGNC.symbol[match(as.character(sig_genes_cncc_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_cncc_subtype$Human_Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID.1[match(as.character(sig_genes_cncc_subtype$gene),annotLookup$MGI.symbol)]

sig_genes_cncc_subtype <- sig_genes_cncc_subtype[, c("gene","Entrez","description","hgnc_symbol", "Human_Entrez", "biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_cncc_subtype)
openxlsx::write.xlsx(sig_genes_cncc_subtype,"tables/SigGenes_Cell_Mouse_CNCC_SubTypes.xlsx")

sig_genes_cncc_subtype_protein <- subset(sig_genes_cncc_subtype, biotype == "protein_coding")

sig_genes_cncc_subtype_lncRNA <- subset(sig_genes_cncc_subtype, biotype == "lncRNA")



top_cncc_subtype <- sig_genes_cncc_subtype %>% group_by(cluster) %>% top_n(5, avg_log2FC);top_cncc_subtype
top100_cncc_subtype <- sig_genes_cncc_subtype %>% group_by(cluster) %>% top_n(100, avg_log2FC);top100_cncc_subtype
top100pval_cncc_subtype <- subset(top100_cncc_subtype, rowSums(top100_cncc_subtype[7] < 1) > 0)

sig_genes_mes_subtype <- FindAllMarkers(mouse_cds_mes, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                        min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                        only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                        latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                        pseudocount.use = 1, return.thresh = 0.01)

sig_genes_mes_subtype <- mutate(sig_genes_mes_subtype,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$MGI.symbol,sig_genes_mes_subtype$gene)
x <- sig_genes_mes_subtype[sig_genes_mes_subtype$gene %in% id,]
sig_genes_mes_subtype$Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID[match(as.character(sig_genes_mes_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_mes_subtype$description <- annotLookup$Gene.description[match(as.character(sig_genes_mes_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_mes_subtype$biotype <- annotLookup$Gene.type[match(as.character(sig_genes_mes_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_mes_subtype$hgnc_symbol <- annotLookup$HGNC.symbol[match(as.character(sig_genes_mes_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_mes_subtype$Human_Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID.1[match(as.character(sig_genes_mes_subtype$gene),annotLookup$MGI.symbol)]

sig_genes_mes_subtype <- sig_genes_mes_subtype[, c("gene","Entrez","description","hgnc_symbol", "Human_Entrez", "biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_mes_subtype)
openxlsx::write.xlsx(sig_genes_mes_subtype,"tables/SigGenes_Cell_Mouse_mes_SubTypes.xlsx")

sig_genes_mes_subtype_protein <- subset(sig_genes_mes_subtype, biotype == "protein_coding")

sig_genes_mes_subtype_lncRNA <- subset(sig_genes_mes_subtype, biotype == "lncRNA")



top_mes_subtype <- sig_genes_mes_subtype_protein %>% group_by(cluster) %>% top_n(5, avg_log2FC);top_mes_subtype
top100_mes_subtype <- sig_genes_mes_subtype_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top100_mes_subtype
top100pval_mes_subtype <- subset(top100_mes_subtype, rowSums(top100_mes_subtype[7] < 0.05) > 0)

sig_genes_ect_subtype <- FindAllMarkers(mouse_cds_ect, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                        min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                        only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                        latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                        pseudocount.use = 1, return.thresh = 0.01)

sig_genes_ect_subtype <- mutate(sig_genes_ect_subtype,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$MGI.symbol,sig_genes_ect_subtype$gene)
x <- sig_genes_ect_subtype[sig_genes_ect_subtype$gene %in% id,]
sig_genes_ect_subtype$Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID[match(as.character(sig_genes_ect_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_ect_subtype$description <- annotLookup$Gene.description[match(as.character(sig_genes_ect_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_ect_subtype$biotype <- annotLookup$Gene.type[match(as.character(sig_genes_ect_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_ect_subtype$hgnc_symbol <- annotLookup$HGNC.symbol[match(as.character(sig_genes_ect_subtype$gene),annotLookup$MGI.symbol)]
sig_genes_ect_subtype$Human_Entrez <- annotLookup$NCBI.gene..formerly.Entrezgene..ID.1[match(as.character(sig_genes_ect_subtype$gene),annotLookup$MGI.symbol)]

sig_genes_ect_subtype <- sig_genes_ect_subtype[, c("gene","Entrez","description","hgnc_symbol", "Human_Entrez", "biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_ect_subtype)
openxlsx::write.xlsx(sig_genes_ect_subtype,"tables/SigGenes_Cell_Mouse_ect_SubTypes.xlsx")

sig_genes_ect_subtype_protein <- subset(sig_genes_ect_subtype, biotype == "protein_coding")

sig_genes_ect_subtype_lncRNA <- subset(sig_genes_ect_subtype, biotype == "lncRNA")



top_ect_subtype <- sig_genes_ect_subtype_protein %>% group_by(cluster) %>% top_n(5, avg_log2FC);top_ect_subtype
top100_ect_subtype <- sig_genes_ect_subtype_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top100_ect_subtype
top100pval_ect_subtype <- subset(top100_ect_subtype, rowSums(top100_ect_subtype[7] < 0.05) > 0)


p <- DotPlot(object = mouse_cds1, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
p1 <- DotPlot(object = mouse_cds2, features = unique(top_subtype$gene))+ coord_flip()+ RotatedAxis()
p2 <- DotPlot(object = mouse_cds_cncc, features = unique(top_cncc_subtype$gene))+ coord_flip()+ RotatedAxis()
p3 <- DotPlot(object = mouse_cds_mes, features = unique(top_mes_subtype$gene))+ coord_flip()+ RotatedAxis()
p4 <- DotPlot(object = mouse_cds_ect, features = unique(top_ect_subtype$gene))+ coord_flip()+ RotatedAxis()

df <- data.frame(Gene = p$data$features.plot, avg.exp = p$data$avg.exp.scaled, pct.exp = p$data$pct.exp, cluster = p$data$id)
df1 <- data.frame(Gene = p1$data$features.plot, avg.exp = p1$data$avg.exp.scaled, pct.exp = p1$data$pct.exp, cluster = p1$data$id)
df2 <- data.frame(Gene = p2$data$features.plot, avg.exp = p2$data$avg.exp.scaled, pct.exp = p2$data$pct.exp, cluster = p2$data$id)
df3 <- data.frame(Gene = p3$data$features.plot, avg.exp = p3$data$avg.exp.scaled, pct.exp = p3$data$pct.exp, cluster = p3$data$id)
df4 <- data.frame(Gene = p4$data$features.plot, avg.exp = p4$data$avg.exp.scaled, pct.exp = p4$data$pct.exp, cluster = p4$data$id)


p2 <- df %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(size = 6, face = "italic" )) +
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p2

p3 <- df1 %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  theme(axis.text.y = element_text(size = 2, face = "italic"))+
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p3

p4 <- df2 %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14)) +
  theme(axis.text.y = element_text(size = 12, face = "italic"))+
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p4

p5 <- df3 %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14)) +
  theme(axis.text.y = element_text(size = 8, face = "italic"))+
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p5

p6 <- df4 %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14)) +
  theme(axis.text.y = element_text(size = 8, face = "italic"))+
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p6
pdf(file = "plots/mouse_marker_gene_dotplots.pdf")
p2
p3
p4
p5
p6
dev.off()



#ClusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("AnnotationHub")
BiocManager::install("ReactomePA")

library("clusterProfiler")
library(enrichplot)
library("org.Hs.eg.db")
library("AnnotationHub")
library(ReactomePA)
library(DOSE)

#gene ontology enrichment across clusters
#need to use Entrez id in clusterprofiler
df2 <- top100pval[,c(5,12)]
df2sample <- split(df2$Human_Entrez,df2$cluster)
length(df2sample)
main_types <- names(df2sample)

df3 <- top100pval_subtype[,c(5,12)]
df3sample <- split(df3$Human_Entrez,df3$cluster)
length(df3sample)
subtypes <- names(df3sample)

df4 <- top100pval_subtype[,c(2,12)]
df4sample <- split(df4$Entrez,df3$cluster)
length(df4sample)
subtypes_mouse_entrez <- names(df4sample)

df5 <- top100pval_cncc_subtype[,c(2,12)]
df5 <- na.omit(df5)
df5sample <- split(df5$Entrez,df5$cluster)
length(df5sample)
subtypes_cncc <- names(df5sample)

df6 <- top100pval_mes_subtype[,c(2,12)]
df6 <- na.omit(df6)
df6sample <- split(df6$Entrez,df6$cluster)
length(df6sample)
subtypes_mes <- names(df6sample)

df7 <- top100pval_ect_subtype[,c(2,12)]
df7 <- na.omit(df7)
df7sample <- split(df7$Entrez,df7$cluster)
length(df7sample)
subtypes_ect <- names(df7sample)

genelist <- as.list(df2sample)

genelist_subtype <- as.list(df3sample)

genelist_subtype_mouse_entrez <- as.list(df4sample)

genelist_subtype_cncc <- as.list(df5sample)
genelist_subtype_mes <- as.list(df6sample)
genelist_subtype_ect <- as.list(df7sample)



BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "hsa")
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
pdf(file = "mouse_main_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
dev.off()

p1 <- emapplot(BPclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))



cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "mouse_main_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

options(enrichplot.colours = c("blue", "grey"))
p1 <- treeplot(BPclusterplot, label_format = 30, hclust_method = "average", nCluster = 7, showCategory = 10, geneClusterPanel = "dotplot") + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- treeplot(CCclusterplot, label_format = 30, hclust_method = "average", nCluster = 7, showCategory = 10, geneClusterPanel = "dotplot") + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- treeplot(MFclusterplot, label_format = 30, hclust_method = "average", nCluster = 7, showCategory = 10, geneClusterPanel = "dotplot") + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- treeplot(Pathwayclusterplot, label_format = 30, hclust_method = "average", nCluster = 7, showCategory = 10, geneClusterPanel = "dotplot") + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- treeplot(KEGGclusterplot, label_format = 30, hclust_method = "average", nCluster = 7, showCategory = 10, geneClusterPanel = "dotplot") + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- treeplot(DOclusterplot, label_format = 30, hclust_method = "average", nCluster = 7, showCategory = 10, geneClusterPanel = "dotplot") + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "mouse_main_cluster_enrichment_treeplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()


#per cluster enrichments
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

openxlsx::write.xlsx(go.ls, file = "mouse_gene_onology_bp_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

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

openxlsx::write.xlsx(go.ls, file = "mouse_gene_onology_cc_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

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

openxlsx::write.xlsx(go.ls, file = "mouse_gene_onology_mf_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "mouse_disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "mouse_kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE,
    organism = "human"
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "mouse_reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)



BPclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot_subtype <- pairwise_termsim(BPclusterplot_subtype)
CCclusterplot_subtype <- pairwise_termsim(CCclusterplot_subtype)
MFclusterplot_subtype <- pairwise_termsim(MFclusterplot_subtype)

Pathwayclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
Pathwayclusterplot_subtype <- pairwise_termsim(Pathwayclusterplot_subtype)

KEGGclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "hsa")
KEGGclusterplot_subtype <- pairwise_termsim(KEGGclusterplot_subtype)


DOclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot_subtype <- pairwise_termsim(DOclusterplot_subtype)


GOclusterplot_subtype_cncc <- enrichGO(
  gene          = genelist_subtype_mouse_entrez$CNCC,
  ont = 'ALL',
  OrgDb = "org.Mm.eg.db",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

GOclusterplot_subtype_cncc <- pairwise_termsim(GOclusterplot_subtype_cncc)
GOclusterplot_subtype_cncc2 <- simplify(GOclusterplot_subtype_cncc, cutoff=0.7, by="p.adjust", select_fun=min)

GOclusterplot_subtype_cncc_3 <- compareCluster(geneCluster = genelist_subtype_cncc, fun = "enrichGO", OrgDb = "org.Mm.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
GOclusterplot_subtype_cncc_3 <- pairwise_termsim(GOclusterplot_subtype_cncc_3)
GOclusterplot_subtype_cncc_4 <- simplify(GOclusterplot_subtype_cncc_3, cutoff=0.7, by="p.adjust", select_fun=min)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot_subtype, showCategory = 5, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot_subtype, showCategory = 5, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot_subtype, showCategory = 5, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot_subtype, showCategory = 5, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot_subtype, showCategory = 5, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot_subtype, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))


cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "mouse_subtype_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
dev.off()

options(enrichplot.colours = c("red", "grey"))
p7 <- dotplot(DOclusterplot_subtype_cncc, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("purple", "grey"))
p8 <- dotplot(GOclusterplot_subtype_cncc_4, showCategory = 1, font.size = 8, title = "GO", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=8),legend.position = "right",legend.box = "vertical",text = element_text(size=7), plot.title = element_text(hjust= 0.5, size = 8))

pdf(file = "plots/mouse_cncc_go_dotplot.pdf", h = 6.9, w = 4.8)
p8
dev.off()

p1 <- emapplot(BPclusterplot_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("MPO") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))



cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "mouse_subtype_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

options(enrichplot.colours = c("blue", "grey"))
p1 <- treeplot(BPclusterplot_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- treeplot(CCclusterplot_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- treeplot(MFclusterplot_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- treeplot(Pathwayclusterplot_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- treeplot(KEGGclusterplot_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- treeplot(DOclusterplot_subtype, label_format = 30, hclust_method = "average", nCluster = 40, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "mouse_subtype_cluster_enrichment_treeplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

go.ls <- genelist_subtype %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "mouse_gene_onology_bp_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_subtype %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "mouse_gene_onology_cc_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_subtype %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "mouse_gene_onology_mf_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist_subtype %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "mouse_disease_onology_subtype_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist_subtype %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "mouse_kegg_subtype_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist_subtype %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE,
    organism = "human"
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "mouse_reactome_subtype_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

cnccMarkers_fulllist <- str_to_title(c("Arpc1b", "Barx1", "Bmi1", "Cald1", "Cd44", "Cd47", "Cd9", "Cfl1", "Col1a1", "Crabp1", "Cyba", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Ebf1", "Epha2", "Ets1", "Foxd3", "Gli3", "Hk2", "Id1", "Id2", "Id3", "Id4", "Itgb1", "Klf4", "Lmo4", "Mef2c", "Msh6", "Msx1", "Myb", "Myc", "Myo10", "Pax3", "Pax7", "Prrx1", "Rac1", "Rarg", "Runx2", "Rxrg", "Snai1", "Snai2", "Sox10", "Sox5", "Sox6", "Sox9", "Spp1", "Tfap2a", "Twist1"))
#smaller list https://www.nature.com/articles/s42003-021-01970-0
cnccMarkers <- str_to_title(c("Ets1", "Foxd3", "Tfap2a", "Tfap2b", "NR2F1", "NR2F2", "Tcf3", "Tcf12"))
migratory_cnccMarkers <- str_to_title(c("Pax7", "Tfap2a"))
premigratory_cnccMarkers <- str_to_title(c("Sox9", "Sox10", "Pax3"))
