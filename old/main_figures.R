Sys.setenv('R_MAX_VSIZE'=64000000000)
library(Seurat)
library(tidyverse)
library(harmony)
library(dplyr)
library(patchwork)
library(viridis)
library(ggforce)
library(gghalves)
library(ggridges)
library(scCustomize)

initial <- readRDS("Initial_data.rds")
cds1 <- readRDS("~/Desktop/scRNA-Seq_GWAS/cds_face_human_temp.rds")
cds2 <- cds1
Idents(cds1) <- "cell_type1"

Idents(cds2) <- "tempanno"
cds_mes <- readRDS("mes_temp1.rds")
Idents(cds_mes) <- "MtoHAlltoALL"
cds_ect <- readRDS("ect_final.rds")
Idents(cds_ect) <- "tempanno"
cds_cncc <- readRDS("cncc_final.rds")
Idents(cds_cncc) <- "tempanno"

DefaultLayer(cds1[["RNA"]]) <- 'data'
DefaultLayer(cds2[["RNA"]]) <- 'data'
require(biomaRt)
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = rownames(cds1),
  uniqueRows=TRUE)



#barplots
ggplot(cds1@meta.data, aes(x=cell_type1, fill=stage)) +  geom_bar(position = "fill") + RotatedAxis()
pdf("plots/initial_cluster_barplot.pdf", h = 10, w = 10)
ggplot(initial@meta.data, aes(x=seurat_clusters, fill=orig.ident)) +  geom_bar(position = "fill") + RotatedAxis()
dev.off()

custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish)

custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-'))

table_samples_by_clusters <- cds1@meta.data %>%
  group_by(stage, cell_type1) %>%
  summarize(count = n()) %>%
  spread(cell_type1, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('stage', 'total_cell_count', everything())) %>%
  arrange(factor(stage, levels = levels(cds1@meta.data$stage)))

knitr::kable(table_samples_by_clusters)

table_clusters_by_samples <- cds1@meta.data %>%
  group_by(cell_type1, stage) %>%
  summarize(count = n()) %>%
  spread(stage, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('cell_type1', 'total_cell_count', everything())) %>%
  arrange(factor(cell_type1, levels = levels(cds1@meta.data$cell_type1)))

knitr::kable(table_clusters_by_samples)

temp_labels <- cds1@meta.data %>%
  group_by(stage) %>%
  tally()

p1 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'stage') %>%
  mutate(stage = factor(stage, levels = levels(cds1@meta.data$stage))) %>%
  ggplot(aes(stage, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = stage, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

temp_labels <- cds1@meta.data %>%
  group_by(cell_type1) %>%
  tally() %>%
  dplyr::rename('cell_type1' = cell_type1)

p2 <- table_clusters_by_samples %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cell_type1') %>%
  mutate(cell_type1 = factor(cell_type1, levels = levels(cds1@meta.data$cell_type1))) %>%
  ggplot(aes(cell_type1, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = cell_type1, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(
  'plots/composition_samples_clusters_by_number.pdf',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      cds1@meta.data$stage %>% unique() %>% length(),
      cds1@meta.data$cell_type1 %>% unique() %>% length()
    )),
  width = 18, height = 8
)


li_genes <- str_to_upper(c("Alx1", "Epcam", "Hemgn", "Cdh5", "Fcer1g"))
neuronal_genes <- c("TUBB3", "MAP2")
features <- str_to_upper(c("Prrx1","Alx1","Epcam","Ttn","Kel","Hemgn","Cdh5","Sox10","Fcer1g", "Irf6", "Pax7", "Sp7", "TCOF1", "TP63", "TFAP2A"))
hox_genes <- c("HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10", "HOXA11", "HOXA13", "HOXB1", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXB9", "HOXB13", "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXC9", "HOXC10", "HOXC11", "HOXC12", "HOXC13", "HOXD1", "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HOXD10", "HOXD11", "HOXD12", "HOXD13")
anterior_hox <- c("HOXA1", "HOXA2", "HOXA3", "HOXB1", "HOXB2", "HOXB3", "HOXD1", "HOXD3")
central_hox <- c("HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXD4", "HOXD8")
posterior_hox <-c("HOXA9", "HOXA10", "HOXA11", "HOXA13", "HOXB9", "HOXB13", "HOXC9", "HOXC10", "HOXC11", "HOXC12", "HOXC13", "HOXD9", "HOXD10", "HOXD11", "HOXD12", "HOXD13")
alx <-c("ALX1", "ALX3", "ALX4")


#results from denovolyzer anlaysis by Leslie lab
gmkf_results <- readxl::read_xlsx("leslie/20241118_geneintersection_and_enrichment.xlsx", sheet = 4)

df1 <- gmkf_results[,c(2,1)]
df1sample <- split(df1$ALL_OFCS,df1$Cluster)
length(df1sample)
gmkf_genes <- names(df1sample)
genelist_gmkf <- as.list(df1sample)

p1 <- VlnPlot(initial, features = c(li_genes,neuronal_genes), stack = T, flip = T)
p2 <- DimPlot(initial, label = T)+NoLegend()
p3 <- DimPlot(initial, group.by = "seurat_clusters")+NoAxes()
p4 <- DimPlot(initial,split.by = "stage",label = F)+NoAxes()+NoLegend()
p5 <- FeaturePlot(initial, features = li_genes, split.by = "stage") & NoAxes()
p6 <- FeaturePlot(initial, features = li_genes) & NoAxes()
p7 <- FeaturePlot(initial, features = neuronal_genes) & NoAxes()

pdf("plots/temp_umaps_neuronal_markers.pdf", h = 8.5, w = 11)
p1
p2
p3
p4
p5
p6
p7
dev.off()

p1 <- VlnPlot(cds1, features = li_genes, stack = T, flip = T)
p2 <- DimPlot(cds1, label = T)+NoLegend()
p3 <- DimPlot(cds1, group.by = "orig.ident")+NoAxes()
p4 <- DimPlot(cds1,split.by = "stage",label = F)+NoAxes()+NoLegend()
p5 <- FeaturePlot(cds1, features = li_genes, split.by = "stage") & NoAxes()
p6 <- FeaturePlot(cds1, features = li_genes) & NoAxes()

pdf("plots/initial_umaps_li_markers.pdf", h = 8.5, w = 11)
p1
p2
p3
p4
p5
p6
dev.off()
#check for M/F with XIST and other genes
male_genes = c("UTY", "DDX3Y", "KDM5D", "NLGN4Y")
female_genes = c("XIST", "TSIX")
sex_determination_genes = c(male_genes, female_genes)

VlnPlot(cds1, features = male_genes, split.by = "orig.ident", group.by = "orig.ident") + RotatedAxis() + NoLegend()
VlnPlot(cds1, features = female_genes, split.by = "orig.ident", group.by = "orig.ident") + RotatedAxis() + NoLegend()
VlnPlot(cds1, features = sex_determination_genes, split.by = "orig.ident", group.by = "orig.ident", stack = T, flip = T) & NoLegend()

pdf(file = "plots/xist_expression_per_sample.pdf")
VlnPlot(cds1, features = sex_determination_genes, split.by = "orig.ident", group.by = "orig.ident", stack = T, flip = T) & NoLegend()
dev.off()

sig_genes_initial <- FindAllMarkers(initial, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                    min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                    only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                    latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                    pseudocount.use = 1, return.thresh = 0.01)

sig_genes_initial <- mutate(sig_genes_initial,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$hgnc_symbol,sig_genes_initial$gene)
x <- sig_genes_initial[sig_genes_initial$gene %in% id,]
sig_genes_initial$Entrez <- annotLookup$entrezgene_id[match(as.character(sig_genes_initial$gene),annotLookup$hgnc_symbol)]
sig_genes_initial$description <- annotLookup$description[match(as.character(sig_genes_initial$gene),annotLookup$hgnc_symbol)]
sig_genes_initial$biotype <- annotLookup$gene_biotype[match(as.character(sig_genes_initial$gene),annotLookup$hgnc_symbol)]

sig_genes_initial <- sig_genes_initial[, c("gene","Entrez","description","biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_initial)
openxlsx::write.xlsx(sig_genes_initial,"tables/SigGenes_Initial_CellTypes.xlsx")

sig_genes_initial_protein <- subset(sig_genes_initial, biotype == "protein_coding")

top_initial <- sig_genes_initial_protein %>% group_by(cluster) %>% top_n(10, avg_log2FC);top_initial
top100_initial <- sig_genes_initial_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top_initial
top100pval_initial <- subset(top100_initial, rowSums(top100[5] < 0.05) > 0)

top50_intial <- sig_genes_initial_protein %>% group_by(cluster) %>% top_n(50, avg_log2FC);top_initial
top50pval_initial <- subset(top50_intial, rowSums(top50[5] < 0.05) > 0)


sig_genes <- FindAllMarkers(cds1, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
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
openxlsx::write.xlsx(sig_genes,"SigGenes_Main_CellTypes.xlsx")

sig_genes_protein <- subset(sig_genes, biotype == "protein_coding")

top <- sig_genes_protein %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top100 <- sig_genes_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

top50 <- sig_genes_protein %>% group_by(cluster) %>% top_n(50, avg_log2FC);top
top50pval <- subset(top50, rowSums(top50[5] < 0.05) > 0)

sig_genes_subtype <- FindAllMarkers(cds2, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                    min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                    only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                    latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                    pseudocount.use = 1, return.thresh = 0.01)

sig_genes_subtype <- mutate(sig_genes_subtype,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$hgnc_symbol,sig_genes_subtype$gene)
x <- sig_genes_subtype[sig_genes_subtype$gene %in% id,]
sig_genes_subtype$Entrez <- annotLookup$entrezgene_id[match(as.character(sig_genes_subtype$gene),annotLookup$hgnc_symbol)]
sig_genes_subtype$description <- annotLookup$description[match(as.character(sig_genes_subtype$gene),annotLookup$hgnc_symbol)]
sig_genes_subtype$biotype <- annotLookup$gene_biotype[match(as.character(sig_genes_subtype$gene),annotLookup$hgnc_symbol)]

sig_genes_subtype <- sig_genes_subtype[, c("gene","Entrez","description","biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_subtype)
openxlsx::write.xlsx(sig_genes_subtype,"SigGenes_CellSubTypes.xlsx")

sig_genes_subtype_protein <- subset(sig_genes_subtype, biotype == "protein_coding")

sig_genes_subtype_lncRNA <- subset(sig_genes_subtype, biotype == "lncRNA")



top_subtype <- sig_genes_subtype_protein %>% group_by(cluster) %>% top_n(5, avg_log2FC);top_subtype
top100_subtype <- sig_genes_subtype_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top_subtype
top100pval_subtype <- subset(top100_subtype, rowSums(top100_subtype[5] < 0.05) > 0)


sig_genes_cncc_subtype <- FindAllMarkers(cds_cncc, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                         min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                         only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                         latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                         pseudocount.use = 1, return.thresh = 0.01)

sig_genes_cncc_subtype <- mutate(sig_genes_cncc_subtype,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$hgnc_symbol,sig_genes_cncc_subtype$gene)
x <- sig_genes_cncc_subtype[sig_genes_cncc_subtype$gene %in% id,]
sig_genes_cncc_subtype$Entrez <- annotLookup$entrezgene_id[match(as.character(sig_genes_cncc_subtype$gene),annotLookup$hgnc_symbol)]
sig_genes_cncc_subtype$description <- annotLookup$description[match(as.character(sig_genes_cncc_subtype$gene),annotLookup$hgnc_symbol)]
sig_genes_cncc_subtype$biotype <- annotLookup$gene_biotype[match(as.character(sig_genes_cncc_subtype$gene),annotLookup$hgnc_symbol)]

sig_genes_cncc_subtype <- sig_genes_cncc_subtype[, c("gene","Entrez","description","biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_cncc_subtype)
openxlsx::write.xlsx(sig_genes_cncc_subtype,"SigGenes_Cell_CNCC_SubTypes.xlsx")

sig_genes_cncc_subtype_protein <- subset(sig_genes_cncc_subtype, biotype == "protein_coding")

sig_genes_cncc_subtype_lncRNA <- subset(sig_genes_cncc_subtype, biotype == "lncRNA")



top_cncc_subtype <- sig_genes_cncc_subtype_protein %>% group_by(cluster) %>% top_n(5, avg_log2FC);top_cncc_subtype
top100_cncc_subtype <- sig_genes_cncc_subtype_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top_cncc_subtype
top100pval_cncc_subtype <- subset(top100_cncc_subtype, rowSums(top100_cncc_subtype[5] < 0.05) > 0)

sig_genes_mes_subtype <- FindAllMarkers(cds_mes, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                        min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                        only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                        latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                        pseudocount.use = 1, return.thresh = 0.01)

sig_genes_mes_subtype <- mutate(sig_genes_mes_subtype,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$hgnc_symbol,sig_genes_mes_subtype$gene)
x <- sig_genes_mes_subtype[sig_genes_mes_subtype$gene %in% id,]
sig_genes_mes_subtype$Entrez <- annotLookup$entrezgene_id[match(as.character(sig_genes_mes_subtype$gene),annotLookup$hgnc_symbol)]
sig_genes_mes_subtype$description <- annotLookup$description[match(as.character(sig_genes_mes_subtype$gene),annotLookup$hgnc_symbol)]
sig_genes_mes_subtype$biotype <- annotLookup$gene_biotype[match(as.character(sig_genes_mes_subtype$gene),annotLookup$hgnc_symbol)]

sig_genes_mes_subtype <- sig_genes_mes_subtype[, c("gene","Entrez","description","biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_mes_subtype)
openxlsx::write.xlsx(sig_genes_mes_subtype,"SigGenes_Cell_mes_SubTypes.xlsx")

sig_genes_mes_subtype_protein <- subset(sig_genes_mes_subtype, biotype == "protein_coding")

sig_genes_mes_subtype_lncRNA <- subset(sig_genes_mes_subtype, biotype == "lncRNA")



top_mes_subtype <- sig_genes_mes_subtype_protein %>% group_by(cluster) %>% top_n(5, avg_log2FC);top_mes_subtype
top100_mes_subtype <- sig_genes_mes_subtype_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top_mes_subtype
top100pval_mes_subtype <- subset(top100_mes_subtype, rowSums(top100_mes_subtype[5] < 0.05) > 0)

sig_genes_ect_subtype <- FindAllMarkers(cds_ect, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                                        min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                                        only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                                        latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                                        pseudocount.use = 1, return.thresh = 0.01)

sig_genes_ect_subtype <- mutate(sig_genes_ect_subtype,pct.diff=pct.1-pct.2)
id <- intersect(annotLookup$hgnc_symbol,sig_genes_ect_subtype$gene)
x <- sig_genes_ect_subtype[sig_genes_ect_subtype$gene %in% id,]
sig_genes_ect_subtype$Entrez <- annotLookup$entrezgene_id[match(as.character(sig_genes_ect_subtype$gene),annotLookup$hgnc_symbol)]
sig_genes_ect_subtype$description <- annotLookup$description[match(as.character(sig_genes_ect_subtype$gene),annotLookup$hgnc_symbol)]
sig_genes_ect_subtype$biotype <- annotLookup$gene_biotype[match(as.character(sig_genes_ect_subtype$gene),annotLookup$hgnc_symbol)]

sig_genes_ect_subtype <- sig_genes_ect_subtype[, c("gene","Entrez","description","biotype","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes_ect_subtype)
openxlsx::write.xlsx(sig_genes_ect_subtype,"SigGenes_Cell_ect_SubTypes.xlsx")

sig_genes_ect_subtype_protein <- subset(sig_genes_ect_subtype, biotype == "protein_coding")

sig_genes_ect_subtype_lncRNA <- subset(sig_genes_ect_subtype, biotype == "lncRNA")



top_ect_subtype <- sig_genes_ect_subtype_protein %>% group_by(cluster) %>% top_n(5, avg_log2FC);top_ect_subtype
top100_ect_subtype <- sig_genes_ect_subtype_protein %>% group_by(cluster) %>% top_n(100, avg_log2FC);top_ect_subtype
top100pval_ect_subtype <- subset(top100_ect_subtype, rowSums(top100_ect_subtype[5] < 0.05) > 0)

p <- DotPlot(object = cds1, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
p1 <- DotPlot(object = cds2, features = unique(top_subtype$gene))+ coord_flip()+ RotatedAxis()
df <- data.frame(Gene = p$data$features.plot, avg.exp = p$data$avg.exp.scaled, pct.exp = p$data$pct.exp, cluster = p$data$id)
df1 <- data.frame(Gene = p1$data$features.plot, avg.exp = p1$data$avg.exp.scaled, pct.exp = p1$data$pct.exp, cluster = p1$data$id)

p4 <- DotPlot(object = initial, features = unique(top_initial$gene))+ coord_flip()+ RotatedAxis()
df4 <- data.frame(Gene = p4$data$features.plot, avg.exp = p4$data$avg.exp.scaled, pct.exp = p4$data$pct.exp, cluster = p4$data$id)

p2 <- df %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(size = 12, face = "italic" )) +
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

p5 <- df4 %>% 
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
p5

pdf(file = "plots/marker_gene_dotplots.pdf", h = 11, w = 8.5)
p2
p3
dev.off()

#ClusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
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
df1 <- top100pval_initial[,c(2,10)]
df1sample <- split(df1$Entrez,df1$cluster)
length(df1sample)
initial_types <- names(df1sample)

df2 <- top100pval[,c(2,10)]
df2sample <- split(df2$Entrez,df2$cluster)
length(df2sample)
main_types <- names(df2sample)

df2_symbol <- top50pval[,c(1,10)]
df2_symbol_sample <- split(df2_symbol$gene,df2_symbol$cluster)
length(df2_symbol_sample)
main_types_symbol <- names(df2_symbol_sample)

df3 <- top100pval_subtype[,c(2,10)]
df3sample <- split(df3$Entrez,df3$cluster)
length(df3sample)
subtypes <- names(df3sample)

df4 <- top100pval_cncc_subtype[,c(2,10)]
df4sample <- split(df4$Entrez,df4$cluster)
length(df4sample)
subtypes <- names(df4sample)

df5 <- top100pval_mes_subtype[,c(2,10)]
df5sample <- split(df5$Entrez,df5$cluster)
length(df5sample)
subtypes <- names(df5sample)

df6 <- top100pval_ect_subtype[,c(2,10)]
df6sample <- split(df6$Entrez,df6$cluster)
length(df6sample)
subtypes <- names(df6sample)

genelist_initial <- as.list(df1sample)
genelist <- as.list(df2sample)
genelist_symbol <- as.list(df2_symbol_sample)
genelist_subtype <- as.list(df3sample)

genelist_cncc_subtype <- as.list(df4sample)

genelist_mes_subtype <- as.list(df5sample)

genelist_ect_subtype <- as.list(df6sample)

background_genes <- unique(sig_genes_initial_protein$Entrez)

BPclusterplot <- compareCluster(geneCluster = genelist_initial, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP', universe = background_genes)
CCclusterplot <- compareCluster(geneCluster = genelist_initial, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC', universe = background_genes)
MFclusterplot <- compareCluster(geneCluster = genelist_initial, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF', universe = background_genes)
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)
options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 10, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 10, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 10, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))

cowplot::plot_grid(p1, p2, p3,ncol=1, labels=LETTERS[1:3])
pdf(file = "plots/initial_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
dev.off()



BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP', universe = background_genes)
CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC', universe = background_genes)
MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF', universe = background_genes)
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)


DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
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
pdf(file = "plots/main_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
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
pdf(file = "plots/main_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
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
pdf(file = "plots/main_cluster_enrichment_treeplot.pdf", h = 17, w = 22)
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

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_bp_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_cc_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_mf_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist %>% map(~{
 
  
 
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "tables/disease_onology_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "tables/kegg_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "tables/reactome_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

barplotTerm <- function(object,
                        x = "Count",
                        color = 'p.adjust',
                        showCategory = 8,
                        font.size = 12,
                        title = "") {

  colorBy <- color
  
  df <- fortify(object, showCategory = showCategory, by = x)
  df$p.adjust <- -log10(df$p.adjust)
  if (colorBy %in% colnames(df)) {
    p <-
      ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue",
        name = color,
        guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(x = x, y = "Description")) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  
  
  p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 p.adjust') + ylab(NULL)
}

pdf(file = "barplots_per_main_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist), function(x){
  
  name <- names(genelist)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()


BPclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP', universe = background_genes)
CCclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC',universe = background_genes)
MFclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF',universe = background_genes,)
BPclusterplot_subtype <- pairwise_termsim(BPclusterplot_subtype)
CCclusterplot_subtype <- pairwise_termsim(CCclusterplot_subtype)
MFclusterplot_subtype <- pairwise_termsim(MFclusterplot_subtype)

Pathwayclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH",universe = background_genes)
Pathwayclusterplot_subtype <- pairwise_termsim(Pathwayclusterplot_subtype)

KEGGclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH",universe = background_genes)
KEGGclusterplot_subtype <- pairwise_termsim(KEGGclusterplot_subtype)


DOclusterplot_subtype <- compareCluster(geneCluster = genelist_subtype, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH",universe = background_genes)
DOclusterplot_subtype <- pairwise_termsim(DOclusterplot_subtype)


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
pdf(file = "plots/subtype_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
dev.off()

pdf(file = "plots/subtype_cluster_disgenet_enrichment_dotplot.pdf", w = 17, h = 22)
p6
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
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))



cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
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
pdf(file = "plots/subtype_cluster_enrichment_treeplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()
pdf(file = "plots/subtype_cluster_disgenet_enrichment_treeplot.pdf", h = 17, w = 22)
p6
dev.off()
go.ls <- genelist_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    universe = background_genes,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_bp_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    universe = background_genes,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_cc_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_mf_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist_subtype %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    universe = background_genes,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "tables/disease_onology_subtype_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist_subtype %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "tables/kegg_subtype_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist_subtype %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "tables/reactome_subtype_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

barplotTerm <- function(object,
                        x = "Count",
                        color = 'p.adjust',
                        showCategory = 8,
                        font.size = 12,
                        title = "") {
  
  colorBy <- color
  
  df <- fortify(object, showCategory = showCategory, by = x)
  df$p.adjust <- -log10(df$p.adjust)
  if (colorBy %in% colnames(df)) {
    p <-
      ggplot(df, aes_string(x = x, y = "Description", fill = colorBy)) +
      theme_dose(font.size) +
      scale_fill_continuous(
        low = "red",
        high = "blue",
        name = color,
        guide = guide_colorbar(reverse = TRUE)
      )
  } else {
    p <- ggplot(df, aes_string(x = x, y = "Description")) +
      theme_dose(font.size) +
      theme(legend.position = "none")
  }
  
  
  p + geom_col(fill = color) + ggtitle(title) + xlab('-log10 p.adjust') + ylab(NULL)
}

pdf(file = "barplots_per_subtype_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist_subtype), function(x){
  
  name <- names(genelist_subtype)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()

BPclusterplot_cncc_subtype <- compareCluster(geneCluster = genelist_cncc_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP', universe = background_genes)
CCclusterplot_cncc_subtype <- compareCluster(geneCluster = genelist_cncc_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC', universe = background_genes)
MFclusterplot_cncc_subtype <- compareCluster(geneCluster = genelist_cncc_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF', universe = background_genes)
BPclusterplot_cncc_subtype <- pairwise_termsim(BPclusterplot_cncc_subtype)
CCclusterplot_cncc_subtype <- pairwise_termsim(CCclusterplot_cncc_subtype)
MFclusterplot_cncc_subtype <- pairwise_termsim(MFclusterplot_cncc_subtype)

Pathwayclusterplot_cncc_subtype <- compareCluster(geneCluster = genelist_cncc_subtype, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
Pathwayclusterplot_cncc_subtype <- pairwise_termsim(Pathwayclusterplot_cncc_subtype)

KEGGclusterplot_cncc_subtype <- compareCluster(geneCluster = genelist_cncc_subtype, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
KEGGclusterplot_cncc_subtype <- pairwise_termsim(KEGGclusterplot_cncc_subtype)


DOclusterplot_cncc_subtype <- compareCluster(geneCluster = genelist_cncc_subtype, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
DOclusterplot_cncc_subtype <- pairwise_termsim(DOclusterplot_cncc_subtype)


options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot_cncc_subtype, showCategory = 5, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot_cncc_subtype, showCategory = 5, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot_cncc_subtype, showCategory = 5, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot_cncc_subtype, showCategory = 5, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot_cncc_subtype, showCategory = 5, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot_cncc_subtype, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))


cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_cncc_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
dev.off()

pdf(file = "plots/subtype_cncc_cluster_disgenet_enrichment_dotplot.pdf", w = 17, h = 22)
p6
dev.off()

p1 <- emapplot(BPclusterplot_cncc_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot_cncc_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot_cncc_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot_cncc_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot_cncc_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot_cncc_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))



cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_cncc_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

options(enrichplot.colours = c("blue", "grey"))
p1 <- treeplot(BPclusterplot_cncc_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- treeplot(CCclusterplot_cncc_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- treeplot(MFclusterplot_cncc_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- treeplot(Pathwayclusterplot_cncc_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- treeplot(KEGGclusterplot_cncc_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- treeplot(DOclusterplot_cncc_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_cncc_cluster_enrichment_treeplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()
pdf(file = "plots/subtype_cncc_cluster_disgenet_enrichment_treeplot.pdf", h = 17, w = 22)
p6
dev.off()
go.ls <- genelist_cncc_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_bp_cncc_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_cncc_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_cc_cncc_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_cncc_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_mf_cncc_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist_cncc_subtype %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    universe = background_genes,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "tables/disease_onology_cncc_subtype_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist_cncc_subtype %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "tables/kegg_cncc_subtype_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist_cncc_subtype %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "tables/reactome_cncc_subtype_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

pdf(file = "barplots_per_cncc_subtype_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist_cncc_subtype), function(x){
  
  name <- names(genelist_cncc_subtype)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()


BPclusterplot_mes_subtype <- compareCluster(geneCluster = genelist_mes_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP', universe = background_genes)
CCclusterplot_mes_subtype <- compareCluster(geneCluster = genelist_mes_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC', universe = background_genes)
MFclusterplot_mes_subtype <- compareCluster(geneCluster = genelist_mes_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF', universe = background_genes)
BPclusterplot_mes_subtype <- pairwise_termsim(BPclusterplot_mes_subtype)
BPclusterplot_mes_subtype_simplify <- clusterProfiler::simplify(x = BPclusterplot_mes_subtype, cutoff = 0.5)
CCclusterplot_mes_subtype <- pairwise_termsim(CCclusterplot_mes_subtype)
MFclusterplot_mes_subtype <- pairwise_termsim(MFclusterplot_mes_subtype)

Pathwayclusterplot_mes_subtype <- compareCluster(geneCluster = genelist_mes_subtype, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
Pathwayclusterplot_mes_subtype <- pairwise_termsim(Pathwayclusterplot_mes_subtype)

KEGGclusterplot_mes_subtype <- compareCluster(geneCluster = genelist_mes_subtype, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
KEGGclusterplot_mes_subtype <- pairwise_termsim(KEGGclusterplot_mes_subtype)


DOclusterplot_mes_subtype <- compareCluster(geneCluster = genelist_mes_subtype, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
DOclusterplot_mes_subtype <- pairwise_termsim(DOclusterplot_mes_subtype)


options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot_mes_subtype, showCategory = 5, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot_mes_subtype, showCategory = 5, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot_mes_subtype, showCategory = 5, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot_mes_subtype, showCategory = 5, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot_mes_subtype, showCategory = 5, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot_mes_subtype, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))


cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_mes_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
dev.off()


options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot_mes_subtype, showCategory = 5, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "bottom",legend.box = "horizontal",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
pdf(file = "plots/subtype_mes_cluster_go_bp_enrichment_dotplot.pdf", w = 6, h = 10)
p1
dev.off()

options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot_mes_subtype, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=6),legend.position = "bottom",legend.box = "horizontal",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))


pdf(file = "plots/subtype_mes_cluster_disgenet_enrichment_dotplot.pdf", w = 6, h = 10)
p6
dev.off()

options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot_mes_subtype, showCategory = c("Craniofacial Abnormalities","Broad forehead","Branchio-Oto-Renal Syndrome","Craniofacial Dysostosis","Cranioschisis","Hypoplasia of the maxilla","Unilateral cleft lip", "Cleft upper lip", "Lip pit","Cleft palate, isolated", "Depressed nasal bridge", "Cleft lip or lips", "Retrognathia", "Cleft Lip with or without Cleft Palate"), font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=6),legend.position = "bottom",legend.box = "horizontal",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))


p1 <- emapplot(BPclusterplot_mes_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot_mes_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot_mes_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot_mes_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot_mes_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot_mes_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))



cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_mes_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

options(enrichplot.colours = c("blue", "grey"))
p1 <- treeplot(BPclusterplot_mes_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- treeplot(CCclusterplot_mes_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- treeplot(MFclusterplot_mes_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- treeplot(Pathwayclusterplot_mes_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- treeplot(KEGGclusterplot_mes_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- treeplot(DOclusterplot_mes_subtype, label_format = 30, hclust_method = "average", nCluster = 40, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_mes_cluster_enrichment_treeplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()
pdf(file = "plots/subtype_mes_cluster_disgenet_enrichment_treeplot.pdf", h = 17, w = 22)
p6
dev.off()
go.ls <- genelist_mes_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_bp_mes_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_mes_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_cc_mes_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_mes_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_mf_mes_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist_mes_subtype %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    universe = background_genes,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "tables/disease_onology_mes_subtype_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist_mes_subtype %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "tables/kegg_mes_subtype_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist_mes_subtype %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "tables/reactome_mes_subtype_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

pdf(file = "barplots_per_mes_subtype_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist_mes_subtype), function(x){
  
  name <- names(genelist_mes_subtype)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()

BPclusterplot_ect_subtype <- compareCluster(geneCluster = genelist_ect_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP',universe = background_genes)
CCclusterplot_ect_subtype <- compareCluster(geneCluster = genelist_ect_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC',universe = background_genes)
MFclusterplot_ect_subtype <- compareCluster(geneCluster = genelist_ect_subtype, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF',universe = background_genes)
BPclusterplot_ect_subtype <- pairwise_termsim(BPclusterplot_ect_subtype)
CCclusterplot_ect_subtype <- pairwise_termsim(CCclusterplot_ect_subtype)
MFclusterplot_ect_subtype <- pairwise_termsim(MFclusterplot_ect_subtype)

Pathwayclusterplot_ect_subtype <- compareCluster(geneCluster = genelist_ect_subtype, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH",universe = background_genes)
Pathwayclusterplot_ect_subtype <- pairwise_termsim(Pathwayclusterplot_ect_subtype)

KEGGclusterplot_ect_subtype <- compareCluster(geneCluster = genelist_ect_subtype, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH",universe = background_genes)
KEGGclusterplot_ect_subtype <- pairwise_termsim(KEGGclusterplot_ect_subtype)


DOclusterplot_ect_subtype <- compareCluster(geneCluster = genelist_ect_subtype, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH",universe = background_genes)
DOclusterplot_ect_subtype <- pairwise_termsim(DOclusterplot_ect_subtype)


options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot_ect_subtype, showCategory = 5, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot_ect_subtype, showCategory = 5, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot_ect_subtype, showCategory = 5, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot_ect_subtype, showCategory = 5, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot_ect_subtype, showCategory = 5, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot_ect_subtype, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))


cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_ect_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])
dev.off()

options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot_ect_subtype, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=6),legend.position = "bottom",legend.box = "horizontal",text = element_text(size=6), axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))


pdf(file = "plots/subtype_ect_cluster_disgenet_enrichment_dotplot.pdf", w = 6, h = 10)
p6
dev.off()

p1 <- emapplot(BPclusterplot_ect_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p1 <- p1 + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p2 <- emapplot(CCclusterplot_ect_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p2 <- p2 + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p3 <- emapplot(MFclusterplot_ect_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p3 <- p3 + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p4 <- emapplot(Pathwayclusterplot_ect_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p4 <- p4 + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p5 <- emapplot(KEGGclusterplot_ect_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p5 <- p5 + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))
p6 <- emapplot(DOclusterplot_ect_subtype, cex_line = 0.1, cex_label_category = 0.6, cex_category = 5, pie="count", repel = TRUE, min_edge = 0.25, legend_n=3, layout="kk")
p6 <- p6 + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", plot.title = element_text(hjust= 0.5, size = 12))



cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_ect_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

options(enrichplot.colours = c("blue", "grey"))
p1 <- treeplot(BPclusterplot_ect_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Biological Process") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- treeplot(CCclusterplot_ect_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Cellular Component") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- treeplot(MFclusterplot_ect_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("GO: Molecular Function") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- treeplot(Pathwayclusterplot_ect_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("Reactome") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- treeplot(KEGGclusterplot_ect_subtype, label_format = 30, hclust_method = "average", nCluster = 20, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("KEGG") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- treeplot(DOclusterplot_ect_subtype, label_format = 30, hclust_method = "average", nCluster = 40, showCategory = 5, geneClusterPanel = "dotplot") + ggtitle("DisGeNet") + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical", axis.text.x=element_text(angle=45, hjust = 1),plot.title = element_text(hjust= 0.5, size = 12))
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
pdf(file = "plots/subtype_ect_cluster_enrichment_treeplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()
pdf(file = "plots/subtype_ect_cluster_disgenet_enrichment_treeplot.pdf", h = 17, w = 22)
p6
dev.off()
go.ls <- genelist_ect_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_bp_ect_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_ect_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_cc_ect_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- genelist_ect_subtype %>% map(~{
  
  
  
  eGO <- enrichGO(
    gene          = .x,
    OrgDb         = org.Hs.eg.db,
    ont           = "MF",
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eGO)
  
})

openxlsx::write.xlsx(go.ls, file = "tables/gene_onology_mf_ect_subtype_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- genelist_ect_subtype %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    universe = background_genes,
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "tables/disease_onology_ect_subtype_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- genelist_ect_subtype %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    universe = background_genes,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "tables/kegg_ect_subtype_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- genelist_ect_subtype %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    universe = background_genes,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "tables/reactome_ect_subtype_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

pdf(file = "barplots_per_ect_subtype_cluster.pdf", h = 11 , w = 8)
lapply(1:length(genelist_ect_subtype), function(x){
  
  name <- names(genelist_ect_subtype)[[x]]
  g1 = barplotTerm(go.ls[[x]], showCategory = 25, title = paste0(name, " GO"), color = 'blue', x = 'p.adjust')
  g3 = barplotTerm(do.ls[[x]], showCategory = 25, title = paste0(name, " DisGeNet"), color = 'red', x = 'p.adjust')
  print(g1)
  print(g3)
  
})
dev.off()

# CNCC module score doi: 10.1242/dev.105445
# https://www.sigmaaldrich.com/US/en/product/mm/scc049
# https://www.liebertpub.com/doi/pdf/10.1089/scd.2012.0155
# full list
cnccMarkers_fulllist <- str_to_upper(c("Arpc1b", "Barx1", "Bmi1", "Cald1", "Cd44", "Cd47", "Cd9", "Cfl1", "Col1a1", "Crabp1", "Cyba", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Ebf1", "Epha2", "Ets1", "Foxd3", "Gli3", "Hk2", "Id1", "Id2", "Id3", "Id4", "Itgb1", "Klf4", "Lmo4", "Mef2c", "Msh6", "Msx1", "Myb", "Myc", "Myo10", "Pax3", "Pax7", "Prrx1", "Rac1", "Rarg", "Runx2", "Rxrg", "Snai1", "Snai2", "Sox10", "Sox5", "Sox6", "Sox9", "Spp1", "Tfap2a", "Twist1"))
#smaller list https://www.nature.com/articles/s42003-021-01970-0
cnccMarkers <- str_to_upper(c("Ets1", "Foxd3", "Tfap2a", "Tfap2b", "NR2F1", "NR2F2", "Tcf3", "Tcf12"))
migratory_cnccMarkers <- str_to_upper(c("Pax7", "Tfap2a"))
premigratory_cnccMarkers <- str_to_upper(c("Sox9", "Sox10", "Pax3"))

dark_green <- read.table(file = "yankee_de_clusters.txt", header = TRUE, row.names = 1, sep = "\t")
dark_green <- subset(dark_green, Cluster == "Dark_Green")
dark_green <- subset(dark_green, SYM != "NA")
dark_green_genes <- dark_green$SYM


hescvcnccdegenes <- read.table(file ="res_h9ncc.txt", header = TRUE, row.names = 1, sep = "\t")
yankeecnccMarkers <- subset(hescvcnccdegenes, padj < 0.05)
yankeecnccMarkers <- subset(yankeecnccMarkers, log2FoldChange > 4)
yankeecnccMarkers <- subset(yankeecnccMarkers, Symbol != "NA")
yankeeCNCC <- yankeecnccMarkers$Symbol

treachercollinsgenes <- c("TCOF1", "POLR1C", "POLR1D")

human_main_markers <- readxl::read_xlsx("SigGenes_Main_CellTypes.xlsx", sheet = 1)
human_main_markers <- subset(human_main_markers, biotype == "protein_coding")
top <- human_main_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top100 <- human_main_markers %>% group_by(cluster) %>% top_n(100, avg_log2FC);top
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)
main_celltypes <- unique(top100pval$cluster)

extract_subtype_genelists <- function(x) {
  main_genelists[[x]] <- subset(top100pval, cluster == x)$gene
}

main_genelists <- list()

main_genelists <- lapply(main_celltypes, extract_genelists)
names(main_genelists) <- main_celltypes

human_subtype_markers <- readxl::read_xlsx("SigGenes_CellSubTypes.xlsx", sheet = 1)
human_subtype_markers <- subset(human_subtype_markers, biotype == "protein_coding")

top <- human_subtype_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC);top
top100 <- human_subtype_markers %>% group_by(cluster) %>% top_n(100, avg_log2FC);top
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)
craniofacial_subtypes <- unique(top100pval$cluster)

extract_subtype_genelists <- function(x) {
  subtype_genelists[[x]] <- subset(top100pval, cluster == x)$gene
}

subtype_genelists <- list()

subtype_genelists <- lapply(craniofacial_subtypes, extract_subtype_genelists)
names(subtype_genelists) <- craniofacial_subtypes


filenames <- list.files(path = "modules",pattern="*_Nodes.csv", full.names = TRUE)

modules <- lapply(setNames(filenames, make.names(gsub("*_Nodes.csv", "", filenames))),read.csv)
modules_symbols <- lapply(modules, function(x){x[,2]})
hubs  <- lapply(modules, function(x) { x[ x$HUB == 1, 2] })
black_hubs <- hubs$modules.black
black_genes <- modules$modules.black$SYMID

all_nsclp_genes <- readxl::read_xlsx("modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 1)
nsclp_genes <- read.table("modules/all_nsclp_genes.txt", header=TRUE, sep ="\t")
clp_genes <- read.table("modules/clp_genes.txt", header=TRUE, sep ="\t")
all_clp_genes <- readxl::read_xlsx("modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 2)
cp_genes <- read.table("modules/cp_genes.txt", header=TRUE, sep ="\t")
all_cp_genes <- readxl::read_xlsx("modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 3)
cleftgenedb <- read.table("modules/cleftgenedb.txtr", header = TRUE)
cleftgenedb$Gene = toupper(cleftgenedb$Gene)
cleftgenedb <- cleftgenedb$Gene
novel_genes <- readxl::read_xlsx("modules/Supplementary_Table_S7_update.xlsx", sheet = 2)
novel_genes <- novel_genes$SYMID
cfse_genes <- read.table("modules/all_cfse_genes.txt", header = FALSE)
gnomad <- read.table("gnomad.v4.1.constraint_metrics.tsv", header = TRUE)
gnomad_dec1 <- subset(gnomad, lof.oe_ci.upper_bin_decile == 1)
gnomad_dec9 <- subset(gnomad, lof.oe_ci.upper_bin_decile == 9)
gnomad_pli <- subset(gnomad, lof.pLI == 1)
gnomad_dec1_genes <- unique(gnomad_dec1$gene)
gnomad_dec9_genes <- unique(gnomad_dec9$gene)
gnomad_pli_genes <- unique(gnomad_pli$gene)

library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

library("rGREAT")

dev_bed_file <- "ALLEUR_sprime_segments_neanMatchingFilter.bed"
bed.gr <- import(dev_bed_file, genome = "hg19")

job2 = submitGreatJob(bed.gr,
                      genome = "hg19",
                      rule = c("oneClosest"))

tbl = getEnrichmentTables(job2)
gene_regions = getRegionGeneAssociations(job2, verbose = great_opt$verbose)
neanderthal_genes <- unique(as.data.frame(gene_regions$annotated_genes)$value)

#table s1 from keough et al science 2023
zoohar_bed_file <- "zooHAR.hg38.bed"
har.gr <- import(zoohar_bed_file, genome = "hg38")

job3 = submitGreatJob(har.gr,
                      genome = "hg38",
                      rule = c("oneClosest"))
har.tbl = getEnrichmentTables(job3)
har_gene_regions = getRegionGeneAssociations(job3, verbose = great_opt$verbose)
all_har_genes <- unique(as.data.frame(har_gene_regions$annotated_genes)$value)

other_gene_lists <- list(cleftgenedb, novel_genes, gnomad_dec1_genes, gnomad_dec9_genes, gnomad_pli_genes, neanderthal_genes, all_har_genes)
names(other_gene_lists) <- c("cleftgenedb", "prioritized_genes", "gnomad_dec1", "gnomad_dec9", "gnomad_pli", "neanderthal", "har")

full_face_1 <- "Hostens_2021/hostens_full_face_1_hg19.bed"
bed.full_face_1 <- import(full_face_1, genome = "hg19")
job_full_face_1 = submitGreatJob(bed.full_face_1,
                                 genome = "hg19",
                                 rule = c("oneClosest"))
full_face_1.tbl = getEnrichmentTables(job_full_face_1)
full_face_1_gene_regions = getRegionGeneAssociations(job_full_face_1, verbose = great_opt$verbose)
full_face_1_genes <- unique(as.data.frame(full_face_1_gene_regions$annotated_genes)$value)


lower_face_III <- "Hostens_2021/hostens_lower_face_III_hg19.bed"
bed.lower_face_III <- import(lower_face_III, genome = "hg19")
job_lower_face_III = submitGreatJob(bed.lower_face_III,
                                    genome = "hg19",
                                    rule = c("oneClosest"))
lower_face_III.tbl = getEnrichmentTables(job_lower_face_III)
lower_face_III_gene_regions = getRegionGeneAssociations(job_lower_face_III, verbose = great_opt$verbose)
lower_face_III_genes <- unique(as.data.frame(lower_face_III_gene_regions$annotated_genes)$value)


midface_2 <- "Hostens_2021/hostens_midface_2_hg19.bed"
bed.midface_2 <- import(midface_2, genome = "hg19")
job_midface_2 = submitGreatJob(bed.midface_2,
                               genome = "hg19",
                               rule = c("oneClosest"))
midface_2.tbl = getEnrichmentTables(job_midface_2)
midface_2_gene_regions = getRegionGeneAssociations(job_midface_2, verbose = great_opt$verbose)
midface_2_genes <- unique(as.data.frame(midface_2_gene_regions$annotated_genes)$value)


nose_II <- "Hostens_2021/hostens_nose_II_hg19.bed"
bed.nose_II <- import(nose_II, genome = "hg19")
job_nose_II = submitGreatJob(bed.nose_II,
                             genome = "hg19",
                             rule = c("oneClosest"))
nose_II.tbl = getEnrichmentTables(job_nose_II)
nose_II_gene_regions = getRegionGeneAssociations(job_nose_II, verbose = great_opt$verbose)
nose_II_genes <- unique(as.data.frame(nose_II_gene_regions$annotated_genes)$value)


outer_face_3 <- "Hostens_2021/hostens_outer_face_3_hg19.bed"
bed.outer_face_3 <- import(outer_face_3, genome = "hg19")
job_outer_face_3 = submitGreatJob(bed.outer_face_3,
                                  genome = "hg19",
                                  rule = c("oneClosest"))
outer_face_3.tbl = getEnrichmentTables(job_outer_face_3)
outer_face_3_gene_regions = getRegionGeneAssociations(job_outer_face_3, verbose = great_opt$verbose)
outer_face_3_genes <- unique(as.data.frame(outer_face_3_gene_regions$annotated_genes)$value)


philtrum_I <- "Hostens_2021/hostens_philtrum_I_hg19.bed"
bed.philtrum_I <- import(philtrum_I, genome = "hg19")
job_philtrum_I = submitGreatJob(bed.philtrum_I,
                                genome = "hg19",
                                rule = c("oneClosest"))
philtrum_I.tbl = getEnrichmentTables(job_philtrum_I)
philtrum_I_gene_regions = getRegionGeneAssociations(job_philtrum_I, verbose = great_opt$verbose)
philtrum_I_genes <- unique(as.data.frame(philtrum_I_gene_regions$annotated_genes)$value)


upper_face_IV <- "Hostens_2021/hostens_upper_face_IV_hg19.bed"
bed.upper_face_IV <- import(upper_face_IV, genome = "hg19")
job_upper_face_IV = submitGreatJob(bed.upper_face_IV,
                                   genome = "hg19",
                                   rule = c("oneClosest"))
upper_face_IV.tbl = getEnrichmentTables(job_upper_face_IV)
upper_face_IV_gene_regions = getRegionGeneAssociations(job_upper_face_IV, verbose = great_opt$verbose)
upper_face_IV_genes <- unique(as.data.frame(upper_face_IV_gene_regions$annotated_genes)$value)

facial_variation_gene_lists <- list(full_face_1_genes, lower_face_III_genes, midface_2_genes, nose_II_genes, outer_face_3_genes,philtrum_I_genes,upper_face_IV_genes)
names(facial_variation_gene_lists) <- c("full_face_1", "lower_face_III", "midface_2", "nose_II", "outer_face_3","philtrum_I","upper_face_IV")

for (n in names(facial_variation_gene_lists))
  cds1 <- AddModuleScore(cds1, list(facial_variation_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))

for (n in names(other_gene_lists))
  cds1 <- AddModuleScore(cds1, list(other_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))


cds1 <- AddModuleScore(cds1, list(cnccMarkers), name = "CNCC_ModuleScore")
cds1 <- AddModuleScore(cds1, list(cnccMarkers_fulllist), name = "CNCC_Full_ModuleScore")
cds1 <- AddModuleScore(cds1, list(yankeeCNCC), name = "YankeeCNCC_ModuleScore")
cds1 <- AddModuleScore(cds1, list(dark_green_genes), name = "darkGreen_ModuleScore")
cds1 <- AddModuleScore(cds1, list(treachercollinsgenes), name = "TreacherCollins_ModuleScore")
cds1 <- AddModuleScore(cds1, list(alx), name = "ALX_ModuleScore")


for (n in names(facial_variation_gene_lists))
  cds2 <- AddModuleScore(cds2, list(facial_variation_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))

for (n in names(other_gene_lists))
  cds2 <- AddModuleScore(cds2, list(other_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))

cds2 <- AddModuleScore(cds2, list(cnccMarkers), name = "CNCC_ModuleScore")
cds2 <- AddModuleScore(cds2, list(anterior_hox), name = "Anterior_ModuleScore")
cds2 <- AddModuleScore(cds2, list(central_hox), name = "Central_ModuleScore")
cds2 <- AddModuleScore(cds2, list(posterior_hox), name = "Posterior_ModuleScore")
cds2 <- AddModuleScore(cds2, list(treachercollinsgenes), name = "TreacherCollins_ModuleScore")
cds2 <- AddModuleScore(cds2, list(alx), name = "ALX_ModuleScore")

for (n in names(facial_variation_gene_lists))
  cds_mes <- AddModuleScore(cds_mes, list(facial_variation_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))

for (n in names(other_gene_lists))
  cds_mes <- AddModuleScore(cds_mes, list(other_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))

cds_mes <- AddModuleScore(cds_mes, list(anterior_hox), name = "Anterior_ModuleScore")
cds_mes <- AddModuleScore(cds_mes, list(central_hox), name = "Central_ModuleScore")
cds_mes <- AddModuleScore(cds_mes, list(posterior_hox), name = "Posterior_ModuleScore")
cds_mes <- AddModuleScore(cds_mes, list(treachercollinsgenes), name = "TreacherCollins_ModuleScore")
cds_mes <- AddModuleScore(cds_mes, list(alx), name = "ALX_ModuleScore")

for (n in names(facial_variation_gene_lists))
  cds_ect <- AddModuleScore(cds_mes, list(facial_variation_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))

for (n in names(other_gene_lists))
  cds_ect <- AddModuleScore(cds_mes, list(other_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))

cds_ect <- AddModuleScore(cds_ect, list(anterior_hox), name = "Anterior_ModuleScore")
cds_ect <- AddModuleScore(cds_ect, list(central_hox), name = "Central_ModuleScore")
cds_ect <- AddModuleScore(cds_ect, list(posterior_hox), name = "Posterior_ModuleScore")
cds_ect <- AddModuleScore(cds_ect, list(treachercollinsgenes), name = "TreacherCollins_ModuleScore")
cds_ect <- AddModuleScore(cds_ect, list(alx), name = "ALX_ModuleScore")

for (n in names(facial_variation_gene_lists))
  cds_cncc <- AddModuleScore(cds_mes, list(facial_variation_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))

for (n in names(other_gene_lists))
  cds_cncc <- AddModuleScore(cds_mes, list(other_gene_lists[[n]]), name = gsub(" ",".",paste(n,"ModuleScore", sep = "_")))
cds_cncc <- AddModuleScore(cds_cncc, list(anterior_hox), name = "Anterior_ModuleScore")
cds_cncc <- AddModuleScore(cds_cncc, list(central_hox), name = "Central_ModuleScore")
cds_cncc <- AddModuleScore(cds_cncc, list(posterior_hox), name = "Posterior_ModuleScore")
cds_cncc <- AddModuleScore(cds_cncc, list(cnccMarkers), name = "CNCC_ModuleScore")
cds_cncc <- AddModuleScore(cds_cncc, list(treachercollinsgenes), name = "TreacherCollins_ModuleScore")
cds_cncc <- AddModuleScore(cds_cncc, list(alx), name = "ALX_ModuleScore")



p1 <- VlnPlot(cds1, "CNCC_ModuleScore1", flip = F) & NoLegend() & coord_flip()
p1

pdf(file = "plots/cncc_module_score.pdf", h = 11, w = 8.5)
p1
dev.off()

DotPlot(cds1, features = cnccMarkers, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds1, features = cnccMarkers_fulllist, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds2, features = hox_genes, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds2, features = treachercollinsgenes, cluster.idents = TRUE ) + RotatedAxis()

pdf(file = "hox_subtype_dot_plot.pdf", h = 8.5, w = 11)
DotPlot(cds2, features = hox_genes, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds2, features = anterior_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds2, features = central_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds2, features = posterior_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_mes, features = hox_genes, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_mes, features = anterior_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_mes, features = central_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_mes, features = posterior_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_ect, features = hox_genes, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_ect, features = anterior_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_ect, features = central_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_ect, features = posterior_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_cncc, features = hox_genes, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_cncc, features = anterior_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_cncc, features = central_hox, cluster.idents = TRUE ) + RotatedAxis()
DotPlot(cds_cncc, features = posterior_hox, cluster.idents = TRUE ) + RotatedAxis()
dev.off()

pdf(file = "hox_subtype_violin_plot.pdf", h = 8.5, w = 11)
VlnPlot(cds2, "Anterior_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds2, "Central_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds2, "Posterior_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds2, features = anterior_hox, stack = TRUE, sort = TRUE, flip = TRUE) + theme(legend.position = "none") + RotatedAxis()
VlnPlot(cds2, "CNCC_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_mes, "Anterior_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_mes, "Central_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_mes, "Posterior_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_mes, features = anterior_hox, stack = TRUE, sort = TRUE, flip = TRUE) + theme(legend.position = "none") + RotatedAxis()
VlnPlot(cds_ect, "Anterior_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_ect, "Central_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_ect, "Posterior_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_ect, features = anterior_hox, stack = TRUE, sort = TRUE, flip = TRUE) + theme(legend.position = "none") + RotatedAxis()
VlnPlot(cds_cncc, "Anterior_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_cncc, "Central_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_cncc, "Posterior_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_cncc, "CNCC_ModuleScore1", sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none")
VlnPlot(cds_cncc, features = anterior_hox, stack = TRUE, sort = TRUE, flip = TRUE) + theme(legend.position = "none") + RotatedAxis()
dev.off()

pdf(file = "hox_subtype_feature_plot.pdf", h = 8.5, w = 11)
FeaturePlot(cds2, "Anterior_ModuleScore1", cols = c("grey90","red"))+NoAxes() + theme(legend.position = "none")
FeaturePlot(cds2, "CNCC_ModuleScore1", cols = c("grey90","red"))+NoAxes() + theme(legend.position = "none")
FeaturePlot(cds_mes, "Anterior_ModuleScore1", cols = c("grey90","red"))+NoAxes() + theme(legend.position = "none")
FeaturePlot(cds_ect, "Anterior_ModuleScore1", cols = c("grey90","red"))+NoAxes() + theme(legend.position = "none")
FeaturePlot(cds_cncc, "Anterior_ModuleScore1", cols = c("grey90","red"))+NoAxes() + theme(legend.position = "none")
FeaturePlot(cds_cncc, "CNCC_ModuleScore1", cols = c("grey90","red"))+NoAxes() + theme(legend.position = "none")
dev.off()

for (n in names(facial_variation_gene_lists)){
  print(VlnPlot(cds1, gsub(" ",".",paste(n,"ModuleScore1", sep = "_")), sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none"))
}
for (n in names(facial_variation_gene_lists)){
  print(VlnPlot(cds2, gsub(" ",".",paste(n,"ModuleScore1", sep = "_")), sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none"))
}
for (n in names(other_gene_lists)){
  print(VlnPlot(cds1, gsub(" ",".",paste(n,"ModuleScore1", sep = "_")), sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none"))
}
for (n in names(other_gene_lists)){
  print(VlnPlot(cds2, gsub(" ",".",paste(n,"ModuleScore1", sep = "_")), sort = TRUE, flip = TRUE) + RotatedAxis() + theme(legend.position = "none"))
}



pdf(file="plots/whole_embryo_main_marker_genes.pdf")
for (n in sort(intersect(rownames(emb.merge),unique(top100pval$gene)))) {
  print(SpatialFeaturePlot(emb.merge, features = n, image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

p1 <- VlnPlot(cds_cncc, features = li_genes, stack = T, flip = T)
p2 <- DimPlot(cds_cncc, label = T)+NoLegend()
p3 <- DimPlot(cds_cncc, group.by = "orig.ident")+NoAxes()
p4 <- DimPlot(cds_cncc,split.by = "stage",label = F)+NoAxes()+NoLegend()
p5 <- FeaturePlot(cds_cncc, features = li_genes, split.by = "stage") & NoAxes()
p6 <- FeaturePlot(cds_cncc, features = li_genes) & NoAxes()

pdf("../plots/cncc_umaps.pdf", h = 8.5, w = 11)
p1
p2
p3
p4
p5
p6
dev.off()



p1 <- DotPlot(object = cds_cncc, features = unique(top_cncc_subtype$gene))+ coord_flip()+ RotatedAxis()
df1 <- data.frame(Gene = p1$data$features.plot, avg.exp = p1$data$avg.exp.scaled, pct.exp = p1$data$pct.exp, cluster = p1$data$id)

p2 <- df1 %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(size = 8, face = "italic" )) +
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p2

pdf(file = "../plots/cncc_marker_gene_dotplots.pdf", h = 11, w = 8.5)
p2
dev.off()


table_samples_by_clusters <- cds_cncc@meta.data %>%
  group_by(stage, tempanno) %>%
  summarize(count = n()) %>%
  spread(tempanno, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('stage', 'total_cell_count', everything())) %>%
  arrange(factor(stage, levels = levels(cds_cncc@meta.data$stage)))

knitr::kable(table_samples_by_clusters)

table_clusters_by_samples <- cds_cncc@meta.data %>%
  group_by(tempanno, stage) %>%
  summarize(count = n()) %>%
  spread(stage, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('tempanno', 'total_cell_count', everything())) %>%
  arrange(factor(tempanno, levels = levels(cds_cncc@meta.data$tempanno)))

knitr::kable(table_clusters_by_samples)

temp_labels <- cds_cncc@meta.data %>%
  group_by(stage) %>%
  tally()

p1 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'stage') %>%
  mutate(stage = factor(stage, levels = levels(cds_cncc@meta.data$stage))) %>%
  ggplot(aes(stage, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = stage, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

temp_labels <- cds_cncc@meta.data %>%
  group_by(tempanno) %>%
  tally() %>%
  dplyr::rename('tempanno' = tempanno)

p2 <- table_clusters_by_samples %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'tempanno') %>%
  mutate(tempanno = factor(tempanno, levels = levels(cds_cncc@meta.data$tempanno))) %>%
  ggplot(aes(tempanno, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = tempanno, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(
  '../plots/cncc_composition_samples_clusters_by_number.pdf',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      cds_cncc@meta.data$stage %>% unique() %>% length(),
      cds_cncc@meta.data$tempanno %>% unique() %>% length()
    )),
  width = 18, height = 8
)

pdf(file="../plots/whole_embryo_cncc_subtype_genes.pdf")
for (n in sort(intersect(rownames(emb.merge),unique(top_cncc_subtype$gene)))) {
  print(SpatialFeaturePlot(emb.merge, features = n, image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

p1 <- VlnPlot(cds_mes, features = li_genes, stack = T, flip = T)
p2 <- DimPlot(cds_mes, label = T)+NoLegend()
p3 <- DimPlot(cds_mes, group.by = "orig.ident")+NoAxes()
p4 <- DimPlot(cds_mes,split.by = "stage",label = F)+NoAxes()+NoLegend()
p5 <- FeaturePlot(cds_mes, features = li_genes, split.by = "stage") & NoAxes()
p6 <- FeaturePlot(cds_mes, features = li_genes) & NoAxes()

pdf("../plots/mes_umaps.pdf", h = 8.5, w = 11)
p1
p2
p3
p4
p5
p6
dev.off()

p1 <- DotPlot(object = cds_mes, features = unique(top_mes_subtype$gene))+ coord_flip()+ RotatedAxis()
df1 <- data.frame(Gene = p1$data$features.plot, avg.exp = p1$data$avg.exp.scaled, pct.exp = p1$data$pct.exp, cluster = p1$data$id)

p2 <- df1 %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(size = 4, face = "italic" )) +
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p2

pdf(file = "../plots/mes_marker_gene_dotplots.pdf", h = 11, w = 8.5)
p2
dev.off()

table_samples_by_clusters <- cds_mes@meta.data %>%
  group_by(stage, MtoHAlltoALL) %>%
  summarize(count = n()) %>%
  spread(MtoHAlltoALL, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('stage', 'total_cell_count', everything())) %>%
  arrange(factor(stage, levels = levels(cds_mes@meta.data$stage)))

knitr::kable(table_samples_by_clusters)

table_clusters_by_samples <- cds_mes@meta.data %>%
  group_by(MtoHAlltoALL, stage) %>%
  summarize(count = n()) %>%
  spread(stage, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('MtoHAlltoALL', 'total_cell_count', everything())) %>%
  arrange(factor(MtoHAlltoALL, levels = levels(cds_mes@meta.data$MtoHAlltoALL)))

knitr::kable(table_clusters_by_samples)

temp_labels <- cds_mes@meta.data %>%
  group_by(stage) %>%
  tally()

p1 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'stage') %>%
  mutate(stage = factor(stage, levels = levels(cds_mes@meta.data$stage))) %>%
  ggplot(aes(stage, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = stage, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

temp_labels <- cds_mes@meta.data %>%
  group_by(MtoHAlltoALL) %>%
  tally() %>%
  dplyr::rename('MtoHAlltoALL' = MtoHAlltoALL)

p2 <- table_clusters_by_samples %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'MtoHAlltoALL') %>%
  mutate(MtoHAlltoALL = factor(MtoHAlltoALL, levels = levels(cds_mes@meta.data$MtoHAlltoALL))) %>%
  ggplot(aes(MtoHAlltoALL, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = MtoHAlltoALL, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(
  '../plots/mes_composition_samples_clusters_by_number.pdf',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      cds_mes@meta.data$stage %>% unique() %>% length(),
      cds_mes@meta.data$MtoHAlltoALL %>% unique() %>% length()
    )),
  width = 18, height = 8
)

pdf(file="plots/whole_embryo_mes_subtype_genes.pdf")
for (n in sort(intersect(rownames(emb.merge),unique(top_mes_subtype$gene)))) {
  print(SpatialFeaturePlot(emb.merge, features = n, image.alpha = 0.3) + theme(legend.position = "right")) 
}
dev.off()

p1 <- DotPlot(object = cds_ect, features = unique(top_ect_subtype$gene))+ coord_flip()+ RotatedAxis()
df1 <- data.frame(Gene = p1$data$features.plot, avg.exp = p1$data$avg.exp.scaled, pct.exp = p1$data$pct.exp, cluster = p1$data$id)

p2 <- df1 %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(size = 4, angle = 45, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(size = 4, face = "italic" )) +
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 
p2

pdf(file = "plots/ect_marker_gene_dotplots.pdf", h = 10, w = 6)
p2
dev.off()

p1 <- VlnPlot(cds_ect, features = li_genes, stack = T, flip = T)
p2 <- DimPlot(cds_ect, label = T)+NoLegend()
p3 <- DimPlot(cds_ect, group.by = "orig.ident")+NoAxes()
p4 <- DimPlot(cds_ect,split.by = "stage",label = F)+NoAxes()+NoLegend()
p5 <- FeaturePlot(cds_ect, features = li_genes, split.by = "stage") & NoAxes()
p6 <- FeaturePlot(cds_ect, features = li_genes) & NoAxes()

pdf("../plots/ect_umaps.pdf", h = 8.5, w = 11)
p1
p2
p3
p4
p5
p6
dev.off()

table_samples_by_clusters <- cds_ect@meta.data %>%
  group_by(stage, tempanno) %>%
  summarize(count = n()) %>%
  spread(tempanno, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('stage', 'total_cell_count', everything())) %>%
  arrange(factor(stage, levels = levels(cds_ect@meta.data$stage)))

knitr::kable(table_samples_by_clusters)

table_clusters_by_samples <- cds_ect@meta.data %>%
  group_by(tempanno, stage) %>%
  summarize(count = n()) %>%
  spread(stage, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('tempanno', 'total_cell_count', everything())) %>%
  arrange(factor(tempanno, levels = levels(cds_ect@meta.data$tempanno)))

knitr::kable(table_clusters_by_samples)

temp_labels <- cds_ect@meta.data %>%
  group_by(stage) %>%
  tally()

p1 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'stage') %>%
  mutate(stage = factor(stage, levels = levels(cds_ect@meta.data$stage))) %>%
  ggplot(aes(stage, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = stage, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_y_continuous(name = 'Number of cells', labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

temp_labels <- cds_ect@meta.data %>%
  group_by(tempanno) %>%
  tally() %>%
  dplyr::rename('tempanno' = tempanno)

p2 <- table_clusters_by_samples %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'tempanno') %>%
  mutate(tempanno = factor(tempanno, levels = levels(cds_ect@meta.data$tempanno))) %>%
  ggplot(aes(tempanno, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = tempanno, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(
  'plots/ect_composition_samples_clusters_by_number.pdf',
  p1 + p2 +
    plot_layout(ncol = 2, widths = c(
      cds_ect@meta.data$stage %>% unique() %>% length(),
      cds_ect@meta.data$tempanno %>% unique() %>% length()
    )),
  width = 18, height = 8
)



#plots for genes identify by denovolylzer anlaysis
pdf("plots/denovlolyzer_nap_clustered_dotplot.pdf", h = 8, w = 12)
Clustered_DotPlot(cds_ect, features = gsub(" ","",strsplit(genelist_gmkf$NaP, ",", )[[1]]), row_label_fontface = "italic", column_label_size = 16, row_label_size = 16, plot_km_elbow = F)
dev.off()

nap_denovo_genes <- gsub(" ","",strsplit(genelist_gmkf$NaP, ",", )[[1]])
nap_denovo_entrez <- annotLookup$entrezgene_id[match(as.character(nap_denovo_genes),annotLookup$hgnc_symbol)]
nap_denovo_go <- enrichGO(nap_denovo_entrez,
                          OrgDb    = org.Hs.eg.db,
                          ont = "ALL",
                          pAdjustMethod = "BH",
                          readable = TRUE)
nap_denovo_dgn <- enrichDGN(nap_denovo_entrez,
                          pAdjustMethod = "BH",
                          readable = TRUE)

nap_denovo_dgn <- pairwise_termsim(nap_denovo_dgn)

human_ect_subtype_markers <- readxl::read_xlsx("tables/SigGenes_Cell_ect_SubTypes.xlsx", sheet = 1)
human_nap <- subset(human_ect_subtype_markers, cluster == "NaP")
foldChange <- human_nap$avg_log2FC
names(foldChange) <- human_nap$Entrez

dotplot(nap_denovo_go, showCategory=30) + ggtitle("CLP Trio De Novo NaP Markers")
cnetplot(nap_denovo_go, categorySize="pvalue", foldChange = foldChange)
cnetplot(nap_denovo_go, categorySize="pvalue", foldChange = foldChange, circular = TRUE, colorEdge = TRUE)
heatplot(nap_denovo_go, foldChange=foldChange, showCategory=25)
heatplot(nap_denovo_dgn, showCategory=25, foldChange=foldChange) + viridis::scale_fill_viridis() + theme(axis.text.x = element_text(face = "italic", size = 12))

pdf("plots/nap_denovo_dgn_heatplot.pdf", h = 10, w = 10)
heatplot(nap_denovo_dgn, showCategory=25, foldChange=foldChange) + viridis::scale_fill_viridis() + theme(axis.text.x = element_text(face = "italic", size = 12))
dev.off()











