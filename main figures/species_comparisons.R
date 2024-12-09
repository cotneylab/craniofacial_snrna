library("ggVennDiagram")
library(eulerr)
library("GeneOverlap")
human_main_markers <- readxl::read_xlsx("tables/SigGenes_Main_CellTypes.xlsx", sheet = 1)
mouse_main_markers <- readxl::read_xlsx("tables/Mouse_SigGenes_Main_CellTypes.xlsx", sheet = 1)
mouse_main_markers_orthologs <- subset(mouse_main_markers, hgnc_symbol != "NA")
#get total number of markers that have human orthologs
unique(mouse_main_markers_orthologs$hgnc_symbol)

human_subtype_markers <- readxl::read_xlsx("tables/SigGenes_CellSubTypes.xlsx", sheet = 1)
mouse_subtype_markers <- readxl::read_xlsx("tables/Mouse_SigGenes_CellSubTypes.xlsx", sheet =1)
mouse_subtype_markers_orthologs <- subset(mouse_subtype_markers, hgnc_symbol != "NA")


#make list of marker genes
human_mesenchyme <- subset(human_main_markers, cluster == "mesenchyme")
human_cncc <- subset(human_main_markers, cluster == "CNCC")
human_ectoderm <- subset(human_main_markers, cluster == "ectoderm")
human_muscle <- subset(human_main_markers, cluster == "muscle")
human_endothelium <- subset(human_main_markers, cluster == "endothelium")
human_blood <- subset(human_main_markers, cluster == "red blood cells")
human_immune <- subset(human_main_markers, cluster == "other blood cells")

mouse_mesenchyme <- subset(mouse_main_markers_orthologs, cluster == "mesenchyme")
mouse_cncc <- subset(mouse_main_markers_orthologs, cluster == "CNCC")
mouse_ectoderm <- subset(mouse_main_markers_orthologs, cluster == "ectoderm")
mouse_muscle <- subset(mouse_main_markers_orthologs, cluster == "muscle")
mouse_endothelium <- subset(mouse_main_markers_orthologs, cluster == "endothelium")
mouse_blood <- subset(mouse_main_markers_orthologs, cluster == "red blood cells")
mouse_immune <- subset(mouse_main_markers_orthologs, cluster == "other blood cells")

mesenchyme_gene_list <- list("Human" = human_mesenchyme$Entrez, "Mouse" = mouse_mesenchyme$Human_Entrez)
cncc_gene_list <- list("Human" = human_cncc$Entrez, "Mouse" = mouse_cncc$Human_Entrez)
ectoderm_gene_list <- list("Human" = human_ectoderm$Entrez, "Mouse" = mouse_ectoderm$Human_Entrez)
muscle_gene_list <- list("Human" = human_muscle$Entrez, "Mouse" = mouse_muscle$Human_Entrez)
endothelium_gene_list <- list("Human" = human_endothelium$Entrez, "Mouse" = mouse_endothelium$Human_Entrez)
blood_gene_list <- list("Human" = human_blood$Entrez, "Mouse" = mouse_blood$Human_Entrez)
immune_gene_list <- list("Human" = human_immune$Entrez, "Mouse" = mouse_immune$Human_Entrez)
human_markers <- list("Human Mesenchyme" = human_mesenchyme$Entrez, "Human CNCC" = human_cncc$Entrez, "Human Muscle" = human_muscle$Entrez, "Human endothelium" = human_endothelium$Entrez, "Human Ectoderm" = human_ectoderm$Entrez, "Human Blood" = human_blood$Entrez, "Human Immune" = human_immune$Entrez)
mouse_markers <- list("Mouse Mesenchyme" = mouse_mesenchyme$Human_Entrez, "Mouse CNCC" = mouse_cncc$Human_Entrez, "Mouse Muscle" = mouse_muscle$Human_Entrez, "Mouse endothelium" = mouse_endothelium$Human_Entrez, "Mouse Ectoderm" = mouse_ectoderm$Human_Entrez, "Mouse Blood" = mouse_blood$Human_Entrez, "Mouse Immune" = mouse_immune$Human_Entrez)
all_markers <- list("Human Mesenchyme" = human_mesenchyme$Entrez, "Human CNCC" = human_cncc$Entrez, "Human Muscle" = human_muscle$Entrez, "Human endothelium" = human_endothelium$Entrez, "Human Ectoderm" = human_ectoderm$Entrez, "Human Blood" = human_blood$Entrez, "Human Immune" = human_immune$Entrez, "Mouse Mesenchyme" = mouse_mesenchyme$Human_Entrez, "Mouse CNCC" = mouse_cncc$Human_Entrez, "Mouse Muscle" = mouse_muscle$Human_Entrez, "Mouse endothelium" = mouse_endothelium$Human_Entrez, "Mouse Ectoderm" = mouse_ectoderm$Human_Entrez, "Mouse Blood" = mouse_blood$Human_Entrez, "Mouse Immune" = mouse_immune$Human_Entrez )
human_marker_background <- unique(human_main_markers$Entrez)


df1 <- human_subtype_markers[,c(1,10)]
df1sample <- split(df1$gene,df1$cluster)
length(df1sample)
human_sub_types <- names(df1sample)
human_subtype_genelist <- as.list(df1sample)
#from above



ggVennDiagram(all_markers)
ggVennDiagram(human_markers)
gs.markers = 7504

go.obj <- newGeneOverlap(human_mesenchyme$Entrez,
                         mouse_mesenchyme$Human_Entrez,
                         genome.size=gs.markers)
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

gom.obj <- newGOM(human_markers, mouse_markers, gs.markers)
drawHeatmap(gom.obj, adj.p=TRUE, grid.col = "Reds", note.col = "white")
pdf("plots/conservation_main_types_marker_gene_overlap.pdf", w = 10, h = 10)
drawHeatmap(gom.obj, adj.p=TRUE, grid.col = "Reds", note.col = "white")
dev.off()
getMatrix(gom.obj, name="pval")
getMatrix(gom.obj, name="odds.ratio")
getMatrix(gom.obj, name="intersection")
conserved_genes <- getNestedList(gom.obj, name="intersection")
conserved_genelist <- list("Mesenchyme" = conserved_genes$`Mouse Mesenchyme`$`Human Mesenchyme`, "CNCC" = conserved_genes$`Mouse CNCC`$`Human CNCC`, "Ectoderm" = conserved_genes$`Mouse Ectoderm`$`Human Ectoderm`, "Muscle" = conserved_genes$`Mouse Muscle`$`Human Muscle`, "Blood" = conserved_genes$`Mouse Blood`$`Human Blood`, "Immune" = conserved_genes$`Mouse Immune`$`Human Immune`)

human_specific_mesenchyme <- setdiff(human_mesenchyme$Entrez, conserved_genelist$Mesenchyme)
human_specific_ectoderm <- setdiff(human_ectoderm$Entrez, conserved_genelist$Ectoderm)
human_specific_muscle <- setdiff(human_muscle$Entrez, conserved_genelist$Muscle)
human_specific_endothelium <- setdiff(human_immune$Entrez, conserved_genelist$Endothelium)
human_specific_blood <- setdiff(human_blood$Entrez, conserved_genelist$Blood)
human_specific_cncc <- setdiff(human_cncc$Entrez, conserved_genelist$CNCC)
human_specific_immune <- setdiff(human_immune$Entrez, conserved_genelist$Immune)

human_specific_genelist <- list("Mesenchyme" = human_specific_mesenchyme, "CNCC" = human_specific_cncc, "Ectoderm" = human_specific_ectoderm, "Endothelium" = human_specific_endothelium, "Muscle" = human_specific_muscle, "Blood" = human_specific_blood, "Immune" = human_specific_immune)

mouse_specific_mesenchyme <- setdiff(mouse_mesenchyme$Human_Entrez, conserved_genelist$Mesenchyme)
mouse_specific_ectoderm <- setdiff(mouse_ectoderm$Human_Entrez, conserved_genelist$Ectoderm)
mouse_specific_muscle <- setdiff(mouse_muscle$Human_Entrez, conserved_genelist$Muscle)
mouse_specific_endothelium <- setdiff(mouse_immune$Human_Entrez, conserved_genelist$Endothelium)
mouse_specific_blood <- setdiff(mouse_blood$Human_Entrez, conserved_genelist$Blood)
mouse_specific_cncc <- setdiff(mouse_cncc$Human_Entrez, conserved_genelist$CNCC)
mouse_specific_immune <- setdiff(mouse_immune$Human_Entrez, conserved_genelist$Immune)

mouse_specific_genelist <- list("Mesenchyme" = mouse_specific_mesenchyme, "CNCC" = mouse_specific_cncc, "Ectoderm" = mouse_specific_ectoderm, "Endothelium" = mouse_specific_endothelium, "Muscle" = mouse_specific_muscle, "Blood" = mouse_specific_blood, "Immune" = mouse_specific_immune)

mesenchyme_comparisons_genelist <- list ("Human Specific" = human_specific_mesenchyme, "Conserved " = conserved_genelist$Mesenchyme, "Mouse Specific" = mouse_specific_mesenchyme)
ectoderm_comparisons_genelist <- list ("Human Specific" = human_specific_ectoderm, "Conserved " = conserved_genelist$Ectoderm, "Mouse Specific" = mouse_specific_ectoderm)
muscle_comparisons_genelist <- list ("Human Specific" = human_specific_muscle, "Conserved " = conserved_genelist$Muscle, "Mouse Specific" = mouse_specific_muscle)
endothelium_comparisons_genelist <- list ("Human Specific" = human_specific_endothelium, "Conserved " = conserved_genelist$Endothelium, "Mouse Specific" = mouse_specific_endothelium)
blood_comparisons_genelist <- list ("Human Specific" = human_specific_blood, "Conserved " = conserved_genelist$Blood, "Mouse Specific" = mouse_specific_blood)
cncc_comparisons_genelist <- list ("Human Specific" = human_specific_cncc, "Conserved " = conserved_genelist$CNCC, "Mouse Specific" = mouse_specific_cncc)
immune_comparisons_genelist <- list ("Human Specific" = human_specific_immune, "Conserved " = conserved_genelist$Immune, "Mouse Specific" = mouse_specific_immune)


library("clusterProfiler")
library(enrichplot)
library("org.Hs.eg.db")
library("AnnotationHub")
library(ReactomePA)
library(DOSE)

BPclusterplot <- compareCluster(geneCluster = conserved_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = conserved_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = conserved_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = conserved_genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = conserved_genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "hsa")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)


DOclusterplot <- compareCluster(geneCluster = conserved_genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
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

pdf(file = "conservation_main_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
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
pdf(file = "conservation_main_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()


go.ls <- conserved_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_bp_conserved_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- conserved_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_cc_conserved_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- conserved_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_mf_conserved_conserved_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- conserved_genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "disease_onology_conserved_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- conserved_genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "kegg_conserved_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- conserved_genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "reactome_conserved_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)



BPclusterplot <- compareCluster(geneCluster = human_specific_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = human_specific_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = human_specific_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = human_specific_genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = human_specific_genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "hsa")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)


DOclusterplot <- compareCluster(geneCluster = human_specific_genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

human_specific_mesenchyme_dgn <- enrichDGN(human_specific_mesenchyme,readable = TRUE, pvalueCutoff=0.05, pAdjustMethod = "BH")
  


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

pdf(file = "human_specific_main_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
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
pdf(file = "human_specific_main_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()


go.ls <- human_specific_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_bp_human_specific_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- human_specific_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_cc_human_specific_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- human_specific_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_mf_human_specific_human_specific_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- human_specific_genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "disease_onology_human_specific_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- human_specific_genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "kegg_human_specific_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- human_specific_genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "reactome_human_specific_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

foldChange <- human_mesenchyme$avg_log2FC
names(foldChange) <- human_mesenchyme$Entrez


sample(human_specific_mesenchyme_dgn$Des)

selected_pathways <- subset(human_specific_mesenchyme_dgn@result, zScore > 4 | grepl('ALX3', geneID))$Description

selected_pathways <- c("Aplasia/Hypoplasia of the iris", "Sloping forehead", "Large nose", "Biparietal narrowing", "Microphthalmos", "Low set ears")
mutate(human_specific_mesenchyme_dgn, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
options(enrichplot.colours = c("red","grey"))
dotplot(human_specific_mesenchyme_dgn, showCategory=30) + ggtitle("Disease Ontology Human Specific Mesenchyme")
options(enrichplot.colours = c("grey","red"))
cnetplot(human_specific_mesenchyme_dgn, categorySize="pvalue", foldChange = foldChange)
cnetplot(human_specific_mesenchyme_dgn, showCategory = selected_pathways, categorySize="pvalue", foldChange = foldChange)
p1 <- cnetplot(human_specific_mesenchyme_dgn, showCategory = selected_pathways, categorySize="pvalue", foldChange = foldChange)

heatplot(human_specific_mesenchyme_dgn, showCategory=5, foldChange = foldChange)

pdf("plots/human_specific_mesenchyme_cnet.pdf", h = 8.5, w = 11)
p1
dev.off()




BPclusterplot <- compareCluster(geneCluster = mouse_specific_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = mouse_specific_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = mouse_specific_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = mouse_specific_genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = mouse_specific_genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "hsa")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)


DOclusterplot <- compareCluster(geneCluster = mouse_specific_genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
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

pdf(file = "mouse_specific_main_cluster_enrichment_dotplot.pdf", w = 17, h = 22)
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
pdf(file = "mouse_specific_main_cluster_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()


go.ls <- mouse_specific_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_bp_mouse_specific_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- mouse_specific_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_cc_mouse_specific_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)

go.ls <- mouse_specific_genelist %>% map(~{
  
  
  
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

openxlsx::write.xlsx(go.ls, file = "gene_onology_mf_mouse_specific_mouse_specific_main_clusters.xlsx", sheetName = names(go.ls), rowNames = FALSE)


do.ls <- mouse_specific_genelist %>% map(~{
  
  
  
  eDO <- enrichDGN(
    gene          = .x,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  return(eDO)
  
})

openxlsx::write.xlsx(do.ls, file = "disease_onology_mouse_specific_main_clusters.xlsx", sheetName = names(do.ls), rowNames = FALSE)

kegg.ls <- mouse_specific_genelist %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.05, 
    organism = 'hsa'
  )
  return(eKEGG)
})

openxlsx::write.xlsx(kegg.ls, file = "kegg_mouse_specific_main_clusters.xlsx", sheetName = names(kegg.ls), rowNames = FALSE)

pathway.ls <- mouse_specific_genelist %>% map(~{
  ePathway <- enrichPathway(
    gene = .x,
    pvalueCutoff = 0.05,
    readable = TRUE, 
  )
  return(ePathway)
})

openxlsx::write.xlsx(pathway.ls, file = "reactome_mouse_specific_main_clusters.xlsx", sheetName = names(pathway.ls), rowNames = FALSE)

BPclusterplot <- compareCluster(geneCluster = mesenchyme_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = mesenchyme_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = mesenchyme_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = mesenchyme_comparisons_genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = mesenchyme_comparisons_genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "hsa")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)


DOclusterplot <- compareCluster(geneCluster = mesenchyme_comparisons_genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 20, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 20, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 20, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 20, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 20, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 20, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))


cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])

pdf(file = "mesenchyme_comparison_enrichment_dotplot.pdf", w = 17, h = 22)
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
pdf(file = "mesenchyme_comparison_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

BPclusterplot <- compareCluster(geneCluster = ectoderm_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = ectoderm_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = ectoderm_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = ectoderm_comparisons_genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = ectoderm_comparisons_genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "hsa")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)


DOclusterplot <- compareCluster(geneCluster = ectoderm_comparisons_genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 20, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 20, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 20, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 20, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 20, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 20, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))


cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])

pdf(file = "ectoderm_comparison_enrichment_dotplot.pdf", w = 17, h = 22)
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
pdf(file = "ectoderm_comparison_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

BPclusterplot <- compareCluster(geneCluster = cncc_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP')
CCclusterplot <- compareCluster(geneCluster = cncc_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC')
MFclusterplot <- compareCluster(geneCluster = cncc_comparisons_genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF')
BPclusterplot <- pairwise_termsim(BPclusterplot)
CCclusterplot <- pairwise_termsim(CCclusterplot)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = cncc_comparisons_genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = cncc_comparisons_genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "hsa")
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)


DOclusterplot <- compareCluster(geneCluster = cncc_comparisons_genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH")
DOclusterplot <- pairwise_termsim(DOclusterplot)

options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 20, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 20, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 20, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 20, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 20, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 20, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))


cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol=2, labels=LETTERS[1:6])

pdf(file = "cncc_comparison_enrichment_dotplot.pdf", w = 17, h = 22)
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
pdf(file = "cncc_comparison_enrichment_emapplot.pdf", h = 17, w = 22)
cowplot::plot_grid(p1, p2, p3, p4,p5,p6, ncol=2, labels=LETTERS[1:6])
dev.off()

#gene symbol version
mesenchyme_gene_list <- list("Human" = human_mesenchyme$gene, "Mouse" = mouse_mesenchyme$hgnc_symbol)
cncc_gene_list <- list("Human" = human_cncc$gene, "Mouse" = mouse_cncc$hgnc_symbol)
ectoderm_gene_list <- list("Human" = human_ectoderm$gene, "Mouse" = mouse_ectoderm$hgnc_symbol)
muscle_gene_list <- list("Human" = human_muscle$gene, "Mouse" = mouse_muscle$hgnc_symbol)
endothelium_gene_list <- list("Human" = human_endothelium$gene, "Mouse" = mouse_endothelium$hgnc_symbol)
blood_gene_list <- list("Human" = human_blood$gene, "Mouse" = mouse_blood$hgnc_symbol)
immune_gene_list <- list("Human" = human_immune$gene, "Mouse" = mouse_immune$hgnc_symbol)
human_markers <- list("Human Mesenchyme" = human_mesenchyme$gene, "Human CNCC" = human_cncc$gene, "Human Muscle" = human_muscle$gene, "Human endothelium" = human_endothelium$gene, "Human Ectoderm" = human_ectoderm$gene, "Human Blood" = human_blood$gene, "Human Immune" = human_immune$gene)
mouse_markers <- list("Mouse Mesenchyme" = mouse_mesenchyme$hgnc_symbol, "Mouse CNCC" = mouse_cncc$hgnc_symbol, "Mouse Muscle" = mouse_muscle$hgnc_symbol, "Mouse endothelium" = mouse_endothelium$hgnc_symbol, "Mouse Ectoderm" = mouse_ectoderm$hgnc_symbol, "Mouse Blood" = mouse_blood$hgnc_symbol, "Mouse Immune" = mouse_immune$hgnc_symbol)
all_markers <- list("Human Mesenchyme" = human_mesenchyme$gene, "Human CNCC" = human_cncc$gene, "Human Muscle" = human_muscle$gene, "Human endothelium" = human_endothelium$gene, "Human Ectoderm" = human_ectoderm$gene, "Human Blood" = human_blood$gene, "Human Immune" = human_immune$gene, "Mouse Mesenchyme" = mouse_mesenchyme$hgnc_symbol, "Mouse CNCC" = mouse_cncc$hgnc_symbol, "Mouse Muscle" = mouse_muscle$hgnc_symbol, "Mouse endothelium" = mouse_endothelium$hgnc_symbol, "Mouse Ectoderm" = mouse_ectoderm$hgnc_symbol, "Mouse Blood" = mouse_blood$hgnc_symbol, "Mouse Immune" = mouse_immune$hgnc_symbol )
human_marker_background <- unique(human_main_markers$gene)

ggVennDiagram(all_markers)

gs.markers = 7504

go.obj <- newGeneOverlap(human_mesenchyme$gene,
                         mouse_mesenchyme$hgnc_symbol,
                         genome.size=gs.markers)
go.obj <- testGeneOverlap(go.obj)
print(go.obj)

gom.obj <- newGOM(human_markers, mouse_markers, gs.markers)
drawHeatmap(gom.obj, adj.p=TRUE, grid.col = "Reds", note.col = "white")
getMatrix(gom.obj, name="pval")
getMatrix(gom.obj, name="odds.ratio")
getMatrix(gom.obj, name="intersection")
conserved_genes <- getNestedList(gom.obj, name="intersection")
conserved_genelist <- list("Mesenchyme" = conserved_genes$`Mouse Mesenchyme`$`Human Mesenchyme`, "CNCC" = conserved_genes$`Mouse CNCC`$`Human CNCC`, "Ectoderm" = conserved_genes$`Mouse Ectoderm`$`Human Ectoderm`, "Muscle" = conserved_genes$`Mouse Muscle`$`Human Muscle`, "Blood" = conserved_genes$`Mouse Blood`$`Human Blood`, "Immune" = conserved_genes$`Mouse Immune`$`Human Immune`)

openxlsx::write.xlsx(conserved_genelist, file = "tables/conserved_marker_genes.xlsx", sheetName = names(conserved_genelist), rowNames = FALSE)


human_specific_mesenchyme <- setdiff(human_mesenchyme$gene, conserved_genelist$Mesenchyme)
human_specific_ectoderm <- setdiff(human_ectoderm$gene, conserved_genelist$Ectoderm)
human_specific_muscle <- setdiff(human_muscle$gene, conserved_genelist$Muscle)
human_specific_endothelium <- setdiff(human_immune$gene, conserved_genelist$Endothelium)
human_specific_blood <- setdiff(human_blood$gene, conserved_genelist$Blood)
human_specific_cncc <- setdiff(human_cncc$gene, conserved_genelist$CNCC)
human_specific_immune <- setdiff(human_immune$gene, conserved_genelist$Immune)

human_specific_genelist <- list("Mesenchyme" = human_specific_mesenchyme, "CNCC" = human_specific_cncc, "Ectoderm" = human_specific_ectoderm, "Endothelium" = human_specific_endothelium, "Muscle" = human_specific_muscle, "Blood" = human_specific_blood, "Immune" = human_specific_immune)

openxlsx::write.xlsx(human_specific_genelist, file = "tables/human_specific_marker_genes.xlsx", sheetName = names(human_specific_genelist), rowNames = FALSE)

mouse_specific_mesenchyme <- setdiff(mouse_mesenchyme$hgnc_symbol, conserved_genelist$Mesenchyme)
mouse_specific_ectoderm <- setdiff(mouse_ectoderm$hgnc_symbol, conserved_genelist$Ectoderm)
mouse_specific_muscle <- setdiff(mouse_muscle$hgnc_symbol, conserved_genelist$Muscle)
mouse_specific_endothelium <- setdiff(mouse_immune$hgnc_symbol, conserved_genelist$Endothelium)
mouse_specific_blood <- setdiff(mouse_blood$hgnc_symbol, conserved_genelist$Blood)
mouse_specific_cncc <- setdiff(mouse_cncc$hgnc_symbol, conserved_genelist$CNCC)
mouse_specific_immune <- setdiff(mouse_immune$hgnc_symbol, conserved_genelist$Immune)

mouse_specific_genelist <- list("Mesenchyme" =mouse_specific_mesenchyme, "CNCC" =mouse_specific_cncc, "Ectoderm" =mouse_specific_ectoderm, "Endothelium" =mouse_specific_endothelium, "Muscle" =mouse_specific_muscle, "Blood" =mouse_specific_blood, "Immune" =mouse_specific_immune)
openxlsx::write.xlsx(mouse_specific_genelist, file = "tables/mouse_specific_marker_genes.xlsx", sheetName = names(mouse_specific_genelist), rowNames = FALSE)
