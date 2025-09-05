library(scToppR)
library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
library(openxlsx)
library(orthogene)

# Mouse ectoderm and mesenchyme subclustering

#subset E11 matching TW E11 data
setwd("/scr1/users/manchela/Data")
cds <- readRDS("face_mouse_final.rds")
Idents(cds) <- "stage"
E11 <- subset(cds, idents = "E11")
Idents(E11) <- "cell_type"

require(biomaRt)
mart_m <- useMart("ENSEMBL_MART_ENSEMBL", host = "useast.ensembl.org")
mart_m <- useDataset("mmusculus_gene_ensembl", mart_m)
mGenes <- getBM(
  mart = mart_m,
  attributes = c(
    "mgi_symbol","description",
    "entrezgene_id",
    "ensembl_gene_id",
    "gene_biotype"),
  filter = "mgi_symbol",
  values = rownames(cds),
  uniqueRows=TRUE)

tf_table <- read.table("human_tf_Lambertetal_PMID29425488.txt",sep="\t",header = T)
mart <-  useMart("ENSEMBL_MART_ENSEMBL") #, host = "https://useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id","ensembl_gene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = tf_table$Human_TF,
  uniqueRows=TRUE)
annotLookup_ortho_one2one <- orthogene::convert_orthologs(unique(annotLookup$hgnc_symbol),input_species="human",output_species = "mouse",
                                                          non121_strategy = "drop_both_species", method="gprofiler") 
annotLookup$mgi_symbol <- rownames(annotLookup_ortho_one2one)[match(annotLookup$hgnc_symbol,annotLookup_ortho_one2one$input_gene)]
mGenes$TF <- "No"
mGenes$TF[which(mGenes$mgi_symbol%in%na.omit(annotLookup$mgi_symbol))] <- "Yes"

### MESENCHYME ###

#subset E11 mesenchyme and transfer subtype annotation from TW
mes <- subset(E11, idents = "mesenchyme")
TW_mes <- readRDS("/Volumes/Extreme SSD/Mouse_cellranger/figures_all/TW/mm.rds")
anchors <- FindTransferAnchors(reference = TW_mes, query = mes, reduction = 'cca', normalization.method = "LogNormalize", dims = 1:30)
predicted.id <- TransferData(anchorset = anchors, refdata = TW_mes$RNA_snn_res.1, weight.reduction = mes[['harmony']],dims = 1:30)
mes$cellType1 <- predicted.id$predicted.id
predicted.id <- TransferData(anchorset = anchors, refdata = TW_mes$cellType2, weight.reduction = mes[['harmony']],dims = 1:30)
mes$cellType2 <- predicted.id$predicted.id
DimPlot(mes, group.by = "cellType2")

#pre-processing
mes <- NormalizeData(mes) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")
mes <- RunUMAP(mes, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1), verbose = FALSE)
Idents(mes) <- "RNA_snn_res.0.8"

#plotting
library(clustree)
p4 <- clustree(mes, prefix = "RNA_snn_res.")
p1 <- DimPlot(mes, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(mes, group.by = "cellType2")+NoAxes()
p3 <- DimPlot(mes,split.by = "stage",label = T)+NoAxes()+NoLegend()
pdf("figures/mes_E11_TWtransferID.pdf",width = 16,height = 8)
p1+p2
p4
dev.off()

# find markers
sig_genes <- FindAllMarkers(mes, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
sig_genes$TF <- "no"
id <- intersect(mGenes$mgi_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$TF[sig_genes$gene %in% id] <- mGenes$TF[match(x$gene,mGenes$mgi_symbol)]
sig_genes$Entrez <- mGenes$entrezgene_id[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes$description <- mGenes$description[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes <- sig_genes[, c("gene","TF","Entrez","description","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
#remove "Rik" and "Gm" genes from markers lists
RikFeatures <- grep("Rik",rownames(count))
GmFeatures <- grep("Gm",rownames(count))
count1 <- count[-RikFeatures,]
count1 <- count1[-GmFeatures,]
sig_genes2 <- filter(sig_genes, gene %in% row.names(count1))
top <- sig_genes2 %>% group_by(cluster) %>% top_n(5, avg_log2FC) #%>% top_n(5, pct.diff);top
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

heatmap(table(mes$cellType1, mes$RNA_snn_res.0.8))
mes <- RenameIdents(mes,c("0"="surface.fusion","1"="pLNP1","2"="palatal shelf1","3"="pLNP2","4"="chondroprogenitors",
                          "5"="palatal shelf2","6"="MxP1","7"="MxP2","8"="MxP.aLNP","9"="ambiguous","10"="fusion.mes2",
                          "11"="surface.mes","12"="palatal.fusion","13"="aLNP","14"="hrmn","15"="oxy","16"="cyc"))


#use scToppR for validating annotations
setwd("/scr1/users/manchela/Data")
mes_degs <- openxlsx::read.xlsx("Craniofacial-scRNAseqManuscript/final-tables/Mouse_SigGenes_mes_SubTypes.xlsx")
mes_toppData <- toppFun(mes_degs,
                    gene_col = "gene",
                    cluster_col = "cluster",
                    p_val_col = "p_val_adj",
                    logFC_col = "avg_log2FC")
head(mes_toppData)

p <- toppBalloon(mes_toppData,
            categories = "MousePheno",x_axis_text_size = 10) + coord_flip() +
  viridis::scale_color_viridis(option = "C",direction=-1) +
  ggtitle("Mouse Mesenchyme Subtypes: Mouse Phenotype ToppGene") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(color="black"),
        plot.title = element_text(color="black",size=14,hjust=0.5))

pdf("Craniofacial-scRNAseqManuscript/final-figs/Mouse_mes_subtypes_scToppR_baloonplot.pdf",height=14,width=8.5)
p
dev.off()


# functional annotation
mmes1 <- readRDS("mmes1_final.rds")
mart <-  useMart("ENSEMBL_MART_ENSEMBL")#, host = "https://asia.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
mes_annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id", "ensembl_gene_id",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = rownames(mmes1),
  uniqueRows=TRUE)

mes_degs$entrez <- mes_annotLookup$entrezgene_id[match(mes_degs$gene,mes_annotLookup$hgnc_symbol)]
mes_degs$biotype <- mes_annotLookup$gene_biotype[match(mes_degs$gene,mes_annotLookup$hgnc_symbol)]

top_mes <- mes_degs %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top100_mes <- mes_degs %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(100, avg_log2FC)
top100pval_mes <- subset(top100_mes, rowSums(top100_mes[5] < 0.05) > 0)

df1 <- top100pval_mes[,c(6,10)]
df1$cluster <- factor(df1$cluster,levels=unique(df1$cluster))
df1sample <- split(df1$entrez,df1$cluster)
length(df1sample)
initial_types <- names(df1sample)

genelist <- as.list(df1sample)
background_genes <- unique(mes_degs$entrez)

BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP', universe = background_genes)
BPclusterplot <- pairwise_termsim(BPclusterplot)

CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC', universe = background_genes)
CCclusterplot <- pairwise_termsim(CCclusterplot)

MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF', universe = background_genes)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)

DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
DOclusterplot <- pairwise_termsim(DOclusterplot)


options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 5, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 5, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 5, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 5, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 5, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))

p1
p2
p3
p4
p5
p6

mes_annot_list <- list(BP =BPclusterplot, CC=CCclusterplot, MF = MFclusterplot, Reactome=Pathwayclusterplot,
                       KEGG=KEGGclusterplot, DisGeNet=DOclusterplot)
saveRDS(mes_annot_list, "mouse_mes_subtype_functional-annotation_lists.rds")

mes_annot_list <- readRDS("mouse_mes_subtype_functional-annotation_lists.rds")
mes_annot_list$BP <- mes_annot_list$BP@compareClusterResult
mes_annot_list$CC <- mes_annot_list$CC@compareClusterResult
mes_annot_list$MF <- mes_annot_list$MF@compareClusterResult
mes_annot_list$Reactome <- mes_annot_list$Reactome@compareClusterResult
mes_annot_list$KEGG <- mes_annot_list$KEGG@compareClusterResult
mes_annot_list$DisGeNet <- mes_annot_list$DisGeNet@compareClusterResult
openxlsx::write.xlsx(mes_annot_list, file = "Craniofacial-scRNAseqManuscript/final-tables/mouse_mes_subtype_functional-annotation.xlsx",
                     sheetName = names(mes_annot_list), rowNames = FALSE)

### ECTODERM ###

#subset E11 ectoderm and transfer subtype annotation from TW
ect <- subset(E11, idents = "ectoderm")
TW_ect <- readRDS("/Volumes/Extreme SSD/Mouse_cellranger/figures_all/TW/ee.rds")
anchors <- FindTransferAnchors(reference = TW_ect, query = ect, reduction = 'cca', normalization.method = "LogNormalize", dims = 1:30)
predicted.id <- TransferData(anchorset = anchors, refdata = TW_ect$cellType2, weight.reduction = ect[['harmony']],dims = 1:30)
ect$cellType2 <- predicted.id$predicted.id
DimPlot(ect, group.by = "cellType2")
ect <- NormalizeData(ect) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")

ect <- RunUMAP(ect, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.5,0.6,0.8,1), verbose = FALSE)
clustree(ect, prefix = "RNA_snn_res.")
Idents(ect) <- "RNA_snn_res.0.8"
p1 <- DimPlot(ect, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(ect, group.by = "cellType2")+NoAxes()
p3 <- DimPlot(ect,split.by = "stage",label = T)+NoAxes()+NoLegend()
p4 <- clustree(ect, prefix = "RNA_snn_res.")

pdf("figures/ect_E11_TWtransferID.pdf",width = 16,height = 8)
p1+p2
p4
dev.off()

sig_genes <- FindAllMarkers(ect, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
sig_genes$TF <- "no"
id <- intersect(mGenes$mgi_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$TF[sig_genes$gene %in% id] <- mGenes$TF[match(x$gene,mGenes$mgi_symbol)]
sig_genes$Entrez <- mGenes$entrezgene_id[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes$description <- mGenes$description[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes <- sig_genes[, c("gene","TF","Entrez","description","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
#remove "Rik" and "Gm" genes from markers lists
# RikFeatures <- grep("Rik",rownames(count))
# GmFeatures <- grep("Gm",rownames(count))
# count1 <- count[-RikFeatures,]
# count1 <- count1[-GmFeatures,]
sig_genes2 <- filter(sig_genes, gene %in% row.names(count1))
top <- sig_genes2 %>% group_by(cluster) %>% top_n(5, avg_log2FC) #%>% top_n(5, pct.diff);top
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
heatmap(table(ect$cellType2, ect$RNA_snn_res.0.8))
ect <- RenameIdents(ect,c("0"="fusion.surface","1"="palate","2"="OE.NaP","3"="ambiguous","4"="surface","5"="OE1",
                          "6"="dental1","7"="periderm2","8"="OE3","9"="periderm","10"="OE2","11"="dental2","12"="NaP","13"="blood"))
ect <- subset(ect,idents = "blood",invert=T)
DimPlot(ect, label = T)
ect$cellType3 <- Idents(ect)
mes$cellType3 <- Idents(mes)

E11_mes_ect <- merge(mes,ect)
Idents(E11_mes_ect) <- "cellType3"
E11$cell_type3 <- E11_mes_ect$cellType3
Idents(cds) <- "cell_type"

cds$E11ectmes1 <- E11_mes_ect$cellType3

for (i in levels(E11_mes_ect)) {
  id <- rownames(E11_mes_ect@meta.data[E11_mes_ect@meta.data$cellType3 == i,])
  cds <- SetIdent(cds,id,value = i)
}        
cds$E11ectmes <- Idents(cds)

DimPlot(cds, group.by = "E11ectmes")+NoLegend()+NoAxes()
#expand the annotation from E11 mesenchyme to combined mesenchyme 
Idents(cds) <- "cell_type"
mes_all <- subset(cds, idents = "mesenchyme")
mes_all <- NormalizeData(mes_all) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")
ElbowPlot(ect,ndims = 50)
mes_all <- RunUMAP(mes_all, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1), verbose = FALSE)
p4 <- clustree(mes_all, prefix = "RNA_snn_res.")
Idents(mes_all) <- "RNA_snn_res.0.8"
p1 <- DimPlot(mes_all, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(mes_all, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(mes_all,split.by = "stage",label = T)+NoAxes()+NoLegend()
p5 <- DimPlot(mes_all,group.by = "RNA_snn_res.0.8", label = T)+NoAxes()+NoLegend()
p1+p2
pdf("figures/mes_all_subtypes.pdf",width = 16,height = 8)
p1+p2
p5+p4
dev.off()

Idents(mes_all) <- "stage"
mes1 <- subset(mes_all, idents = "E11") 
p1 <- DimPlot(mes_all, group.by = "RNA_snn_res.0.8",label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(mes1, group.by = "E11ectmes",label = T)+NoAxes()
p1+p2
heatmap(table(mes1$E11ectmes1, mes1$RNA_snn_res.0.8))
barplot(prop.table(table(mes_all$stage,mes_all$RNA_snn_res.0.8),margin = 2),col = c("navy","blue","purple","brown","red","orange", "yellow"),legend.text = T)
Idents(mes_all) <- "RNA_snn_res.0.8"
mes_all <- RenameIdents(mes_all,c("0"="palatalShelf1","1"="MxP.surface","2"="LNP","3"="MxP.aLNP","4"="MandibularArch",
                                  "5"="pLNP.fusion","6"="e.osteoblast","7"="palatalShelf2","8"="MxP2","9"="cycling",
                                  "10"="pLNP2","11"="muscle.mes","12"="fusion.mes2","13"="palatal.fusion","14"="palatalShelf2",
                                  "15"="chondrocyte","16"="cartilage","17"="neural","18"="l.Osteoblast","19"="immune"))
mes_all$cellType4 <- Idents(mes_all)
heatmap(table(mes_all$E11ectmes1, mes_all$cellType4))
sig_genes <- FindAllMarkers(mes_all, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
sig_genes$TF <- "no"
id <- intersect(mGenes$mgi_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$TF[sig_genes$gene %in% id] <- mGenes$TF[match(x$gene,mGenes$mgi_symbol)]
sig_genes$Entrez <- mGenes$entrezgene_id[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes$description <- mGenes$description[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes <- sig_genes[, c("gene","TF","Entrez","description","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
#remove "Rik" and "Gm" genes from markers lists
# RikFeatures <- grep("Rik",rownames(count))
# GmFeatures <- grep("Gm",rownames(count))
# count1 <- count[-RikFeatures,]
# count1 <- count1[-GmFeatures,]
sig_genes2 <- filter(sig_genes, gene %in% row.names(count1))
top <- sig_genes2 %>% group_by(cluster) %>% top_n(10, avg_log2FC) #%>% top_n(5, pct.diff);top
p <- DotPlot(object = mes_all, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
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
openxlsx::write.xlsx(sig_genes,"tables/SigGenes_mes_all_res0.8.xlsx")

#expand the annotation from E11 ectoderm to combined ectoderm 
ect_all <- subset(cds, idents = "ectoderm")
ect_all <- NormalizeData(ect_all) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")
ElbowPlot(ect,ndims = 50)
ect_all <- RunUMAP(ect_all, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1), verbose = FALSE)
p4 <- clustree(ect_all, prefix = "RNA_snn_res.")
Idents(ect_all) <- "RNA_snn_res.0.8"

p1 <- DimPlot(ect_all, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(ect_all, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(ect_all,split.by = "stage",label = T)+NoAxes()+NoLegend()
p5 <- DimPlot(ect_all,group.by = "RNA_snn_res.0.8", label = T)+NoAxes()+NoLegend()
p1+p2
pdf("figures/ect_all_subtypes.pdf",width = 16,height = 8)
p1+p2
p5+p4
dev.off()

Idents(ect_all) <- "stage"
ect1 <- subset(ect_all, idents = "E11") 
p1 <- DimPlot(ect_all, group.by = "RNA_snn_res.0.8",label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(ect1, group.by = "E11ectmes",label = T)+NoAxes()
p1+p2
heatmap(table(ect1$E11ectmes1, ect1$RNA_snn_res.0.8))
barplot(prop.table(table(ect_all$stage,ect_all$RNA_snn_res.0.8),margin = 2),col = c("navy","blue","purple","brown","red","orange", "yellow"),legend.text = T)
Idents(ect_all) <- "RNA_snn_res.0.8"
sig_genes <- FindAllMarkers(ect_all, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
sig_genes$TF <- "no"
id <- intersect(mGenes$mgi_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$TF[sig_genes$gene %in% id] <- mGenes$TF[match(x$gene,mGenes$mgi_symbol)]
sig_genes$Entrez <- mGenes$entrezgene_id[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes$description <- mGenes$description[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes <- sig_genes[, c("gene","TF","Entrez","description","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
head(sig_genes)
#remove "Rik" and "Gm" genes from markers lists
# RikFeatures <- grep("Rik",rownames(count))
# GmFeatures <- grep("Gm",rownames(count))
# count1 <- count[-RikFeatures,]
# count1 <- count1[-GmFeatures,]
sig_genes2 <- filter(sig_genes, gene %in% row.names(count1))
top <- sig_genes2 %>% group_by(cluster) %>% top_n(10, avg_log2FC) #%>% top_n(5, pct.diff);top
p <- DotPlot(object = ect_all, features = unique(top$gene))+ coord_flip()+ RotatedAxis()
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

ect_all <- RenameIdents(ect_all,c("0"="fusion.surface","1"="ambiguous","2"="palate","3"="periderm2","4"="dental1",
                                  "5"="OE2","6"="dental2","7"="periderm","8"="surface.ect1","9"="NaP",
                                  "10"="OE.NaP","11"="OE1.Pax7","12"="OE3","13"="surface.ect2","14"="ciliated",
                                  "15"="collagen","16"="cardiac","17"="auditory","18"="pituitary","19"="myeloid"))
ect_all$cellType4 <- Idents(ect_all)
heatmap(table(ect_all$E11ectmes1, ect_all$cellType4))
openxlsx::write.xlsx(sig_genes,"tables/SigGenes_ect_all_res0.8.xlsx")



#use scToppR for validating annotations
ect_degs <- openxlsx::read.xlsx("Craniofacial-scRNAseqManuscript/final-tables/Mouse_SigGenes_ect_SubTypes.xlsx")
ect_toppData <- toppFun(ect_degs,
                        gene_col = "gene",
                        cluster_col = "cluster",
                        p_val_col = "p_val_adj",
                        logFC_col = "avg_log2FC")
head(ect_toppData)

p <- toppBalloon(ect_toppData,
                 categories = "MousePheno",x_axis_text_size = 10) + coord_flip() +
  viridis::scale_color_viridis(option = "C",direction=-1) +
  ggtitle("Mouse Ectoderm Subtypes: Mouse Phenotype ToppGene") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(color="black"),
        plot.title = element_text(color="black",size=14,hjust=0.5))

p

pdf("Craniofacial-scRNAseqManuscript/final-figs/Mouse_ect_subtypes_scToppR_baloonplot.pdf",height=14,width=8.5)
p
dev.off()

## functional annotation
mect1 <- readRDS("mect1_final.rds")

ect_annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol","description",
    "entrezgene_id", "ensembl_gene_id",
    "gene_biotype"),
  filter = "hgnc_symbol",
  values = rownames(mect1),
  uniqueRows=TRUE)

ect_degs$entrez <- ect_annotLookup$entrezgene_id[match(ect_degs$gene,ect_annotLookup$hgnc_symbol)]
ect_degs$biotype <- ect_annotLookup$gene_biotype[match(ect_degs$gene,ect_annotLookup$hgnc_symbol)]

top_mes <- ect_degs %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top100_mes <- ect_degs %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(100, avg_log2FC)
top100pval_mes <- subset(top100_mes, rowSums(top100_mes[5] < 0.05) > 0)

df1 <- top100pval_mes[,c(6,10)]
df1$cluster <- factor(df1$cluster,levels=unique(df1$cluster))
df1sample <- split(df1$entrez,df1$cluster)
length(df1sample)
initial_types <- names(df1sample)

genelist <- as.list(df1sample)
background_genes <- unique(ect_degs$entrez)

BPclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'BP', universe = background_genes)
BPclusterplot <- pairwise_termsim(BPclusterplot)

CCclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'CC', universe = background_genes)
CCclusterplot <- pairwise_termsim(CCclusterplot)

MFclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod = "BH",  ont = 'MF', universe = background_genes)
MFclusterplot <- pairwise_termsim(MFclusterplot)

Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
Pathwayclusterplot <- pairwise_termsim(Pathwayclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
KEGGclusterplot <- pairwise_termsim(KEGGclusterplot)

DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDGN", pvalueCutoff=0.05, pAdjustMethod = "BH", universe = background_genes)
DOclusterplot <- pairwise_termsim(DOclusterplot)


options(enrichplot.colours = c("blue", "grey"))
p1 <- dotplot(BPclusterplot, showCategory = 5, font.size = 6, title = "GO: Biological Process", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("orange", "grey"))
p2 <- dotplot(CCclusterplot, showCategory = 5, font.size = 6, title = "GO: Cellular Component", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("green", "grey"))
p3 <- dotplot(MFclusterplot, showCategory = 5, font.size = 6, title = "GO: Molecular Function", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("brown", "grey"))
p4 <- dotplot(Pathwayclusterplot, showCategory = 5, font.size = 6, title = "Reactome Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("yellow", "grey"))
p5 <- dotplot(KEGGclusterplot, showCategory = 5, font.size = 6, title = "KEGG Pathway", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))
options(enrichplot.colours = c("red", "grey"))
p6 <- dotplot(DOclusterplot, showCategory = 5, font.size = 6, title = "DisGeNet", color = "p.adjust", label_format = 50) + theme( legend.text=element_text(size=5),legend.position = "right",legend.box = "vertical",text = element_text(size=6), plot.title = element_text(hjust= 0.5, size = 12))

p1
p2
p3
p4
p5
p6

ect_annot_list <- list(BP =BPclusterplot, CC=CCclusterplot, MF = MFclusterplot, Reactome=Pathwayclusterplot,
                       KEGG=KEGGclusterplot, DisGeNet=DOclusterplot)
saveRDS(ect_annot_list, "mouse_ect_subtype_functional-annotation_lists.rds")

ect_annot_list <- readRDS("mouse_ect_subtype_functional-annotation_lists.rds")
ect_annot_list$BP <- ect_annot_list$BP@compareClusterResult
ect_annot_list$CC <- ect_annot_list$CC@compareClusterResult
ect_annot_list$MF <- ect_annot_list$MF@compareClusterResult
ect_annot_list$Reactome <- ect_annot_list$Reactome@compareClusterResult
ect_annot_list$KEGG <- ect_annot_list$KEGG@compareClusterResult
ect_annot_list$DisGeNet <- ect_annot_list$DisGeNet@compareClusterResult
openxlsx::write.xlsx(ect_annot_list, file = "Craniofacial-scRNAseqManuscript/final-tables/mouse_ect_subtype_functional-annotation.xlsx",
                     sheetName = names(ect_annot_list), rowNames = FALSE)





