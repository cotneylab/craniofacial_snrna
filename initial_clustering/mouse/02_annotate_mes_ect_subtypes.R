library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
#subset E11 matching TW E11 data
cds <- readRDS("data/cds.rds")
Idents(cds) <- "stage"
E11 <- subset(cds, idents = "E11")
Idents(E11) <- "cell_type"
#subset E11 mesenchyme and transfer subtype annotation from TW
mes <- subset(E11, idents = "mesenchyme")
TW_mes <- readRDS("/Volumes/Extreme SSD/Mouse_cellranger/figures_all/TW/mm.rds")
anchors <- FindTransferAnchors(reference = TW_mes, query = mes, reduction = 'cca', normalization.method = "LogNormalize", dims = 1:30)
predicted.id <- TransferData(anchorset = anchors, refdata = TW_mes$RNA_snn_res.1, weight.reduction = mes[['harmony']],dims = 1:30)
mes$cellType1 <- predicted.id$predicted.id
predicted.id <- TransferData(anchorset = anchors, refdata = TW_mes$cellType2, weight.reduction = mes[['harmony']],dims = 1:30)
mes$cellType2 <- predicted.id$predicted.id
DimPlot(mes, group.by = "cellType2")
mes <- NormalizeData(mes) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% 
  RunPCA(verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident")
mes <- RunUMAP(mes, reduction = "harmony", dims = 1:30, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1,0.2,0.4,0.6,0.8,1), verbose = FALSE)
Idents(mes) <- "RNA_snn_res.0.8"
library(clustree)
p4 <- clustree(mes, prefix = "RNA_snn_res.")
p1 <- DimPlot(mes, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(mes, group.by = "cellType2")+NoAxes()
p3 <- DimPlot(mes,split.by = "stage",label = T)+NoAxes()+NoLegend()
pdf("figures/mes_E11_TWtransferID.pdf",width = 16,height = 8)
p1+p2
p4
dev.off()

mGenes <- readRDS(file="/Users/naghamkhourifarah/Documents/Annotation1.rds")
sig_genes <- FindAllMarkers(mes, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
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

mGenes <- readRDS(file="/Users/naghamkhourifarah/Documents/Annotation1.rds")
sig_genes <- FindAllMarkers(ect, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
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
mGenes <- readRDS(file="/Users/naghamkhourifarah/Documents/Annotation1.rds")
sig_genes <- FindAllMarkers(mes_all, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
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
mGenes <- readRDS(file="/Users/naghamkhourifarah/Documents/Annotation1.rds")
sig_genes <- FindAllMarkers(ect_all, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
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

ectmes <- merge(ect_all,mes_all)
Idents(cds) <- "cell_type"
for (i in levels(ectmes)) {
  id <- rownames(ectmes@meta.data[ectmes@meta.data$cellType4 == i,])
  cds <- SetIdent(cds,id,value = i)
}        
cds$ectmes <- Idents(cds)
DimPlot(cds, group.by = "ectmes")
saveRDS(cds,"data/cds.rds")
saveRDS(ect_all,"data/ect_all.rds")
saveRDS(mes_all,"data/mes_all.rds")
saveRDS(E11,"data/E11.rds")
saveRDS(ect,"data/E11_ect.rds")
saveRDS(mes,"data/E11_mes.rds")


