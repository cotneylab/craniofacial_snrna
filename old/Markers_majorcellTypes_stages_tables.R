library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
#subset E11 matching TW E11 data
cds <- readRDS("data/cds.rds")
Idents(cds) <- "cell_type"
dir.create("tables/E10_E15_majorcellTypes_stages")

for (i in levels(cds)) {
  seu_temp <- subset(cds, idents = i)
  Idents(seu_temp) <- "stage"
  sig_genes <- FindAllMarkers(seu_temp, assay = "RNA" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                              min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                              only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                              latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                              pseudocount.use = 1, return.thresh = 0.01)
  sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
  openxlsx::write.xlsx(sig_genes,paste0("tables/E10_E15_majorcellTypes_stages/SigGenes_",i,".xlsx"))
                       
}
