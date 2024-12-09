library(Seurat)
library(ShinyCell)
library(stringr)

cds1 <- readRDS("cds_face_human_temp.rds")
columns <- colnames(cds1@meta.data)
columns <- gsub("tempanno", "annotation", columns)
colnames(cds1@meta.data) <- columns
cds2 <- cds1
Idents(cds1) <- "cell_type1"
Idents(cds2) <- "annotation"
cds_mes <- readRDS("mes_temp1.rds")
columns <- colnames(cds_mes@meta.data)
columns <- gsub("MtoHAlltoALL", "annotation", columns)
colnames(cds_mes@meta.data) <- columns
Idents(cds_mes) <- "annotation"
cds_ect <- readRDS("ect_final.rds")
columns <- colnames(cds_ect@meta.data)
columns <- gsub("tempanno", "annotation", columns)
colnames(cds_ect@meta.data) <- columns
Idents(cds_ect) <- "annotation"
cds_cncc <- readRDS("cncc_final.rds")
columns <- colnames(cds_cncc@meta.data)
columns <- gsub("tempanno", "annotation", columns)
colnames(cds_cncc@meta.data) <- columns
Idents(cds_cncc) <- "annotation"

mcds1 <- readRDS("face_mouse_final.rds")
mcds2 <- mcds1
Idents(mcds1) <- "cell_type1"
Idents(mcds2) <- "annotation"
cds_mcncc <- readRDS("mcncc_final.rds")
columns <- colnames(cds_mcncc@meta.data)
columns <- gsub("tempanno", "annotation", columns)
colnames(cds_mcncc@meta.data) <- columns
Idents(cds_mcncc) <- "annotation"
cds_mmes <- readRDS("mmes_final.rds")
Idents(cds_mmes) <- "annotation"
cds_mect <- readRDS("mect_final.rds")
Idents(cds_mect) <- "annotation"

main_conf <- createConfig(cds1, legendCols = 10, maxLevels = 100)
subtypes_conf <- createConfig(cds2, legendCols = 10, maxLevels = 100)
mes_conf <- createConfig(cds_mes, legendCols = 10, maxLevels = 100)
ect_conf <- createConfig(cds_ect, legendCols = 10, maxLevels = 100)
cncc_conf <- createConfig(cds_cncc, legendCols = 10, maxLevels = 100)


main_conf = modDefault(main_conf, "cell_type1", "stage")
subtypes_conf = modDefault(subtypes_conf, "annotation", "stage")
mes_conf = modDefault(mes_conf, "annotation", "stage")
ect_conf = modDefault(ect_conf, "annotation", "stage")
cncc_conf = modDefault(cncc_conf, "annotation", "stage")

m_main_conf <- createConfig(mcds1, legendCols = 10, maxLevels = 100)
m_subtypes_conf <- createConfig(mcds2, legendCols = 10, maxLevels = 100)
mmes_conf <- createConfig(cds_mmes, legendCols = 10, maxLevels = 100)
mect_conf <- createConfig(cds_mect, legendCols = 10, maxLevels = 100)
mcncc_conf <- createConfig(cds_mcncc, legendCols = 10, maxLevels = 100)


m_main_conf = modDefault(m_main_conf, "cell_type1", "stage")
m_subtypes_conf = modDefault(m_subtypes_conf, "annotation", "stage")
mmes_conf = modDefault(mmes_conf, "annotation", "stage")
mect_conf = modDefault(mect_conf, "annotation", "stage")
mcncc_conf = modDefault(mcncc_conf, "annotation", "stage")

prioritized_genes <- c("ABCA1","ACTB","AGO1","ALX4","ATXN1","BCL3","CACNA1D","CACNA1G","CDH1","CLIP1","COL11A1","COL2A1","CRISPLD2","CTNND1","DGKH","DOCK2","EPB41","EPHB3","ERF","EYA1","FGFR1","GRHL3","HMCN1","HSPG2","IFT122","INTS1","IRF6","KCNQ5","KIF7","MAST4","MEIS2","MSI1","NEDD4L","NIPBL","PRTG","PTCH1","SACS","SATB2","SH3PXD2B","SMAD4","TBX3","TFAP2A","TGFBR2","TMEM132D","TNS3","TRIM27","WDR35","ZFHX4")
prioritized_genes_mouse <- str_to_title(prioritized_genes)

makeShinyFiles(cds1, main_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc1",
               shiny.dir = "shinyMulti",
               default.gene1 = "MSX1", default.gene2 = "TP63",
               default.multigene = prioritized_genes)

makeShinyFiles(cds2, subtypes_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc2",
               shiny.dir = "shinyMulti",
               default.gene1 = "MSX1", default.gene2 = "TP63",
               default.multigene = prioritized_genes)

makeShinyFiles(cds_mes, mes_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc3",
               shiny.dir = "shinyMulti",
               default.gene1 = "TWIST1", default.gene2 = "PRRX1",
               default.multigene = prioritized_genes)

makeShinyFiles(cds_ect, ect_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc4",
               shiny.dir = "shinyMulti",
               default.gene1 = "TP63", default.gene2 = "IRF6",
               default.multigene = prioritized_genes)

makeShinyFiles(cds_cncc, cncc_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc5",
               shiny.dir = "shinyMulti",
               default.gene1 = "SOX10", default.gene2 = "TFAP2A",
               default.multigene = prioritized_genes)

makeShinyFiles(mcds1, m_main_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc6",
               shiny.dir = "shinyMulti",
               default.gene1 = "Msx1", default.gene2 = "Trp63",
               default.multigene = prioritized_genes_mouse)

makeShinyFiles(mcds2, m_subtypes_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc7",
               shiny.dir = "shinyMulti",
               default.gene1 = "Msx1", default.gene2 = "Trp63",
               default.multigene = prioritized_genes_mouse)

makeShinyFiles(cds_mmes, mmes_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc8",
               shiny.dir = "shinyMulti",
               default.gene1 = "Twist1", default.gene2 = "Prrx1",
               default.multigene = prioritized_genes_mouse)

makeShinyFiles(cds_mect, mect_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc9",
               shiny.dir = "shinyMulti",
               default.gene1 = "Trp63", default.gene2 = "Irf6",
               default.multigene = prioritized_genes_mouse)

makeShinyFiles(cds_mcncc, mcncc_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc10",
               shiny.dir = "shinyMulti",
               default.gene1 = "Sox10", default.gene2 = "Tfap2a",
               default.multigene = prioritized_genes_mouse)



citation = list(
  author  = "Farah et al.",
  title   = "Dynamics of cell-type level gene expression during patterning of the human face",
  journal = "biorxiv.org",
  volume  = "",
  page    = "",
  year    = "2024", 
  doi     = "TBD",
  link    = "TBD")

makeShinyCodesMulti(
  shiny.title = "Dynamics of cell-type level gene expression during patterning of human and mouse faces", shiny.footnotes = citation,
  shiny.prefix = c("sc1","sc2","sc3","sc4","sc5","sc6","sc7","sc8","sc9","sc10"),
  defPtSiz = c(0.5,0.5,0.5,0.5,0.5),
  shiny.headers = c("Human Main Types", "Human Subtypes", "Human Mesenchymal Subclustering", "Human Epithelial Subclustering", "Human CNCC Subclustering", "Mouse Main Types", "Mouse Subtypes", "Mouse Mesenchymal Subclustering", "Mouse Epithelial Subclustering", "Mouse CNCC Subclustering"),
  shiny.dir = "shinyMulti")






