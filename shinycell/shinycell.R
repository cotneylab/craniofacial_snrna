library(Seurat)
library(ShinyCell2)
library(stringr)
library(grid)
library(scales)
library(scCustomize)
library(SeuratDisk)

# Code to generate Shiny App

setwd("/scr1/users/manchela/Data")

seu <- readRDS("human_face_no-neuro_clustering_cellrangerARC-raw_emptyDrops_singlets_finalannot_27Mar.rds")
seu$subtype_reduced[which(seu$subtype_reduced=="other blood cells")] <- "immune"
seu$subtype_reduced[which(seu$subtype_reduced=="erythrocytes")] <- "erythrocytes"
Idents(seu) <- seu$celltype

main_conf <- createConfig(seu, legendCols = 10, maxLevels = 100)
main_conf = modColours(main_conf, meta.to.mod = "celltype", 
                       new.colours= hue_pal()(8))
# grid.draw(showLegend(main_conf))

main_conf <- delMeta(main_conf, c("old.ident","seurat_clusters","RBC_ModuleScore1","WNT_ModuleScore1",
                                  "EarlyMigration_ModuleScore1","Delam_ModuleScore1","Delam_ModuleScore1",
                                  "NeuralTube_ModuleScore1","initial_cluster","colors_celltype","colors_subtype",
                                  "monocle3_pseudotime","Pre.Delam_ModuleScore1","data","rep",
                                  setdiff(grep("RNA_snn_res", main_conf$ID,value=T),"RNA_snn_res.0.2"))) 


seu2 <- seu
Idents(seu2) <- seu2$subtype
subtypes_conf <- createConfig(seu2, legendCols = 10, maxLevels = 100)

seu_mes <- readRDS("human_face_mes_cellrangerARC-raw_emptyDrops_singlets_8Apr.rds")
Idents(seu_mes) <- seu_mes$subtype_final
mes_conf <- createConfig(seu_mes, legendCols = 10, maxLevels = 100)
mes_colors <- DiscretePalette_scCustomize(num_colors = 24,palette = "stepped")
mes_colors <- c(mes_colors,"#E47D00")
mes_conf = modColours(mes_conf, meta.to.mod = "subtype_final", 
                       new.colours= mes_colors)

mes_conf <- delMeta(mes_conf, c("old.ident","seurat_clusters","RBC_ModuleScore1","WNT_ModuleScore1",
                                  "EarlyMigration_ModuleScore1","Delam_ModuleScore1","Delam_ModuleScore1",
                                  "NeuralTube_ModuleScore1","Pre.Delam_ModuleScore1","data","rep",
                                  setdiff(grep("RNA_snn_res", mes_conf$ID,value=T),"RNA_snn_res.0.8"),
                                grep("Pluripotency_ModuleScore",mes_conf$ID,value=T),"initial_cluster",
                                "predicted","color","S.Score.adjust","G2M.Score.adjust","cellcycle.ratio")) 


seu_ect <- readRDS("human_face_ect_cellrangerARC-raw_emptyDrops_singlets_27Mar.rds")
Idents(seu_ect) <- seu_ect$subtype_final
ect_conf <- createConfig(seu_ect, legendCols = 10, maxLevels = 100)
ect_colors <- DiscretePalette_scCustomize(num_colors = length(unique(as.character(seu_ect$subtype_final))),
                                          palette = "alphabet")
ect_conf = modColours(ect_conf, meta.to.mod = "subtype_final", 
                      new.colours= ect_colors)

ect_conf <- delMeta(ect_conf, c("old.ident","seurat_clusters","data","rep","subtype1","predicted.score","subtype2",
                                setdiff(grep("RNA_snn_res", ect_conf$ID,value=T),"RNA_snn_res.0.8"),"predicted.score1",
                                "predicted","color","S.Score.adjust","G2M.Score.adjust","cellcycle.ratio")) 



seu_cncc <- readRDS("human_face_cncc_cellrangerARC-raw_emptyDrops_singlets_27Mar.rds")
Idents(seu_cncc) <- seu_cncc$subtype_final
cncc_conf <- createConfig(seu_cncc, legendCols = 10, maxLevels = 100)
cncc_colors <- DiscretePalette_scCustomize(num_colors = length(unique(as.character(seu_cncc$RNA_snn_res.1))),
                                          palette = "ditto_seq")
cncc_conf = modColours(cncc_conf, meta.to.mod = "subtype_final", 
                      new.colours= cncc_colors)

cncc_conf <- delMeta(cncc_conf, c("old.ident","seurat_clusters","data","rep","subtype1","predicted.score","subtype2",
                                setdiff(grep("RNA_snn_res", cncc_conf$ID,value=T),"RNA_snn_res.1"),"predicted.score1",
                                "predicted","color","S.Score.adjust","G2M.Score.adjust","cellcycle.ratio",
                                grep("HOX",cncc_conf$ID,value=T),grep("zebrafish",cncc_conf$ID,value=T),
                                grep("Solfatov",cncc_conf$ID,value=T),"Schwann_allorgan_ModuleScore1","Placode_ModuleScore1",
                                "EarlyMigration_ModuleScore1","Delam_ModuleScore1","Delam_ModuleScore1",
                                "NeuralTube_ModuleScore1","Pre.Delam_ModuleScore1","Migration1_ModuleScore1","Migration2_ModuleScore1",
                                "subtype3","subtype_stage","subtype_function","subtype_final2")) 


seu_prog <- readRDS("human_face_prog_cellrangerARC-raw_emptyDrops_singlets_27Mar.rds")
Idents(seu_prog) <- seu_prog$subtype_final
prog_conf <- createConfig(seu_prog, legendCols = 10, maxLevels = 100)
prog_conf <- delMeta(prog_conf, c("old.ident","seurat_clusters","data","rep","subtype1","colors_subtype",
                                  setdiff(grep("RNA_snn_res", prog_conf$ID,value=T),"RNA_snn_res.1"),"color",
                                  "S.Score.adjust","G2M.Score.adjust","cellcycle.ratio","subtype","colors_celltype",
                                  "monocle3_pseudotime","RBC_ModuleScore1","WNT_ModuleScore1","initial_cluster",
                                  "EarlyMigration_ModuleScore1","Pre.Delam_ModuleScore1","Delam_ModuleScore1",
                                  "Germ_ModuleScore1","Hematopoietic_ModuleScore1","Chondrocyte_ModuleScore1",
                                  "MyocyteProg_ModuleScore1","Lymphoblast_ModuleScore1","Osteoblast_ModuleScore1",
                                  "Fibroblast_ModuleScore1","Hypoblast_ModuleScore1","Mes_ModuleScore1",
                                  "OsteoblastProg_ModuleScore1","RetinalProg_ModuleScore1","NeuralTube1_ModuleScore1",
                                  "NeuralTube2_ModuleScore1","NeuralTube3_ModuleScore1","IntMesoderm_ModuleScore1",
                                  "LateralPlateMesoderm1_ModuleScore1","LateralPlateMesoderm2_ModuleScore1",
                                  "ParaxialMesoderm_ModuleScore1","Mesoderm_ModuleScore1","NMP_ModuleScore1",
                                  "NT_ModuleScore1","NeuralTube_ModuleScore1","Endoderm_ModuleScore1","Notochord_ModuleScore1",
                                  "Erythroid_ModuleScore1","PresomiticMesoderm_ModuleScore1","Somite_ModuleScore1",
                                  "sctype_classification","Astrocyte_ModuleScore1","Microglia_ModuleScore1",
                                  "Endothelium_ModuleScore1")) 



mcds1 <- readRDS("face_mouse_final.rds")
mcds1$celltype <- mcds1$cell_type1
Idents(mcds1) <- mcds1$celltype
m_main_conf <- createConfig(mcds1, legendCols = 10, maxLevels = 100)
m_main_conf = modColours(m_main_conf, meta.to.mod = "celltype", 
                       new.colours= hue_pal()(8)[1:7])
m_main_conf <- delMeta(m_main_conf, c("rep","CNCC_ModuleScore1","seurat_clusters","cell_type1","cell_type",
                                  setdiff(grep("RNA_snn_res", m_main_conf$ID,value=T),"RNA_snn_res.0.1"),"annotation")) 

mcds2 <- mcds1
Idents(mcds2) <- mcds2$subtype
m_subtypes_conf <- createConfig(mcds2, legendCols = 10, maxLevels = 100)
m_subtypes_conf <- delMeta(m_subtypes_conf, c("rep","CNCC_ModuleScore1","seurat_clusters","cell_type1","cell_type",
                                      setdiff(grep("RNA_snn_res", m_subtypes_conf$ID,value=T),"RNA_snn_res.0.1"),"annotation")) 

cds_mcncc <- readRDS("mcncc_final.rds")
Idents(cds_mcncc) <- cds_mcncc$subtype
cds_mcncc$celltype <- cds_mcncc$cell_type1
mcncc_conf <- createConfig(cds_mcncc, legendCols = 10, maxLevels = 100)
mcncc_conf = modColours(mcncc_conf, meta.to.mod = "subtype", 
                         new.colours= hue_pal()(length(unique(cds_mcncc$subtype))))
mcncc_conf <- delMeta(mcncc_conf, c("rep","seurat_clusters","cell_type1","cell_type","cellType","tempanno",
                                    setdiff(grep("RNA_snn_res", mcncc_conf$ID,value=T),"RNA_snn_res.1"))) 


cds_mmes <- readRDS("mmes_final.rds")
Idents(cds_mmes) <- cds_mmes$subtype
cds_mmes$celltype <- cds_mmes$cell_type1
mmes_conf <- createConfig(cds_mmes, legendCols = 10, maxLevels = 100)
mmes_conf = modColours(mmes_conf, meta.to.mod = "subtype", 
                        new.colours= hue_pal()(length(unique(cds_mmes$subtype))))
mmes_conf <- delMeta(mmes_conf, c("rep","seurat_clusters","cell_type1","cell_type","cellType",
                                    setdiff(grep("RNA_snn_res", mmes_conf$ID,value=T),"RNA_snn_res.0.8"))) 


cds_mect <- readRDS("mect_final.rds")
Idents(cds_mect) <- cds_mect$subtype
cds_mect$celltype <- cds_mect$cell_type1
mect_conf <- createConfig(cds_mect, legendCols = 10, maxLevels = 100)
mect_conf = modColours(mect_conf, meta.to.mod = "subtype", 
                       new.colours= hue_pal()(length(unique(cds_mect$subtype))))
mect_conf <- delMeta(mect_conf, c("rep","seurat_clusters","cell_type1","cell_type","cellType",
                                  setdiff(grep("RNA_snn_res", mect_conf$ID,value=T),"RNA_snn_res.0.8"))) 



main_conf = modDefault(main_conf, "celltype", "stage")
subtypes_conf = modDefault(subtypes_conf, "subtype", "stage")
mes_conf = modDefault(mes_conf, "subtype_final", "stage")
ect_conf = modDefault(ect_conf, "subtype_final", "stage")
cncc_conf = modDefault(cncc_conf, "subtype_final", "stage")
prog_conf = modDefault(prog_conf, "subtype_final", "stage")


m_main_conf = modDefault(m_main_conf, "celltype", "stage")
m_subtypes_conf = modDefault(m_subtypes_conf, "subtype", "stage")
mmes_conf = modDefault(mmes_conf, "subtype", "stage")
mect_conf = modDefault(mect_conf, "subtype", "stage")
mcncc_conf = modDefault(mcncc_conf, "subtype", "stage")

prioritized_genes <- c("ABCA1","ACTB","AGO1","ALX4","ATXN1","BCL3","CACNA1D","CACNA1G","CDH1","CLIP1","COL11A1","COL2A1","CRISPLD2","CTNND1","DGKH","DOCK2","EPB41","EPHB3","ERF","EYA1","FGFR1","GRHL3","HMCN1","HSPG2","IFT122","INTS1","IRF6","KCNQ5","KIF7","MAST4","MEIS2","MSI1","NEDD4L","NIPBL","PRTG","PTCH1","SACS","SATB2","SH3PXD2B","SMAD4","TBX3","TFAP2A","TGFBR2","TMEM132D","TNS3","TRIM27","WDR35","ZFHX4")
prioritized_genes_mouse <- str_to_title(prioritized_genes)

makeShinyFiles(seu, main_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc1",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "MSX1", default.gene2 = "TP63",
               default.multigene = prioritized_genes)

makeShinyFiles(seu2, subtypes_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc2",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "MSX1", default.gene2 = "TP63",
               default.multigene = prioritized_genes)

makeShinyFiles(seu_mes, mes_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc3",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "TWIST1", default.gene2 = "PRRX1",
               default.multigene = prioritized_genes)

makeShinyFiles(seu_ect, ect_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc4",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "TP63", default.gene2 = "IRF6",
               default.multigene = prioritized_genes)

makeShinyFiles(seu_cncc, cncc_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc5",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "SOX10", default.gene2 = "TFAP2A",
               default.multigene = prioritized_genes)

makeShinyFiles(seu_prog, prog_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc11",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "SOX2", default.gene2 = "PAX6",
               default.multigene = prioritized_genes)

makeShinyFiles(mcds1, m_main_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc6",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "Msx1", default.gene2 = "Trp63",
               default.multigene = prioritized_genes_mouse)

makeShinyFiles(mcds2, m_subtypes_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc7",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "Msx1", default.gene2 = "Trp63",
               default.multigene = prioritized_genes_mouse)

makeShinyFiles(cds_mmes, mmes_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc8",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "Twist1", default.gene2 = "Prrx1",
               default.multigene = prioritized_genes_mouse)

makeShinyFiles(cds_mect, mect_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc9",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "Trp63", default.gene2 = "Irf6",
               default.multigene = prioritized_genes_mouse)

makeShinyFiles(cds_mcncc, mcncc_conf, gene.mapping = TRUE,
               gex.assay = "RNA", gex.slot = "data",
               shiny.prefix = "sc10",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "Sox10", default.gene2 = "Tfap2a",
               default.multigene = prioritized_genes_mouse)

# spatial objects
human_spatial <- readRDS("human_emb.merge_final.rds")
human_spatial$celltype <- human_spatial@active.ident
DefaultAssay(human_spatial) <- "SCT"

cds <- readRDS("human_face_no-neuro_clustering_cellrangerARC-raw_emptyDrops_singlets_finalannot_27Mar.rds")

human_spatial_image1 <- subset(human_spatial, cells = colnames(human_spatial)[which(human_spatial$orig.ident=="s1_spatial")])
colnames(human_spatial_image1@meta.data)[na.omit(match(unique(cds$subtype),colnames(human_spatial_image1@meta.data)))] <- paste0(colnames(human_spatial_image1@meta.data)[na.omit(match(unique(cds$subtype),colnames(human_spatial_image1@meta.data)))],"_ModuleScore")

scConf_hspatial1 <- createConfig(human_spatial_image1, legendCols = 10, maxLevels = 100)
scConf_hspatial1 <- delMeta(scConf_hspatial1, c("slice", "region","SCT_snn_res.0.8","seurat_clusters","CNCC_Full_ModuleScore1",
                                     "YankeeCNCC_ModuleScore1","TreacherCollins_ModuleScore1","ALX_ModuleScore1","orig.ident"))
scConf_hspatial1 = modDefault(scConf_hspatial1, "celltype","nCount_SCT")


makeShinyFilesGEX(human_spatial_image1, scConf_hspatial1, gex.assay = "SCT", gex.slot = "data",
                  shiny.prefix = "sc12", shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell/",
                  default.gene1 = "MSX1", default.gene2 = "TP63", 
                  default.multigene = prioritized_genes
                  )

sc1image = list()
sc1image$coord <- GetTissueCoordinates(human_spatial_image1)
colnames(sc1image$coord)[1:2] <- c("imagerow","imagecol")
bg_image <- GetImage(human_spatial_image1, mode = "raster")
sc1image$bg_image <- bg_image
sc1image$lowres <- human_spatial_image1@images[[1]]@scale.factors$lowres
saveRDS(sc1image, file = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell/sc12image.rds")


human_spatial_image2 <- subset(human_spatial, cells = colnames(human_spatial)[which(human_spatial$orig.ident=="s2_spatial")])

colnames(human_spatial_image2@meta.data)[na.omit(match(unique(cds$subtype),colnames(human_spatial_image2@meta.data)))] <- paste0(colnames(human_spatial_image2@meta.data)[na.omit(match(unique(cds$subtype),colnames(human_spatial_image2@meta.data)))],"_ModuleScore")

scConf_hspatial2 <- createConfig(human_spatial_image2, legendCols = 10, maxLevels = 100)
scConf_hspatial2 <- delMeta(scConf_hspatial2, c("slice", "region","SCT_snn_res.0.8","seurat_clusters","CNCC_Full_ModuleScore1",
                                                "YankeeCNCC_ModuleScore1","TreacherCollins_ModuleScore1","ALX_ModuleScore1","orig.ident"))
scConf_hspatial2 = modDefault(scConf_hspatial2, "celltype","nCount_SCT")

makeShinyFilesGEX(human_spatial_image2, scConf_hspatial2, gex.assay = "SCT", gex.slot = "data",
                  shiny.prefix = "sc13", shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell/",
                  default.gene1 = "MSX1", default.gene2 = "TP63", 
                  default.multigene = prioritized_genes
)

sc2image = list()
sc2image$coord <- GetTissueCoordinates(human_spatial_image2)
colnames(sc2image$coord)[1:2] <- c("imagerow","imagecol")
bg_image <- GetImage(human_spatial_image2, mode = "raster")
sc2image$bg_image <- bg_image
sc2image$lowres <- human_spatial_image2@images[[1]]@scale.factors$lowres
saveRDS(sc2image, file = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell/sc13image.rds")



dsouza_mouse_spatial <- readRDS("DSOUZA_anno_spatial_E15_final.rds")
dsouza_mouse_spatial$subtype <- dsouza_mouse_spatial$subtype_cca
dsouza_mouse_spatial$celltype <- dsouza_mouse_spatial$celltype_cca

scConf_mspatial <- createConfig(dsouza_mouse_spatial, legendCols = 10, maxLevels = 100)
scConf_mspatial <- delMeta(scConf_mspatial, c("orig.ident", "centroid_1","centroid_2","seurat_clusters","Noharmony",
                                                "tempanno","cell_type","subtype_cca","celltype_cca","subtype_cca_toplot"))
scConf_mspatial = modDefault(scConf_mspatial, "celltype","subtype")

dsouza_mouse_spatial <- JoinLayers(dsouza_mouse_spatial)
makeShinyFiles(dsouza_mouse_spatial, scConf_mspatial, gene.mapping = TRUE,
               gex.assay = "Spatial", gex.slot = "data",
               shiny.prefix = "sc14",
               shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell",
               default.gene1 = "Msx1", default.gene2 = "Trp63",
               default.multigene = prioritized_genes_mouse)



citation = list(
  author  = "Khouri-Farah, N., Manchel, A., Winchester, E.W., Schilder B.M., Robinson, K., Curtis, S.W., Skene, N.G., Leslie-Clarkson, E.J., Cotney, J.",
  title   = "Gene expression patterns of the developing human face at single cell resolution reveal cell type contributions to normal facial variation and disease risk",
  journal = "biorxiv.org",
  volume  = "",
  page    = "",
  year    = "2025", 
  doi     = "10.1101/2025.01.18.633396",
  link    = "https://www.biorxiv.org/content/10.1101/2025.01.18.633396v2.full")

makeShinyCodes(
  shiny.title = "Dynamics of cell-type level gene expression during patterning of human and mouse faces", shiny.footnotes = citation,
  shiny.prefix = c("sc1","sc2","sc3","sc4","sc5","sc11","sc6","sc7","sc8","sc9","sc10","sc12","sc13","sc14"),
  defPtSiz = c(0.5,0.5,0.5,0.5,0.5),
  shiny.headers = c("Human Main Types", "Human Subtypes", "Human Mesenchymal Subclustering", "Human Ectodermal Subclustering", "Human CNCC Subclustering", "Human Progenitor Subclustering",
                    "Mouse Main Types", "Mouse Subtypes", "Mouse Mesenchymal Subclustering", "Mouse Ectodermal Subclustering", "Mouse CNCC Subclustering", "Human Spatial Slice1","Human Spatial Slice2",
                    "Mouse Spatial"),
  shiny.dir = "/scr1/users/manchela/Data/Multiome/Humanface/GEX/shinyCell")



## Add these lines of code to the shinyFunc.R script ##

# #new code
# img_with_alpha <- apply(inpImg[["bg_image"]], 2, adjustcolor, alpha.f = 0.1)
# inpImgGrob <- rasterGrob(img_with_alpha, interpolate = FALSE, 
#                          width=unit(1,"npc"), height=unit(1,"npc")) 
# #end of new code
# 
# 
# #new code
# inpsiz = inpsiz*2
# #end of new code