library(sceasy)
library(reticulate)
library(SingleCellExperiment)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(scater)

setwd("/scr1/users/manchela/Data")

#change to ENSEMBL ID for cellxgene

# convert rds to h5ad for human objects
seu <- readRDS("human_face_no-neuro_clustering_cellrangerARC-raw_emptyDrops_singlets_finalannot_27Mar.rds")
geneid <- mapIds(org.Hs.eg.db, keys=rownames(seu), keytype="SYMBOL", column="ENSEMBL")
rownames(seu) <- uniquifyFeatureNames(rownames(seu), geneid)
seu[["RNA"]] <- as(seu[["RNA"]],Class = "Assay")
human_main_adata <- sceasy:::seurat2anndata(seu,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                                            outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_human_final.h5ad")


seu_mes <- readRDS("human_face_mes_cellrangerARC-raw_emptyDrops_singlets_8Apr.rds")
geneid <- mapIds(org.Hs.eg.db, keys=rownames(seu_mes), keytype="SYMBOL", column="ENSEMBL")
rownames(seu_mes) <- uniquifyFeatureNames(rownames(seu_mes), geneid)
seu_mes[["RNA"]] <- as(seu_mes[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(seu_mes,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_human_mes_final.h5ad")


seu_ect <- readRDS("human_face_ect_cellrangerARC-raw_emptyDrops_singlets_27Mar.rds")
geneid <- mapIds(org.Hs.eg.db, keys=rownames(seu_ect), keytype="SYMBOL", column="ENSEMBL")
rownames(seu_ect) <- uniquifyFeatureNames(rownames(seu_ect), geneid)
seu_ect[["RNA"]] <- as(seu_ect[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(seu_ect,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_human_ect_final.h5ad")


seu_cncc <- readRDS("human_face_cncc_cellrangerARC-raw_emptyDrops_singlets_27Mar.rds")
geneid <- mapIds(org.Hs.eg.db, keys=rownames(seu_cncc), keytype="SYMBOL", column="ENSEMBL")
rownames(seu_cncc) <- uniquifyFeatureNames(rownames(seu_cncc), geneid)
seu_cncc[["RNA"]] <- as(seu_cncc[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(seu_cncc,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_human_cncc_final.h5ad")


seu_prog <- readRDS("human_face_prog_cellrangerARC-raw_emptyDrops_singlets_27Mar.rds")
geneid <- mapIds(org.Hs.eg.db, keys=rownames(seu_prog), keytype="SYMBOL", column="ENSEMBL")
rownames(seu_prog) <- uniquifyFeatureNames(rownames(seu_prog), geneid)
seu_prog[["RNA"]] <- as(seu_prog[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(seu_prog,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_human_prog_final.h5ad")


# convert rds to h5ad for mouse objects

mcds1 <- readRDS("face_mouse_final.rds")
geneid <- mapIds(org.Mm.eg.db, keys=rownames(mcds1), keytype="SYMBOL", column="ENSEMBL")
rownames(mcds1) <- uniquifyFeatureNames(rownames(mcds1), geneid)
mcds1[["RNA"]] <- as(mcds1[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(mcds1,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                                            outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_mouse_final.h5ad")


cds_mcncc <- readRDS("mcncc_final.rds")
geneid <- mapIds(org.Mm.eg.db, keys=rownames(cds_mcncc), keytype="SYMBOL", column="ENSEMBL")
rownames(cds_mcncc) <- uniquifyFeatureNames(rownames(cds_mcncc), geneid)
cds_mcncc[["RNA"]] <- as(cds_mcncc[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(cds_mcncc,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_mouse_cncc_final.h5ad")


cds_mmes <- readRDS("mmes_final.rds")
geneid <- mapIds(org.Mm.eg.db, keys=rownames(cds_mmes), keytype="SYMBOL", column="ENSEMBL")
rownames(cds_mmes) <- uniquifyFeatureNames(rownames(cds_mmes), geneid)
cds_mmes[["RNA"]] <- as(cds_mmes[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(cds_mmes,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_mouse_mes_final.h5ad")


cds_mect <- readRDS("mect_final.rds")
geneid <- mapIds(org.Mm.eg.db, keys=rownames(cds_mect), keytype="SYMBOL", column="ENSEMBL")
rownames(cds_mect) <- uniquifyFeatureNames(rownames(cds_mect), geneid)
cds_mect[["RNA"]] <- as(cds_mect[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(cds_mect,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_mouse_ect_final.h5ad")


seu_ect <- readRDS("human_face_ect_cellrangerARC-raw_emptyDrops_singlets_27Mar.rds")
geneid <- mapIds(org.Hs.eg.db, keys=rownames(seu_ect), keytype="SYMBOL", column="ENSEMBL")
rownames(seu_ect) <- uniquifyFeatureNames(rownames(seu_ect), geneid)
seu_ect[["RNA"]] <- as(seu_ect[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(seu_ect,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_human_ect_final.h5ad")


seu_cncc <- readRDS("human_face_cncc_cellrangerARC-raw_emptyDrops_singlets_27Mar.rds")
geneid <- mapIds(org.Hs.eg.db, keys=rownames(seu_cncc), keytype="SYMBOL", column="ENSEMBL")
rownames(seu_cncc) <- uniquifyFeatureNames(rownames(seu_cncc), geneid)
seu_cncc[["RNA"]] <- as(seu_cncc[["RNA"]],Class = "Assay")
sceasy:::seurat2anndata(seu_cncc,transfer_layers='counts', assay='RNA',drop_single_values=FALSE,
                        outFile="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/cellxgene/face_human_cncc_final.h5ad")

