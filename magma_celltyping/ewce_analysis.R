library(EWCE)
library(SingleCellExperiment)
library(scRNAseq)
library(scater)
library(Seurat)
library(tidyr)
library(pheatmap)
# Sys.setenv('R_MAX_VSIZE'=64000000000)
options(future.globals.maxSize = 10 * 1024^3)
setwd("/scr1/users/manchela/Data")

cds_subtype <- readRDS("/scr1/users/manchela/Data/human_face_no-neuro_clustering_cellrangerARC-raw_emptyDrops_singlets_finalannot_27Mar.rds")

## for figure 8A

#normalize data for use in EWCE
# cds_subtype <- PercentageFeatureSet(cds_subtype, pattern = "^MT-", col.name = "percent.mt")
cds_subtype <- SCTransform(cds_subtype, vars.to.regress = "percent.mt", verbose = FALSE)

# saveRDS(cds_subtype,"/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/cds_subtype.rds")
cds_subtype <- readRDS("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/cds_subtype.rds")
cds_subtype$subtype_noprog <- as.character(cds_subtype$subtype)
cds_subtype$subtype_noprog[which(cds_subtype$celltype=="SOX2+_SOX10-")] <- "SOX2+_SOX10-"
cds_subtype$subtype_noprog <- factor(cds_subtype$subtype_noprog,levels=c(levels(cds_subtype$subtype)[1:70],"SOX2+_SOX10-"))

sce_subtype <- as.SingleCellExperiment(cds_subtype, assay = 'SCT')

annotLevels = list(level1class=sce_subtype$subtype,
                   level2class=sce_subtype$celltype,
                   level3class=sce_subtype$subtype_noprog)

Face_Subtype_SCtransform_CTD <- generate_celltype_data(exp=assays(sce_subtype)$counts,
                                                       annotLevels=annotLevels,
                                                       specificity_quantiles = TRUE,
                                                       groupName="face_subtype_sctransform_ctd",
                                                       no_cores = 10,
                                                       input_species = "human",
                                                       output_species = "human",
                                                       return_ctd = TRUE,
                                                       savePath="/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/")


face_subtype_ctd <- EWCE::load_rdata("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/ctd_face_subtype_sctransform_ctd.rda") 

setwd("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS")
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

gmfk_all <- read.table("leslie/20241024_CPSeqGMFKDDD_allDNs.txt", header=TRUE, sep ="\t")
gmfk <- subset(gmfk_all, study == "GMKF_CPSeq")
gmfk_genes <- gmfk$Gene
gmfk_lof_genes <- subset(gmfk, variant.class == "lof")$Gene
gmfk_mis_genes <- subset(gmfk, variant.class == "mis")$Gene
gmfk_syn_genes <- subset(gmfk, variant.class == "syn")$Gene
gmfk_nosyn_genes <- subset(gmfk, variant.class != "syn")$Gene
ddd<- subset(gmfk_all, study == "DDD")
ddd_genes <- ddd$Gene
ddd_lof_genes <- subset(ddd, variant.class == "lof")$Gene
ddd_mis_genes <- subset(ddd, variant.class == "mis")$Gene
ddd_syn_genes <- subset(ddd, variant.class == "syn")$Gene
ddd_nosyn_genes <- subset(ddd, variant.class != "syn")$Gene

ofc_panels <- read.table("leslie/20241028_OFC_GenePanelList.txt", header=TRUE, sep ="\t")
ofc_panel_genes <- ofc_panels$approved.symbol

novel_gene_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                      sctSpecies = "human",
                                                      genelistSpecies = "human",
                                                      hits = novel_genes, 
                                                      reps = 10000,
                                                      geneSizeControl = TRUE,
                                                      annotLevel = 2)

cleftgenedb_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                       sctSpecies = "human",
                                                       genelistSpecies = "human",
                                                       hits = cleftgenedb, 
                                                       reps = 10000,
                                                       geneSizeControl = TRUE,
                                                       annotLevel = 2)

black_hubs_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                      sctSpecies = "human",
                                                      genelistSpecies = "human",
                                                      hits = black_hubs, 
                                                      reps = 10000,
                                                      geneSizeControl = TRUE,
                                                      annotLevel = 2)

black_module_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                        sctSpecies = "human",
                                                        genelistSpecies = "human",
                                                        hits = black_genes, 
                                                        reps = 10000,
                                                        geneSizeControl = TRUE,
                                                        annotLevel = 2)

gnomad_dec1_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                       sctSpecies = "human",
                                                       genelistSpecies = "human",
                                                       hits = gnomad_dec1_genes, 
                                                       reps = 10000,
                                                       geneSizeControl = TRUE,
                                                       annotLevel = 2)

gnomad_dec9_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                       sctSpecies = "human",
                                                       genelistSpecies = "human",
                                                       hits = gnomad_dec9_genes, 
                                                       reps = 10000,
                                                       geneSizeControl = TRUE,
                                                       annotLevel = 2)

gnomad_pli_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                      sctSpecies = "human",
                                                      genelistSpecies = "human",
                                                      hits = gnomad_pli_genes, 
                                                      reps = 10000,
                                                      geneSizeControl = TRUE,
                                                      annotLevel = 2)

gmkf_all_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                    sctSpecies = "human",
                                                    genelistSpecies = "human",
                                                    hits = gmfk_genes, 
                                                    reps = 10000,
                                                    geneSizeControl = TRUE,
                                                    annotLevel = 2)

gmkf_lof_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                    sctSpecies = "human",
                                                    genelistSpecies = "human",
                                                    hits = gmfk_lof_genes, 
                                                    reps = 10000,
                                                    geneSizeControl = TRUE,
                                                    annotLevel = 2)

gmkf_mis_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                    sctSpecies = "human",
                                                    genelistSpecies = "human",
                                                    hits = gmfk_mis_genes, 
                                                    reps = 10000,
                                                    geneSizeControl = TRUE,
                                                    annotLevel = 2)

gmkf_nosyn_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                      sctSpecies = "human",
                                                      genelistSpecies = "human",
                                                      hits = gmfk_nosyn_genes, 
                                                      reps = 10000,
                                                      geneSizeControl = TRUE,
                                                      annotLevel = 2)

gmkf_syn_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                    sctSpecies = "human",
                                                    genelistSpecies = "human",
                                                    hits = gmfk_syn_genes, 
                                                    reps = 10000,
                                                    geneSizeControl = TRUE,
                                                    annotLevel = 2)
ddd_all_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                   sctSpecies = "human",
                                                   genelistSpecies = "human",
                                                   hits = ddd_genes, 
                                                   reps = 10000,
                                                   geneSizeControl = TRUE,
                                                   annotLevel = 2)
ddd_lof_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                   sctSpecies = "human",
                                                   genelistSpecies = "human",
                                                   hits = ddd_lof_genes, 
                                                   reps = 10000,
                                                   geneSizeControl = TRUE,
                                                   annotLevel = 2)

ddd_mis_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                   sctSpecies = "human",
                                                   genelistSpecies = "human",
                                                   hits = ddd_mis_genes, 
                                                   reps = 10000,
                                                   geneSizeControl = TRUE,
                                                   annotLevel = 2)

ddd_nosyn_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                     sctSpecies = "human",
                                                     genelistSpecies = "human",
                                                     hits = ddd_nosyn_genes, 
                                                     reps = 10000,
                                                     geneSizeControl = TRUE,
                                                     annotLevel = 2)

ddd_syn_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                   sctSpecies = "human",
                                                   genelistSpecies = "human",
                                                   hits = ddd_syn_genes, 
                                                   reps = 10000,
                                                   geneSizeControl = TRUE,
                                                   annotLevel = 2)

ofc_panel_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                     sctSpecies = "human",
                                                     genelistSpecies = "human",
                                                     hits = ofc_panel_genes, 
                                                     reps = 10000,
                                                     geneSizeControl = TRUE,
                                                     annotLevel = 2)

human.bg <- ewceData::mouse_to_human_homologs()$HGNC.symbol
gene.list.2 <- sample(human.bg,size = 539)
second_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                  sctSpecies = "human",
                                                  sctSpecies_origin = "human",
                                                  genelistSpecies = "human",
                                                  hits = gene.list.2, 
                                                  reps = 10000,
                                                  annotLevel = 2,
                                                  geneSizeControl = TRUE)
novel_res2 = data.frame(novel_gene_results$results,
                        list="novel_genes")
cleftgenedb_res2 = data.frame(cleftgenedb_results$results,
                              list="cleftgenedb")
black_genes_res2 = data.frame(black_module_results$results,
                              list="black_module")
rando_res2 = data.frame(second_results$results,
                        list="Random")
gnomad_dec1_res2 = data.frame(gnomad_dec1_results$results,
                              list="Gnomad Decile 1")

gnomad_dec9_res2 = data.frame(gnomad_dec9_results$results,
                              list="Gnomad Decile 9")

gmkf_all_res2 = data.frame(gmkf_all_results$results,
                           list="GMKF ALL Genes")
gmkf_lof_res2 = data.frame(gmkf_lof_results$results,
                           list="GMKF LOF Genes")
gmkf_mis_res2 = data.frame(gmkf_mis_results$results,
                           list="GMKF MIS Genes")
gmkf_nosyn_res2 = data.frame(gmkf_nosyn_results$results,
                             list="GMKF NO SYN Genes")
gmkf_syn_res2 = data.frame(gmkf_syn_results$results,
                           list="GMKF SYN Genes")
ddd_all_res2 = data.frame(ddd_all_results$results,
                          list="ddd ALL Genes")
ddd_lof_res2 = data.frame(ddd_lof_results$results,
                          list="ddd LOF Genes")
ddd_mis_res2 = data.frame(ddd_mis_results$results,
                          list="ddd MIS Genes")
ddd_syn_res2 = data.frame(ddd_syn_results$results,
                          list="ddd SYN Genes")
ofc_panel_res2 = data.frame(ofc_panel_results$results,
                            list="OFC Panel Genes")


ofc_panel_res2 = data.frame(ofc_panel_results$results,
                            list="OFC Panel Genes")

library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)
#load Denisovan/Neanderthal segments
#obtain from Capra Lab https://github.com/emcarthur/trait-h2-neanderthals
library("rGREAT")

dev_bed_file <- "ALLEUR_sprime_segments_neanMatchingFilter.bed"
bed.gr <- import(dev_bed_file, genome = "hg19")

job2 = submitGreatJob(bed.gr,
                      genome = "hg19",
                      rule = c("oneClosest"))

tbl = getEnrichmentTables(job2)
gene_regions = getRegionGeneAssociations(job2, verbose = great_opt$verbose)
neanderthal_genes <- unique(as.data.frame(gene_regions$annotated_genes)$value)

zoohar_bed_file <- "zooHAR.hg38.bed"
har.gr <- import(zoohar_bed_file, genome = "hg38")

job3 = submitGreatJob(har.gr,
                      genome = "hg38",
                      rule = c("oneClosest"))
har.tbl = getEnrichmentTables(job3)
har_gene_regions = getRegionGeneAssociations(job3, verbose = great_opt$verbose)
all_har_genes <- unique(as.data.frame(har_gene_regions$annotated_genes)$value)

neanderthal_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                       sctSpecies = "human",
                                                       sctSpecies_origin = "human",
                                                       genelistSpecies = "human",
                                                       hits = neanderthal_genes, 
                                                       reps = 10000,
                                                       geneSizeControl = TRUE,
                                                       annotLevel = 2)

neanderthal_genes_res2 = data.frame(neanderthal_results$results,
                                    list="neanderthal_genes")

har_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                               sctSpecies = "human",
                                               sctSpecies_origin = "human",
                                               genelistSpecies = "human",
                                               hits = all_har_genes, 
                                               reps = 10000,
                                               geneSizeControl = TRUE,
                                               annotLevel = 2)

har_genes_res2 = data.frame(har_results$results,
                            list="har_genes")


merged_results = rbind(novel_res2, cleftgenedb_res2, black_genes_res2, neanderthal_genes_res2, gnomad_dec1_res2, gnomad_dec9_res2, gmkf_all_res2, gmkf_lof_res2, gmkf_mis_res2, gmkf_syn_res2, gmkf_nosyn_res2, ddd_all_res2, ddd_lof_res2, ddd_mis_res2, ddd_syn_res2, ofc_panel_res2, rando_res2)


merged_results$CellType[which(merged_results$CellType=="SOX2_SOX10+")] <-"SOX2-_SOX10+"
merged_results$CellType[which(merged_results$CellType=="SOX2+_SOX10_")] <-"SOX2+_SOX10-"


# saveRDS(merged_results, file = "/scr1/users/manchela/Data/EWCE_merged_results_all_celltype.rds")

merged_results <- readRDS("/scr1/users/manchela/Data/EWCE_merged_results_all_celltype.rds")

# plot_list <- EWCE::ewce_plot(total_res = merged_results,
#                              mtc_method = "BH",
#                              ctd = face_subtype_ctd,
#                              heights = c(0.1, 1.5),
#                              annotLevel = 1)



total_res = merged_results[-which(merged_results$list=="OFC Panel Genes"),]
q_threshold = 0.05
ast_q <- rep("", dim(total_res)[1])
ast_q[total_res$q < q_threshold] <- "*"
total_res$ast_q <- ast_q
total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
graph_theme <- ggplot2::theme_bw(base_size = 12, base_family = "Helvetica") + 
  ggplot2::theme(text = ggplot2::element_text(size = 14), 
                 axis.title.y = ggplot2::element_text(vjust = 0.6), 
                 strip.background = ggplot2::element_rect(fill = "white"), 
                 strip.text = ggplot2::element_text(color = "black"))
upperLim <- max(abs(total_res$sd_from_mean), na.rm = TRUE)
total_res$y_ast <- total_res$sd_from_mean * 1.05
total_res$abs_sd <- abs(total_res$sd_from_mean)

the_plot <- ggplot2::ggplot(total_res) + 
    ggplot2::geom_bar(ggplot2::aes_string(x = "CellType", y = "abs_sd", fill = "abs_sd"), stat = "identity") + 
    ggplot2::scale_fill_gradient(low = "blue", high = "red") + 
    graph_theme + ggplot2::theme(legend.position = "none") +
    theme(plot.margin = ggplot2::unit(c(0.5, 0, 0, 0), "mm"), 
          axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
          strip.text.y = ggplot2::element_text(angle = 0)) + 
    xlab("Cell type") + ylab("Std.Devs. from the mean") +
    scale_y_continuous(breaks = c(0, ceiling(upperLim * 0.66)), expand = expansion(mult = 0.4)) + 
    geom_text(ggplot2::aes_string(label = "ast_q",x = "CellType", y = "y_ast"), size = 10) +
    facet_grid("list ~ .", scales = "free_y", space = "free_x")

the_plot

pdf(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/ewce_disease_genes_main_types.pdf", h = 16, w = 11)
the_plot
dev.off()

total_res_heatmap_data <- as.data.frame(pivot_wider(total_res[,c(1,7,10)],names_from = CellType,values_from = abs_sd))
rownames(total_res_heatmap_data) <- total_res_heatmap_data[,1]
total_res_heatmap_data <- total_res_heatmap_data[,-1]

total_res_heatmap_sig <- as.data.frame(pivot_wider(total_res[,c(1,7,8)],names_from = CellType,values_from = ast_q))
rownames(total_res_heatmap_sig) <- total_res_heatmap_sig[,1]
total_res_heatmap_sig <- total_res_heatmap_sig[,-1]

p <- pheatmap(total_res_heatmap_data,cluster_rows = T,cluster_cols = T,silent = TRUE)


number_color <- ifelse(total_res_heatmap_data > 9, "white", "black")
number_color <- number_color[p$tree_row$order, p$tree_col$order]

annot_col <- data.frame(main_type=factor(colnames(total_res_heatmap_data),levels=c("mesenchyme","ectoderm","muscle","erythrocytes","endothelium",
                                                                                   "SOX2-_SOX10+","immune","SOX2+_SOX10-")))
rownames(annot_col) <- annot_col[,1]

library(scales)
annot_colors <-  hue_pal()(length(colnames(total_res_heatmap_data)))
annot_colors <- list(main_type = annot_colors)
names(annot_colors$main_type) <-  c("mesenchyme","ectoderm","muscle","erythrocytes","endothelium",
                                    "SOX2-_SOX10+","immune","SOX2+_SOX10-")

# Generate final heatmap with correct number colors
heat <- pheatmap(t(total_res_heatmap_data),
         display_numbers = t(total_res_heatmap_sig),
         angle_col = 45,
         color = colorRampPalette(c("white", "#950606"))(100),
         number_color = t(number_color),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_number = 12,
         annotation_row = annot_col,
         annotation_colors = annot_colors)

pdf(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/ewce_disease_genes_main_types_heatmap.pdf", h = 6, w = 8)
heat
dev.off()

#repeat analysis with subtype data
novel_gene_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                      sctSpecies = "human",
                                                      sctSpecies_origin = "human",
                                                      genelistSpecies = "human",
                                                      hits = novel_genes, 
                                                      reps = 10000,
                                                      geneSizeControl = TRUE,
                                                      annotLevel = 1)

cleftgenedb_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                       sctSpecies = "human",
                                                       sctSpecies_origin = "human",
                                                       genelistSpecies = "human",
                                                       hits = cleftgenedb, 
                                                       reps = 10000,
                                                       geneSizeControl = TRUE,
                                                       annotLevel = 1)

black_hubs_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                      sctSpecies = "human",
                                                      sctSpecies_origin = "human",
                                                      genelistSpecies = "human",
                                                      hits = black_hubs, 
                                                      reps = 10000,
                                                      geneSizeControl = TRUE,
                                                      annotLevel = 1)

black_module_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                        sctSpecies = "human",
                                                        sctSpecies_origin = "human",
                                                        genelistSpecies = "human",
                                                        hits = black_genes, 
                                                        reps = 10000,
                                                        geneSizeControl = TRUE,
                                                        annotLevel = 1)

gnomad_dec1_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                       sctSpecies = "human",
                                                       sctSpecies_origin = "human",
                                                       genelistSpecies = "human",
                                                       hits = gnomad_dec1_genes, 
                                                       reps = 10000,
                                                       geneSizeControl = TRUE,
                                                       annotLevel = 1)

gnomad_dec9_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                       sctSpecies = "human",
                                                       sctSpecies_origin = "human",
                                                       genelistSpecies = "human",
                                                       hits = gnomad_dec9_genes, 
                                                       reps = 10000,
                                                       geneSizeControl = TRUE,
                                                       annotLevel = 1)

gnomad_pli_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                      sctSpecies = "human",
                                                      sctSpecies_origin = "human",
                                                      genelistSpecies = "human",
                                                      hits = gnomad_pli_genes, 
                                                      reps = 10000,
                                                      geneSizeControl = TRUE,
                                                      annotLevel = 1)

gmkf_all_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                    sctSpecies = "human",
                                                    sctSpecies_origin = "human",
                                                    genelistSpecies = "human",
                                                    hits = gmfk_genes, 
                                                    reps = 10000,
                                                    geneSizeControl = TRUE,
                                                    annotLevel = 1)
gmkf_lof_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                    sctSpecies = "human",
                                                    sctSpecies_origin = "human",
                                                    genelistSpecies = "human",
                                                    hits = gmfk_lof_genes, 
                                                    reps = 10000,
                                                    geneSizeControl = TRUE,
                                                    annotLevel = 1)

gmkf_mis_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                    sctSpecies = "human",
                                                    sctSpecies_origin = "human",
                                                    genelistSpecies = "human",
                                                    hits = gmfk_mis_genes, 
                                                    reps = 10000,
                                                    geneSizeControl = TRUE,
                                                    annotLevel = 1)

gmkf_nosyn_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                      sctSpecies = "human",
                                                      sctSpecies_origin = "human",
                                                      genelistSpecies = "human",
                                                      hits = gmfk_nosyn_genes, 
                                                      reps = 10000,
                                                      geneSizeControl = TRUE,
                                                      annotLevel = 1)

gmkf_syn_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                    sctSpecies = "human",
                                                    sctSpecies_origin = "human",
                                                    genelistSpecies = "human",
                                                    hits = gmfk_syn_genes, 
                                                    reps = 10000,
                                                    geneSizeControl = TRUE,
                                                    annotLevel = 1)

ofc_panel_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                     sctSpecies = "human",
                                                     sctSpecies_origin = "human",
                                                     genelistSpecies = "human",
                                                     hits = ofc_panel_genes, 
                                                     reps = 10000,
                                                     geneSizeControl = TRUE,
                                                     annotLevel = 1)

ddd_all_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                   sctSpecies = "human",
                                                   sctSpecies_origin = "human",
                                                   genelistSpecies = "human",
                                                   hits = ddd_genes, 
                                                   reps = 10000,
                                                   geneSizeControl = TRUE,
                                                   annotLevel = 1)
ddd_lof_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                   sctSpecies = "human",
                                                   sctSpecies_origin = "human",
                                                   genelistSpecies = "human",
                                                   hits = ddd_lof_genes, 
                                                   reps = 10000,
                                                   geneSizeControl = TRUE,
                                                   annotLevel = 1)

ddd_mis_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                   sctSpecies = "human",
                                                   sctSpecies_origin = "human",
                                                   genelistSpecies = "human",
                                                   hits = ddd_mis_genes, 
                                                   reps = 10000,
                                                   geneSizeControl = TRUE,
                                                   annotLevel = 1)

ddd_syn_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                   sctSpecies = "human",
                                                   sctSpecies_origin = "human",
                                                   genelistSpecies = "human",
                                                   hits = gmfk_syn_genes, 
                                                   reps = 10000,
                                                   geneSizeControl = TRUE,
                                                   annotLevel = 1)

neanderthal_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                       sctSpecies = "human",
                                                       sctSpecies_origin = "human",
                                                       genelistSpecies = "human",
                                                       hits = neanderthal_genes, 
                                                       reps = 10000,
                                                       geneSizeControl = TRUE,
                                                       annotLevel = 1)

har_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                               sctSpecies = "human",
                                               sctSpecies_origin = "human",
                                               genelistSpecies = "human",
                                               hits = all_har_genes, 
                                               reps = 10000,
                                               geneSizeControl = TRUE,
                                               annotLevel = 1)



second_results <- EWCE::bootstrap_enrichment_test(sct_data = face_subtype_ctd,
                                                  sctSpecies = "human",
                                                  sctSpecies_origin = "human",
                                                  genelistSpecies = "human",
                                                  hits = gene.list.2, 
                                                  reps = 10000,
                                                  annotLevel = 1,
                                                  geneSizeControl = TRUE)
novel_res2 = data.frame(novel_gene_results$results,
                        list="03_novel_genes")
cleftgenedb_res2 = data.frame(cleftgenedb_results$results,
                              list="01_cleftgenedb")
black_genes_res2 = data.frame(black_module_results$results,
                              list="02_black_module")
rando_res2 = data.frame(second_results$results,
                        list="14_Random")
gnomad_dec1_res2 = data.frame(gnomad_dec1_results$results,
                              list="04_Gnomad Decile 1")

gnomad_dec9_res2 = data.frame(gnomad_dec9_results$results,
                              list="05_Gnomad Decile 9")

gmkf_all_res2 = data.frame(gmkf_all_results$results,
                           list="06_GMKF ALL Genes")
gmkf_lof_res2 = data.frame(gmkf_lof_results$results,
                           list="07_GMKF LOF Genes")
gmkf_mis_res2 = data.frame(gmkf_mis_results$results,
                           list="08_GMKF MIS Genes")
gmkf_nosyn_res2 = data.frame(gmkf_nosyn_results$results,
                             list="09_GMKF No SYN Genes")
gmkf_syn_res2 = data.frame(gmkf_syn_results$results,
                           list="10_GMKF SYN Genes")
ddd_all_res2 = data.frame(ddd_all_results$results,
                          list="ddd ALL Genes")
ddd_lof_res2 = data.frame(ddd_lof_results$results,
                          list="ddd LOF Genes")
ddd_mis_res2 = data.frame(ddd_mis_results$results,
                          list="ddd MIS Genes")
ddd_syn_res2 = data.frame(ddd_syn_results$results,
                          list="ddd SYN Genes")
ofc_panel_res2 = data.frame(ofc_panel_results$results,
                            list="11_OFC Panel Genes")
neanderthal_genes_res2 = data.frame(neanderthal_results$results,
                                    list="12_neanderthal_genes")
har_genes_res2 = data.frame(har_results$results,
                            list="13_har_genes")

merged_results = rbind(novel_res2, cleftgenedb_res2, black_genes_res2, neanderthal_genes_res2, har_genes_res2, gnomad_dec1_res2, gnomad_dec9_res2, gmkf_all_res2, gmkf_lof_res2, gmkf_mis_res2, gmkf_nosyn_res2, gmkf_syn_res2, ofc_panel_res2, rando_res2)

# saveRDS(merged_results, file = "/scr1/users/manchela/Data/EWCE_merged_results_all_subtype.rds")

merged_results <- readRDS("/scr1/users/manchela/Data/EWCE_merged_results_all_subtype.rds")

# plot_list <- EWCE::ewce_plot(total_res = merged_results,
#                              mtc_method = "BH",
#                              ctd = face_subtype_ctd,
#                              heights = c(0.1, 1.5),
#                              annotLevel = 1)

total_res = merged_results[-which(merged_results$list=="11_OFC Panel Genes"),]
q_threshold = 0.05
ast_q <- rep("", dim(total_res)[1])
ast_q[total_res$q < q_threshold] <- "*"
total_res$ast_q <- ast_q
total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
graph_theme <- ggplot2::theme_bw(base_size = 12, base_family = "Helvetica") + 
  ggplot2::theme(text = ggplot2::element_text(size = 14), 
                 axis.title.y = ggplot2::element_text(vjust = 0.6), 
                 strip.background = ggplot2::element_rect(fill = "white"), 
                 strip.text = ggplot2::element_text(color = "black"))
upperLim <- max(abs(total_res$sd_from_mean), na.rm = TRUE)
total_res$y_ast <- total_res$sd_from_mean * 1.05
total_res$abs_sd <- abs(total_res$sd_from_mean)

the_plot <- ggplot2::ggplot(total_res) + 
  ggplot2::geom_bar(ggplot2::aes_string(x = "CellType", y = "abs_sd", fill = "abs_sd"), stat = "identity") + 
  ggplot2::scale_fill_gradient(low = "blue", high = "red") + 
  graph_theme + ggplot2::theme(legend.position = "none") +
  theme(plot.margin = ggplot2::unit(c(0.5, 0, 0, 0), "mm"), 
        axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.text.y = ggplot2::element_text(angle = 0)) + 
  xlab("Cell type") + ylab("Std.Devs. from the mean") +
  scale_y_continuous(breaks = c(0, ceiling(upperLim * 0.66)), expand = expansion(mult = 0.6)) + 
  geom_text(ggplot2::aes_string(label = "ast_q",x = "CellType", y = "y_ast"), size = 10) +
  facet_grid("list ~ .", scales = "free_y", space = "free_x")

the_plot


pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/ewce_disease_genes_subtypes.pdf", h = 12, w = 18)
the_plot
dev.off()

total_res_heatmap_data <- as.data.frame(pivot_wider(total_res[,c(1,7,10)],names_from = CellType,values_from = abs_sd))
rownames(total_res_heatmap_data) <- total_res_heatmap_data[,1]
total_res_heatmap_data <- total_res_heatmap_data[,-1]

total_res_heatmap_sig <- as.data.frame(pivot_wider(total_res[,c(1,7,8)],names_from = CellType,values_from = ast_q))
rownames(total_res_heatmap_sig) <- total_res_heatmap_sig[,1]
total_res_heatmap_sig <- total_res_heatmap_sig[,-1]

p <- pheatmap(total_res_heatmap_data,cluster_rows = T,cluster_cols = T,silent = TRUE)


number_color <- ifelse(total_res_heatmap_data > 9, "white", "black")
number_color <- number_color[p$tree_row$order, p$tree_col$order]

annot_col <- data.frame(subtype=colnames(total_res_heatmap_data))
rownames(annot_col) <- annot_col[,1]

cds_subtype$subtype <- as.character(cds_subtype$subtype)
cds_subtype$subtype[which(cds_subtype$subtype=="other blood cells")] <- "immune"
cds_subtype$subtype[which(cds_subtype$subtype=="red blood cells")] <- "erythrocytes"
annot_col$main_type <- cds_subtype$celltype[match(rownames(annot_col),gsub("\\.","_",cds_subtype$subtype))]
annot_col$subtype <- NULL

library(scales)
annot_colors <-  hue_pal()(length(levels(annot_col$main_type)))
annot_colors <- list(main_type = annot_colors)
names(annot_colors$main_type) <-  c("mesenchyme","ectoderm","muscle","erythrocytes","endothelium",
                                    "SOX2-_SOX10+","immune","SOX2+_SOX10-")

# Generate final heatmap with correct number colors
library(limma)
rownames(total_res_heatmap_data) <- strsplit2(rownames(total_res_heatmap_data),split="[0-9][0-9]_")[,2]
heat <- pheatmap(t(total_res_heatmap_data),
                 display_numbers = t(total_res_heatmap_sig),
                 angle_col = 45,
                 color = colorRampPalette(c("white", "#950606"))(100),
                 number_color = t(number_color),  # Now matches clustered data
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 fontsize_number = 12,
                 annotation_row = annot_col,
                 annotation_colors = annot_colors,
                 fontsize_col = 10)

pdf(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/ewce_disease_genes_subtypes_heatmap.pdf", h = 16, w = 12)
heat
dev.off()


heat2 <- pheatmap(t(log2(total_res_heatmap_data+1)),
                 display_numbers = t(total_res_heatmap_sig),
                 angle_col = 45,
                 color = colorRampPalette(c("white", "#950606"))(100),
                 number_color = t(number_color),  # Now matches clustered data
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 fontsize_number = 12,
                 annotation_row = annot_col,
                 annotation_colors = annot_colors,
                 fontsize_col = 10)

pdf(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/ewce_disease_genes_subtypes_heatmap_log2.pdf", h = 16, w = 12)
heat2
dev.off()

#shortend list
merged_results = rbind(novel_res2, cleftgenedb_res2, black_genes_res2, neanderthal_genes_res2, har_genes_res2, gnomad_dec1_res2, gnomad_dec9_res2, gmkf_nosyn_res2, gmkf_syn_res2, ofc_panel_res2, rando_res2)

# plot_list <- EWCE::ewce_plot(total_res = merged_results,
#                              mtc_method = "BH",
#                              q_threshold = 0.05,
#                              ctd = face_subtype_ctd,
#                              heights = c(0.1, 2),
#                              annotLevel = 1)

# print(plot_list$plain)
total_res = merged_results[-which(merged_results$list%in%c("06_GMKF ALL Genes","07_GMKF LOF Genes","08_GMKF MIS Genes","11_OFC Panel Genes")),]
q_threshold = 0.05
ast_q <- rep("", dim(total_res)[1])
ast_q[total_res$q < q_threshold] <- "*"
total_res$ast_q <- ast_q
total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
graph_theme <- ggplot2::theme_bw(base_size = 12, base_family = "Helvetica") + 
  ggplot2::theme(text = ggplot2::element_text(size = 14), 
                 axis.title.y = ggplot2::element_text(vjust = 0.6), 
                 strip.background = ggplot2::element_rect(fill = "white"), 
                 strip.text = ggplot2::element_text(color = "black"))
upperLim <- max(abs(total_res$sd_from_mean), na.rm = TRUE)
total_res$y_ast <- total_res$sd_from_mean * 1.05
total_res$abs_sd <- abs(total_res$sd_from_mean)

the_plot <- ggplot2::ggplot(total_res) + 
  ggplot2::geom_bar(ggplot2::aes_string(x = "CellType", y = "abs_sd", fill = "abs_sd"), stat = "identity") + 
  ggplot2::scale_fill_gradient(low = "blue", high = "red") + 
  graph_theme + ggplot2::theme(legend.position = "none") +
  theme(plot.margin = ggplot2::unit(c(0.5, 0, 0, 0), "mm"), 
        axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.text.y = ggplot2::element_text(angle = 0)) + 
  xlab("Cell type") + ylab("Std.Devs. from the mean") +
  scale_y_continuous(breaks = c(0, ceiling(upperLim * 0.66)), expand = expansion(mult = 0.4)) + 
  geom_text(ggplot2::aes_string(label = "ast_q",x = "CellType", y = "y_ast"), size = 10) +
  facet_grid("list ~ .", scales = "free_y", space = "free_x")

the_plot

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/ewce_disease_genes_subtypes_shortened.pdf", h = 10, w = 16)
the_plot
dev.off()


total_res_heatmap_data <- as.data.frame(pivot_wider(total_res[,c(1,7,10)],names_from = CellType,values_from = abs_sd))
rownames(total_res_heatmap_data) <- total_res_heatmap_data[,1]
total_res_heatmap_data <- total_res_heatmap_data[,-1]

total_res_heatmap_sig <- as.data.frame(pivot_wider(total_res[,c(1,7,8)],names_from = CellType,values_from = ast_q))
rownames(total_res_heatmap_sig) <- total_res_heatmap_sig[,1]
total_res_heatmap_sig <- total_res_heatmap_sig[,-1]

p <- pheatmap(total_res_heatmap_data,cluster_rows = T,cluster_cols = T,silent = TRUE)


number_color <- ifelse(total_res_heatmap_data > 9, "white", "black")
number_color <- number_color[p$tree_row$order, p$tree_col$order]

annot_col <- data.frame(subtype=colnames(total_res_heatmap_data))
rownames(annot_col) <- annot_col[,1]

cds_subtype$subtype <- as.character(cds_subtype$subtype)
cds_subtype$subtype[which(cds_subtype$subtype=="other blood cells")] <- "immune"
cds_subtype$subtype[which(cds_subtype$subtype=="red blood cells")] <- "erythrocytes"
annot_col$main_type <- cds_subtype$celltype[match(rownames(annot_col),gsub("\\.","_",cds_subtype$subtype))]
annot_col$subtype <- NULL

library(scales)
annot_colors <-  hue_pal()(length(levels(annot_col$main_type)))
annot_colors <- list(main_type = annot_colors)
names(annot_colors$main_type) <-  c("mesenchyme","ectoderm","muscle","erythrocytes","endothelium",
                                    "SOX2-_SOX10+","immune","SOX2+_SOX10-")

# Generate final heatmap with correct number colors
library(limma)
rownames(total_res_heatmap_data) <- strsplit2(rownames(total_res_heatmap_data),split="[0-9][0-9]_")[,2]
heat <- pheatmap(t(total_res_heatmap_data),
                 display_numbers = t(total_res_heatmap_sig),
                 angle_col = 45,
                 color = colorRampPalette(c("white", "#950606"))(100),
                 number_color = t(number_color),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 fontsize_number = 12,
                 annotation_row = annot_col,
                 annotation_colors = annot_colors,
                 fontsize_col = 10)


pdf(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/ewce_disease_genes_subtypes_shortened_heatmap.pdf", h = 16, w = 10)
heat
dev.off()

heat2 <- pheatmap(t(log2(total_res_heatmap_data+1)),
                 display_numbers = t(total_res_heatmap_sig),
                 angle_col = 45,
                 color = colorRampPalette(c("white", "#950606"))(100),
                 number_color = t(number_color),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 fontsize_number = 12,
                 annotation_row = annot_col,
                 annotation_colors = annot_colors,
                 fontsize_col = 10)


pdf(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/ewce_disease_genes_subtypes_shortened_heatmap_log2.pdf", h = 16, w = 10)
heat2
dev.off()
