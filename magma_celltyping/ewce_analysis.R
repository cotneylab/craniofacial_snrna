library(EWCE)
library(SingleCellExperiment)
library(scRNAseq)
library(scater)
library(Seurat)

cds_subtype <- readRDS("../cds_face_human_temp.rds")
#normalize data for use in EWCE
cds_subtype <- PercentageFeatureSet(cds_subtype, pattern = "^MT-", col.name = "percent.mt")
cds_subtype <- SCTransform(cds_subtype, vars.to.regress = "percent.mt", verbose = FALSE)

sce_subtype <- as.SingleCellExperiment(cds_subtype, assay = 'SCT')

annotLevels = list(level1class=sce_subtype$tempanno,
                   level2class=sce_subtype$cell_type1)

Face_Subtype_SCtransform_CTD <- generate_celltype_data(exp=assays(sce_subtype)$counts,
                                           annotLevels=annotLevels,
                                           specificity_quantiles = TRUE,
                                           groupName="face_subtype_sctransform_ctd",
                                           no_cores = 4,
                                           input_species = "human",
                                           output_species = "human",
                                           return_ctd = TRUE,
                                           savePath="~/Desktop/scRNA-Seq_GWAS")

face_subtype_ctd <- EWCE::load_rdata("../ctd_face_subtype_sctransform_ctd.rda") 

filenames <- list.files(path = "../modules",pattern="*_Nodes.csv", full.names = TRUE)

modules <- lapply(setNames(filenames, make.names(gsub("*_Nodes.csv", "", filenames))),read.csv)
modules_symbols <- lapply(modules, function(x){x[,2]})
hubs  <- lapply(modules, function(x) { x[ x$HUB == 1, 2] })
black_hubs <- hubs$...modules.black
black_genes <- modules$...modules.black$SYMID

all_nsclp_genes <- readxl::read_xlsx("../modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 1)
nsclp_genes <- read.table("../modules/all_nsclp_genes.txt", header=TRUE, sep ="\t")
clp_genes <- read.table("../modules/clp_genes.txt", header=TRUE, sep ="\t")
all_clp_genes <- readxl::read_xlsx("../modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 2)
cp_genes <- read.table("../modules/cp_genes.txt", header=TRUE, sep ="\t")
all_cp_genes <- readxl::read_xlsx("../modules/20230110_CLP_DNM_list_synonymous_removed.xlsx", sheet = 3)
cleftgenedb <- read.table("../modules/cleftgenedb.txtr", header = TRUE)
cleftgenedb$Gene = toupper(cleftgenedb$Gene)
cleftgenedb <- cleftgenedb$Gene
novel_genes <- readxl::read_xlsx("../modules/Supplementary_Table_S7_update.xlsx", sheet = 2)
novel_genes <- novel_genes$SYMID
cfse_genes <- read.table("../modules/all_cfse_genes.txt", header = FALSE)
gnomad <- read.table("../gnomad.v4.1.constraint_metrics.tsv", header = TRUE)
gnomad_dec1 <- subset(gnomad, lof.oe_ci.upper_bin_decile == 1)
gnomad_dec9 <- subset(gnomad, lof.oe_ci.upper_bin_decile == 9)
gnomad_pli <- subset(gnomad, lof.pLI == 1)
gnomad_dec1_genes <- unique(gnomad_dec1$gene)
gnomad_dec9_genes <- unique(gnomad_dec9$gene)
gnomad_pli_genes <- unique(gnomad_pli$gene)

gmfk_all <- read.table("20241024_CPSeqGMFKDDD_allDNs.txt", header=TRUE, sep ="\t")
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

ofc_panels <- read.table("20241028_OFC_GenePanelList.txt", header=TRUE, sep ="\t")
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

dev_bed_file <- "../ALLEUR_sprime_segments_neanMatchingFilter.bed"
bed.gr <- import(dev_bed_file, genome = "hg19")

job2 = submitGreatJob(bed.gr,
                      genome = "hg19",
                      rule = c("oneClosest"))

tbl = getEnrichmentTables(job2)
gene_regions = getRegionGeneAssociations(job2, verbose = great_opt$verbose)
neanderthal_genes <- unique(as.data.frame(gene_regions$annotated_genes)$value)

zoohar_bed_file <- "../zooHAR.hg38.bed"
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

plot_list <- EWCE::ewce_plot(total_res = merged_results,
                             mtc_method = "BH",
                             ctd = face_subtype_ctd,
                             heights = c(0.1, 1.5),
                             annotLevel = s)

print(plot_list$plain)
pdf(file = "ewce_disease_genes_main_types.pdf", h = 8.5, w = 11)
print(plot_list$plain)
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
plot_list <- EWCE::ewce_plot(total_res = merged_results,
                             mtc_method = "BH",
                             ctd = face_subtype_ctd,
                             heights = c(0.1, 1.5),
                             annotLevel = 1)

print(plot_list$plain)
pdf(file = "ewce_disease_genes_subtypes_suppfigure.pdf", h = 11, w = 8.5)
print(plot_list$plain)
dev.off()

#shortend list
merged_results = rbind(novel_res2, cleftgenedb_res2, black_genes_res2, neanderthal_genes_res2, har_genes_res2, gnomad_dec1_res2, gnomad_dec9_res2, gmkf_nosyn_res2, gmkf_syn_res2, ofc_panel_res2, rando_res2)

plot_list <- EWCE::ewce_plot(total_res = merged_results,
                             mtc_method = "BH",
                             q_threshold = 0.05,
                             ctd = face_subtype_ctd,
                             heights = c(0.1, 2),
                             annotLevel = 1)

print(plot_list$plain)
pdf(file = "ewce_disease_genes_subtypes_mainfigure.pdf", h = 11, w = 8.5)
print(plot_list$plain)
dev.off()

saveRDS(merged_results, file = "merged_results_short.rds")
