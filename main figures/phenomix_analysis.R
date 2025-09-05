library(EWCE)
library(HPOExplorer)
library(MSTExplorer)
library(ggplot2)
library(readr)
library(data.table)
library(dplyr)
library(tidyr)
library(scales)

# for figure 7C

results <- read_csv("/mnt/isilon/cotney_lab_hpc/manchela/misc/phenomix_results_sig.csv")
results <- setDT(results) 

hpo <- HPOExplorer::get_hpo()
ymat <- HPOExplorer::hpo_to_matrix()

test_method <- "glm"
metric <- "specificity"

source("/home/manchela/Downloads/plot_bar_dendrogram.R")
source("/home/manchela/Downloads/run_prop_test.R")
source("/home/manchela/Downloads/plot_bar.R")

ctd <- EWCE::load_rdata("/home/manchela/Downloads/ctd_face_ctd.rda")

ctd_name <- "ctd_face"

results_face <- MSTExplorer::run_phenomix(ymat = ymat,
                                          ctd = ctd,
                                          ctd_name = ctd_name,
                                          test_method = test_method,
                                          metric = metric,
                                          save_path = paste0("/scr1/users/manchela/Data/Multiome/Humanface/GEX/phenomix/phenomix_glm_specificity",
                                          "phenomix_",ctd_name,"_results.tsv.gz")
                                          )
saveRDS(results_face,"/scr1/users/manchela/Data/Multiome/Humanface/GEX/phenomix/phenomix_glm_specificity/results_face.rds")

ctd <- EWCE::load_rdata("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/ctd_face_subtype_sctransform_ctd.rda")
ctd_name <- "ctd_face_subtype"

results_face_subtype <- MSTExplorer::run_phenomix(ymat = ymat, 
                                                  ctd = ctd, 
                                                  ctd_name = ctd_name, 
                                                  test_method = test_method, 
                                                  metric = metric,
                                                  save_path = paste0("/scr1/users/manchela/Data/Multiome/Humanface/GEX/phenomix/phenomix_glm_specificity",
                                                                     "phenomix_",ctd_name,"_results.tsv.gz")
                                                  )
# saveRDS(results_face_subtype,"/scr1/users/manchela/Data/Multiome/Humanface/GEX/phenomix/phenomix_glm_specificity/results_face_subtype.rds")
results_face_subtype <- readRDS("/scr1/users/manchela/Data/Multiome/Humanface/GEX/phenomix/phenomix_glm_specificity/results_face_subtype.rds")

results_face_subtype <- results_face_subtype[-which(results_face_subtype$annotLevel==3),]

use_both_ctd <- FALSE
if(use_both_ctd){
  results <- data.table::rbindlist(list(face=results_face, 
                                        face_subtype=results_face_subtype), 
                                   idcol = "ctd")
} else {
  results <-  data.table::rbindlist(list(face_subtype=results_face_subtype), 
                                    idcol = "ctd")
}

results <- HPOExplorer::add_hpo_name(results)

results <- HPOExplorer::add_ancestor(results, 
                                     keep_descendants = "Phenotypic abnormality")

# data.table::fwrite(results,"/scr1/users/manchela/Data/Multiome/Humanface/GEX/phenomix/phenomix_glm_specificity/phenomix_results_all.csv.gz")


cell_ordering <- c(ctd[[1]]$plotting$cell_ordering, 
                   ctd[[2]]$plotting$cell_ordering) 
overlapping_celltypes <- intersect(ctd[[1]]$plotting$cell_ordering, ctd[[2]]$plotting$cell_ordering)

results[,N_celltypes:=data.table::uniqueN(CellType),
        by= c("ctd","annotLevel")]  
results <- results[hpo_name!="Phenotypic abnormality"]

results$annotLevel <- as.numeric(results$annotLevel)

gg_bar <- plot_bar_dendrogram(results,
                              hpo = hpo,
                              ctd = ctd,
                              barplot_margins = list(c(0,0,0,0), c(0,0,0.4,0))
)

print(gg_bar$plot)

pdf("/scr1/users/manchela/Data/Multiome/Humanface/GEX/phenomix/phenomix_subtypes_barplot.pdf",height=20,width=24)
print(gg_bar$plot)
dev.off()


mc.cores = 1
formula=ancestor_label~ctd_label
hpo=HPOExplorer::get_hpo()
prune=TRUE
run_tests=TRUE
test="prop_test"
y_var=c("N",
        "N_normalised",
        "Proportion of phenotypes")[2]
y_lab=if(y_var=="N"){
  "Significant associations"
} else if (y_var=="N_normalised"){
  "Significant associations (normalised)"
} else if (y_var=="Proportion of phenotypes"){
  "Proportion of significant associations"
}
results[,N_phenotypes:=data.table::uniqueN(hpo_id),
        by=c("ancestor_name")] 
results <- KGExplorer::prune_ancestors(results, 
                                       id_col = "ancestor",
                                       ont=hpo)
test_res <- if(run_tests){
  suppressWarnings(
    run_prop_test(results = results, 
                  test = test,
                  mc.cores = mc.cores)
  )
}
results_agg <- results[q<0.05,
                       list(N=.N,
                            N_phenotypes,
                            `Proportion of phenotypes`=.N/N_phenotypes,
                            ancestor_label=paste0(ancestor_name,
                                                  "\n(",N_phenotypes," phenotypes)"),
                            ctd_label=paste0(ctd,"\n(level ",annotLevel,")\n",
                                             N_celltypes," cell types")),
                       by=c("ancestor_name","ctd","annotLevel","CellType")] |>
  unique()
results_agg[,N_normalised:=N/max(N, na.rm = TRUE), by=c("ancestor_name"), drop=FALSE]

# Sort rows by N
results_agg <- results_agg[order(N, decreasing = TRUE)] 

if(!is.null(test_res)){
  results_agg <- merge(results_agg,
                       test_res[q<0.05,, drop=FALSE], all.x = TRUE)
} 
results_agg[,ancestor_label:=paste0(ancestor_name, "\n(",N_phenotypes," phenotypes)"),]

# Make ancestor labels factors
results_agg$ancestor_label <- factor(results_agg$ancestor_label, 
                                     levels=unique(results_agg$ancestor_label),
                                     ordered = TRUE)

max_xtick_len <- sapply(ctd, FUN=function(x) max(nchar(as.character(x$plotting$cell_ordering)))) |> max()
# Relabel all xaxis ticks so that they have exactly the same number of characters
results_agg$CellType <- sapply(results_agg$CellType, 
                               FUN=function(x) stringr::str_pad(as.character(x), 
                                                                width = max_xtick_len, 
                                                                side = "left")) |> as.character()

cell_ordering_list <- lapply(ctd, FUN=function(x) stringr::str_pad(as.character(x$plotting$cell_ordering), 
                                                                   width = max_xtick_len, 
                                                                   side = "left"))


lvl = 1
results_agg_lvl <- results_agg[annotLevel==lvl,, drop=FALSE] |> data.table::copy()
# Order the cell types 
results_agg_lvl$CellType <- factor(results_agg_lvl$CellType, 
                                   levels=cell_ordering_list[[lvl]],
                                   ordered = TRUE)

ancestor_name = unique(results_agg_lvl$ancestor_name)
celltype = unique(results_agg_lvl$CellType)

sig_mat = results_agg_lvl[,c("ancestor_name","CellType","q_sig")]

sig_mat_wide <- sig_mat %>%
  pivot_wider(names_from = CellType, values_from = q_sig, values_fill = NA)

sig_mat_wide <- as.matrix(sig_mat_wide)
rownames(sig_mat_wide) <- sig_mat_wide[,1]
sig_mat_wide <- sig_mat_wide[,-1]

sig_mat_wide[which(is.na(sig_mat_wide))] <- ""

sobj <- readRDS("/scr1/users/manchela/Data/human_face_no-neuro_clustering_cellrangerARC-raw_emptyDrops_singlets_finalannot_27Mar.rds")
results_agg_lvl$main_type <- sobj$celltype[match(gsub(" ","",as.character(results_agg_lvl$CellType)),sobj$subtype)]

annot_colors <-  hue_pal()(8)
annot_colors <- list(main_type = annot_colors)
names(annot_colors$main_type) <- levels(results_agg_lvl$main_type)

library(tidyheatmaps)
p1 <- tidyheatmap(df = results_agg_lvl,
                  rows = ancestor_name,
                  columns = CellType,
                  values = N_normalised,
                  scale = "none",
                  annotation_col = main_type,
                  cluster_rows = FALSE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  fontsize_number = 14,
                  fontsize_col = 12,
                  fontsize_row = 12,
                  border_color = "grey",
                  color_na = "white",
                  display_numbers = sig_mat_wide,
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  colors = c("white","#950606"),
                  annotation_colors = annot_colors
                  )

pdf("/scr1/users/manchela/Data/Multiome/Humanface/GEX/phenomix/phenomix_subtypes_heatmap.pdf",height=8,width=24)
p1
dev.off()

