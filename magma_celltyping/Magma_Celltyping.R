#neurogenomicslab/magma.celltyping
Sys.setenv(GITHUB_PAT="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")
install.packages("piggyback")
Sys.setenv('R_MAX_VSIZE'=64000000000)
# if(!require("remotes")) install.packages("remotes")


# for figure 7

# face_ctd <- EWCE::load_rdata("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/ctd_face_ctd.rda") 

face_subtype_ctd <- EWCE::load_rdata("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/ctd_face_subtype_sctransform_ctd.rda") 

storage_dir <- c("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA")

setwd("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/")
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_ANKYLOGLOSSIA.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_APVR.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_ASD.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 449747,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_AND_CARIES.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_HARD_PALATE_NOS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453074,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_CLEFT_PALATE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_OR_LIP_AND_PALATE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453074,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_PALATE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453539,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_DEFORMITI_FEET.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451335,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_LENS_MALFO.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451518,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_ANTER_SEGMENT_EYE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451621,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_CIRCULATO_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_DEFORMAT_MUSCULOS_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_EAR_CAUSI_IMPAIRM_HEARING.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451480,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GALLB_BILE_DUCTS_LIVER.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452554,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GENITAL_ORGANS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GREAT_ARTERIES.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 449138,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GREAT_VEINS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 448543,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_NERVOUS_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_OVARIES_FALLOP_TUBES_BROAD_LIGAM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 253663,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_POSTERI_SEGMENT_EYE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451518,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_PULMONARY_TRICU_VALVES.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 448693,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_RETI.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451442,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_TRACHEA_BRONCHUS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453484,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_URINARY_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONOTR_DEFEC.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CYSTIC_KIDNEY_DISEA.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452843,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_LVOTO_BROAD.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_LVOTO_NARROW.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OSTEOCHONDROD_W_DEFECTS_GROWTH_TUBULAR_BONES_SPINE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450566,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CHROMOSOME_ABNORMALITI_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453278,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_BRAIN.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453157,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_CIRCULATO_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 449240,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_DIGES_SYSTEM1.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_DIGES_SYSTEM2.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452357,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_EAR.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451803,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_FACE_NECK.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452348,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_FEMALE_GENITALIA.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450456,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_HEART.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 448657,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_INTESTINE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452555,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_KIDNEY.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451936,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_LIMBS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450750,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_MALE_GENITAL_ORGANS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450370,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_NERVOUS_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453326,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450857,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_PERIP_VASCULAR_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 448862,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_RESPI_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453246,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_SKIN.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451959,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_SKULL_FACE_BONES.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450601,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_TONGUE_MOUTH_PHARYNX.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452740,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_UPPER_ALIME_TRACT.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452503,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_URINARY_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451965,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_OSTEOCHONDROD.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450779,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_SPECIFE_CONGEN_MALFO_SYNDR_AFFECTING_MULTIPLE_SYSTEMS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451126,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_PHAKOMATOSES_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450945,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_RVOTO.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_SEPTA_DEFEC.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_UNDES_TESTICLE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450490,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_RX_CROHN_1STLINE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_RX_CROHN_2NDLINE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_SLE_FG.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 339484,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_CHRONLARGE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 434250,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_CHRONNAS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 434744,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_CHRONOTH.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 434933,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
magma_dirs_finngen_immune = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_CHRONLARGE.formatted.tsv.bgz",
                              "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_SLE_FG.formatted.tsv.bgz")
magma_dirs_finngen_other = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_ANKYLOGLOSSIA.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_APVR.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_ASD.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_DEFORMITI_FEET.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_LENS_MALFO.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_ANTER_SEGMENT_EYE.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_CIRCULATO_SYSTEM.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_DEFORMAT_MUSCULOS_SYSTEM.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_EAR_CAUSI_IMPAIRM_HEARING.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GALLB_BILE_DUCTS_LIVER.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GENITAL_ORGANS.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GREAT_ARTERIES.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GREAT_VEINS.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_NERVOUS_SYSTEM.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_OVARIES_FALLOP_TUBES_BROAD_LIGAM.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_POSTERI_SEGMENT_EYE.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_PULMONARY_TRICU_VALVES.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_RETI.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_TRACHEA_BRONCHUS.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_URINARY_SYSTEM.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CONOTR_DEFEC.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CYSTIC_KIDNEY_DISEA.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_LVOTO_BROAD.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_LVOTO_NARROW.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OSTEOCHONDROD_W_DEFECTS_GROWTH_TUBULAR_BONES_SPINE.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CHROMOSOME_ABNORMALITI_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_DIGES_SYSTEM1.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_EAR.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_FACE_NECK.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_HEART.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_SKIN.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_TONGUE_MOUTH_PHARYNX.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_UPPER_ALIME_TRACT.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_OSTEOCHONDROD.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_SPECIFE_CONGEN_MALFO_SYNDR_AFFECTING_MULTIPLE_SYSTEMS.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_PHAKOMATOSES_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_RVOTO.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_SEPTA_DEFEC.formatted.tsv.bgz",
                             "/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_UNDES_TESTICLE.formatted.tsv.bgz")


#spot check manhattan plots
library(data.table)
library(dplyr)
library(manhattanly)

finngen_cleft = fread("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_CLEFT_PALATE.formatted.tsv.bgz")
finngen_cleft.subset <- finngen_cleft %>% 
  filter(-log10(P)>1)

head(finngen_cleft.subset)

qqman::manhattan(finngen_cleft.subset, chr="CHR", bp="BP", snp="SNP", p="P", annotatePval = 0.0001)

#interactive version
manhattanly(finngen_cleft.subset, chr="CHR", bp="BP", snp="SNP", p="P", gene = "NEAREST_GENES")

library(locuszoomr)
library(LDlinkR)
library(EnsDb.Hsapiens.v86)
#API token aba03b266910
token <- c("aba03b266910")


loc <- locus(data = finngen_cleft, chrom = "CHR", labs = "SNP", pos = "BP", p = "P", gene = 'IRF6', index_snp = "rs642961", flank = 1e5, ens_db = "EnsDb.Hsapiens.v86")
loc <- link_LD(loc, pop = "FIN", genome_build = "grch38", r2d = "r2", method = "proxy", token = token)
locus_plot(loc, label_x = c(4, -5), labels = c("index", "rs570516915"))

magma_dirs_finngen_cleft  = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_AND_CARIES.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_HARD_PALATE_NOS.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_CLEFT_PALATE.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_OR_LIP_AND_PALATE.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_PALATE.formatted.tsv.bgz")
magma_dirs_finngen = c(magma_dirs_finngen_cleft, magma_dirs_finngen_other)
magma_dirs_finngen_alt = c(magma_dirs_finngen_cleft, magma_dirs_finngen_other,magma_dirs_finngen_immune)
magma_dirs_xiong = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.AlL-AlR.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.AlL-Ls.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.AlL-Sn.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.AlR-ChL.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.EnL-AlL.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.EnL-Sn.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.EnR-AlR.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.EnR-ChL.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.EnR-ChR.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.EnR-N.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.EnR-Prn.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.EnR-Sn.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.ExL-AlL.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.ExR-ChL.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.ExR-ChR.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.ExR-N.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.N-En.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.N-ExL.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.N-Prn.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.Prn-AlL.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.Prn-AlR.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.Prn-EnL.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.Prn-Ls.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/XiongZ_31763980.Sn-AlR.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/GCST011096.h.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/28067908-GCST004132-EFO_0000384.crohns.h.formatted.tsv.bgz")
magma_dirs_bonfontae = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Brow_Ridge_Protrusion_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Brow_Ridge_Protrusion_2.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Brow_Ridge_Protrusion_3.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Chin_Protrusion_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Chin_Protrusion_2.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Columella_Inclination.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Columella_Size.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Ear_Size.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Eye_Position_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Forehead_Protrusion_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Jaw_Protrusion_2.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Jaw_Protrusion_5.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Jaw_Slope_2.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Lip_Protrusion.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Lip_Thickness_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Lip_Thickness_Ratio_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Lip_Thickness_Ratio_2.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Lower_Face_Flatness.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Lower_Lip_Protrusion.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Lower_Lip_Thickness_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Lower_Lip_Thickness_2.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nasion_Depth_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nasion_Depth_2.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nasion_Position_2.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nose_Height.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nose_Protrusion.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nose_Roundness_1.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nose_Roundness_3.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nose_Size.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Nostril_Size.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Philtrum_Length.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/Upper_Face_Flatness.formatted.tsv.bgz")
magma_dirs_all_varition = c(magma_dirs_xiong, magma_dirs_bonfontae,magma_dirs_finngen_immune)
magma_dirs_nsclp = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/GMKF_CPseq_allpopulations_TDT_plinkSE_050124.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.meta.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/GCST011096.h.formatted.tsv.bgz")
#magma_dirs_nsclp = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/GMKF_CPseq_allpopulations_TDT_plinkSE_050124.com.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/GMKF_CPseq_allpopulations_TDT_plinkSE_050124.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/GMKF_CPseq_allpopulations_TDT_plinkSE_050124.par.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.beaty.tdt.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.marazita.cc.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.marazita.tdt.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.meta.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.meta.marazita.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/GCST011096.h.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/28067908-GCST004132-EFO_0000384.crohns.h.formatted.tsv.bgz")
magma_dirs_microsomia = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/craniofacial_microsomia.formatted.tsv.bgz")
magma_dirs_immune = c("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/GCST011096.h.formatted.tsv.bgz","/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/28067908-GCST004132-EFO_0000384.crohns.h.formatted.tsv.bgz")
magma_dirs_craniofacial_disease =c(magma_dirs_nsclp,magma_dirs_microsomia,magma_dirs_immune)
magma_dirs_craniofacial_disease_variation = c(magma_dirs_nsclp,magma_dirs_microsomia,magma_dirs_xiong)
magma_dirs_craniofacial_disease_variation2 = c(magma_dirs_nsclp,magma_dirs_microsomia,magma_dirs_immune,magma_dirs_bonfontae)
# test = c(magma_dirs_craniofacial_disease, magma_dirs_finngen_subset)


MAGMA_results_finngen_subtypes <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs =   magma_dirs_finngen_alt,
  ctd = face_subtype_ctd,
  ctd_species = "human", 
  ctd_name = "face_subtype_ctd",
  upstream_kb = 100,
  downstream_kb = 100,
  nThread = 8,
  force_new = TRUE,
  run_linear = TRUE,
  run_conditional = FALSE)

names(MAGMA_results_finngen_subtypes) <- gsub("finngen_R11_","", names(MAGMA_results_finngen_subtypes))
names(MAGMA_results_finngen_subtypes) <- gsub("Q17_","", names(MAGMA_results_finngen_subtypes))
names(MAGMA_results_finngen_subtypes) <- gsub(".formatted.tsv.bgz","", names(MAGMA_results_finngen_subtypes))
# saveRDS(MAGMA_results_finngen_subtypes, "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/MAGMA_results_finngen_subtypes.rds")

MAGMA_results_finngen_subtypes <- readRDS("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/MAGMA_results_finngen_subtypes.rds")




library(ggplot2)
library(tidyverse)
library(tidyheatmaps)
library(scales)

#subtype finngen analysis
merged_res <- MAGMA.Celltyping::merge_results(MAGMA_results_finngen_subtypes, level = 1, filetype = "ctAssocsLinear")

heat <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res,
                                          fdr_thresh = 1, fill_var = "-log10p")
merged_res_2 <- MAGMA.Celltyping::merge_results(MAGMA_results_finngen_subtypes, level = 2, filetype = "ctAssocsLinear")

heat_2 <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res_2,
                                            fdr_thresh = 1, fill_var = "-log10p")

merged_res_3 <- MAGMA.Celltyping::merge_results(MAGMA_results_finngen_subtypes, level = 3, filetype = "ctAssocsLinear")

heat_3 <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res_3,
                                            fdr_thresh = 1, fill_var = "-log10p")

# setwd("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS")
# cell_type_metadata <- readxl::read_xlsx("annotation_expanded.xlsx", sheet = 1)
initial.filt_singlet_noneuro <- readRDS("/scr1/users/manchela/Data/human_face_no-neuro_clustering_cellrangerARC-raw_emptyDrops_singlets_finalannot_27Mar.rds")
initial.filt_singlet_noneuro$subtype_noprog <- as.character(initial.filt_singlet_noneuro$subtype)
initial.filt_singlet_noneuro$subtype_noprog[which(initial.filt_singlet_noneuro$celltype=="SOX2+_SOX10-")] <- "SOX2+_SOX10-"
# initial.filt_singlet_noneuro$subtype_noprog <- factor(initial.filt_singlet_noneuro$subtype_noprog,levels=c(levels(initial.filt_singlet_noneuro$subtype)[1:70],"SOX2+_SOX10-"))

cell_type_metadata <- data.frame(subtype= initial.filt_singlet_noneuro@meta.data[,"subtype"],
                                 main_type=initial.filt_singlet_noneuro@meta.data[,"celltype"],
                                 subtype_noprog= initial.filt_singlet_noneuro@meta.data[,"subtype_noprog"])
cell_type_metadata <- cell_type_metadata[-which(duplicated(cell_type_metadata)),]
cell_type_metadata$subtype <- gsub("\\.","_",cell_type_metadata$subtype)
cell_type_metadata$subtype <- gsub(" ","_",cell_type_metadata$subtype)
cell_type_metadata$subtype_noprog <- gsub("\\.","_",cell_type_metadata$subtype_noprog)
cell_type_metadata$subtype_noprog <- gsub(" ","_",cell_type_metadata$subtype_noprog)



# cell_type_metadata$subtype[which(cell_type_metadata$subtype=="other_blood_cells")] <- "immune"
# cell_type_metadata$subtype[which(cell_type_metadata$subtype=="red_blood_cells")] <- "erythrocytes"

heat_tidy <- as_tibble(cbind(heat$data,cell_type_metadata[match(heat$data$Celltype_id,cell_type_metadata$subtype),"main_type"]))
colnames(heat_tidy)[ncol(heat_tidy)] <- "main_type"

#generate matrix of significance symbols for plot overlay
gwas = unique(heat_tidy$GWAS)
celltype = unique(heat_tidy$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy))
{
  sig_mat[heat_tidy$GWAS[i], heat_tidy$Celltype_id[i]] = heat_tidy$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,  .001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))


sig_mat_longer <- as.data.frame(heat_tidy)
sig_mat_longer <- sig_mat_longer[,which(colnames(sig_mat_longer)%in%c("GWAS","Celltype_id","P","main_type"))]

sig_counts_maintype <- sig_mat_longer %>%
  group_by(main_type) %>%
  summarize(count_sig = sum(P < 0.05))%>% 
  arrange(desc(count_sig))

sig_counts_subtype <- sig_mat_longer %>%
  group_by(Celltype_id) %>%
  summarize(count_sig = sum(P < 0.05)) %>% 
  arrange(desc(count_sig))

sig_counts_subtype$main_type <- sig_mat_longer$main_type[match(sig_counts_subtype$Celltype_id,sig_mat_longer$Celltype_id)]


# cleft_celltypes <- unique(sig_mat_longer$Celltype_id[intersect(grep("cleft",sig_mat_longer$GWAS,ignore.case = T),which(sig_mat_longer$P<0.05))])

annot_colors <-  hue_pal()(length(unique(heat_tidy$main_type)))
annot_colors <- list(main_type = annot_colors)
names(annot_colors$main_type) <- levels(heat_tidy$main_type)

# heat_tidy$Celltype_id[which(heat_tidy$Celltype_id=="other_blood_cells")] <- "immune"
# heat_tidy$Celltype_id[which(heat_tidy$Celltype_id=="red_blood_cells")] <- "erythrocytes"

p1 <- tidyheatmap(df = heat_tidy,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -2.5,
                  color_legend_max = 0,
                  # cutree_rows = 12,
                  # cutree_cols = 5,
                  angle_col = 315,
                  fontsize_number = 14,
                  fontsize_col = 12,
                  fontsize_row = 12,
                  border_color = "grey",
                  display_numbers = sig_mat,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  colors = c("darkblue", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_subtype.pdf", h = 18, w = 32)
p1
dev.off()

p1 <- tidyheatmap(df = heat_tidy[-which(heat_tidy$GWAS%in%c("CHRONLARGE","SLE_FG")),],
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -2.5,
                  color_legend_max = 0,
                  # cutree_rows = 12,
                  # cutree_cols = 5,
                  angle_col = 315,
                  fontsize_number = 14,
                  fontsize_col = 12,
                  fontsize_row = 12,
                  border_color = "grey",
                  display_numbers = sig_mat[-which(gwas%in%c("CHRONLARGE","SLE_FG")),],
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  colors = c("darkblue", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_subtype_NoChrons.pdf", h = 18, w = 32)
p1
dev.off()


finngen_subtype_data <- heat_tidy

sig_df <- data.frame(matrix(as.character(sig_mat), nrow = 45, ncol = 89))
rownames(sig_df) <- rownames(sig_mat)
colnames(sig_df) <- unique(heat_tidy$Celltype_id)
sig_df <-melt(as.matrix(sig_df))

celltype_colors <- hue_pal()(8)
names(celltype_colors) <- levels(heat_tidy$main_type)

GWAS_to_plot <- unique(finngen_subtype_data$GWAS)
p <- list()
p_top20 <- list()
for (name in GWAS_to_plot){
  data.to.plot <- finngen_subtype_data[which(finngen_subtype_data$GWAS==name),]
  data.to.plot <- data.to.plot[order(data.to.plot$log10p),]
  data.to.plot$Celltype_id <- factor(data.to.plot$Celltype_id,levels=data.to.plot$Celltype_id)
  
  sig <- sig_df[which(sig_df$Var1==name),]
  data.to.plot$significance <- ""
  data.to.plot$significance[match(sig$Var2,data.to.plot$Celltype_id)] <- sig$value

p[[name]] <- ggplot(data.to.plot,
       aes(x=Celltype_id,y=-log10p,fill=main_type))+
  geom_bar(stat = "identity") +  # `stat = "identity"` ensures bars reflect values
  geom_text(aes(label = significance, y = -log10p),
            size = 6, vjust = -0.25) +  # Adjust y to place text above bars
  geom_hline(yintercept = -log10(0.05),linetype=2,linewidth = 0.25) +
  labs(title = name, x = "Cell Subtype", y = "-log10p") +
  theme_minimal() + theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size=8),  # Angles X-axis text
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color="black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

p_top20[[name]] <- ggplot(data.to.plot[1:20,],
                    aes(x=Celltype_id,y=-log10p,fill=main_type))+
  geom_bar(stat = "identity") +  # `stat = "identity"` ensures bars reflect values
  geom_text(aes(label = significance, y = -log10p),
            size = 6, vjust = -0.15) +  # Adjust y to place text above bars
  geom_hline(yintercept = -log10(0.05),linetype=2,linewidth = 0.25) +
  labs(title = name, x = "Cell Subtype", y = "-log10p") +
  scale_fill_manual(values=celltype_colors) +
  theme_minimal() + theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size=8),  # Angles X-axis text
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color="black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
}

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_subtype_barplots.pdf",height=8.5,width=18)
p
dev.off()

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_subtype_top20_barplots.pdf",height=6,width=8)
p_top20
dev.off()

#generate matrix of significance symbols for plot overlay
heat_tidy_2 <- as_tibble(heat_2$data)
gwas = unique(heat_tidy_2$GWAS)
celltype = unique(heat_tidy_2$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_2))
{
  sig_mat[heat_tidy_2$GWAS[i], heat_tidy_2$Celltype_id[i]] = heat_tidy_2$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,  .001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))


heat_tidy_2$Celltype_id[which(heat_tidy_2$Celltype_id=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_2$Celltype_id[which(heat_tidy_2$Celltype_id=="SOX2_SOX10+")] <- "SOX2-_SOX10+"

# heat_tidy_2$Celltype[which(heat_tidy_2$Celltype=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
# heat_tidy_2$Celltype[which(heat_tidy_2$Celltype=="SOX2_SOX10+")] <- "SOX2-_SOX10+"

heat_tidy_2$Celltype <- heat_tidy_2$Celltype_id

p2 <- tidyheatmap(df = heat_tidy_2,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = Celltype,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -2.5,
                  color_legend_max = 0,
                  # cutree_rows = 10,
                  # cutree_cols = 4,
                  fontsize_number = 12,
                  border_color = "grey",
                  display_numbers = sig_mat,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  colors = c("darkblue", "white"),
                  annotation_colors = annot_colors)



pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_main.pdf", h = 8.5, w = 11)
p2
dev.off()


sig_df <- data.frame(matrix(as.character(sig_mat), nrow = 47, ncol = 8))
rownames(sig_df) <- rownames(sig_mat)
colnames(sig_df) <- unique(heat_tidy_2$Celltype_id)
sig_df <-melt(as.matrix(sig_df))


GWAS_to_plot <-  unique(heat_tidy_2$GWAS)
p <- list()
for (name in GWAS_to_plot){
  data.to.plot <- heat_tidy_2[which(heat_tidy_2$GWAS==name),]
  data.to.plot$fill <- factor(data.to.plot$Celltype_id,levels=data.to.plot$Celltype_id)
  data.to.plot <- data.to.plot[order(data.to.plot$log10p),]
  data.to.plot$Celltype_id <- factor(data.to.plot$Celltype_id,levels=data.to.plot$Celltype_id)
  
  sig <- sig_df[which(sig_df$Var1==name),]
  data.to.plot$significance <- ""
  data.to.plot$significance[match(sig$Var2,data.to.plot$Celltype_id)] <- sig$value
  
  p[[name]] <- ggplot(data.to.plot,
                      aes(x=Celltype_id,y=-log10p,fill=fill))+
    geom_bar(stat = "identity") +  # `stat = "identity"` ensures bars reflect values
    geom_text(aes(label = significance, y = -log10p),
              size = 6, vjust = -0.1) +  # Adjust y to place text above bars
    geom_hline(yintercept = -log10(0.05),linetype=2,linewidth = 0.25) +
    labs(title = name, x = "Cell Type", y = "-log10p") +
    theme_minimal() + theme(
      axis.text.x = element_text(angle = 45, hjust = 1,size=8),  # Angles X-axis text
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color="black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
}

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_main_barplots.pdf",height=6,width=8)
p
dev.off()

## for subtype_noprog
heat_3$data$Celltype_id[which(heat_3$data$Celltype_id=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_3 <- as_tibble(cbind(heat_3$data,cell_type_metadata[match(heat_3$data$Celltype_id,cell_type_metadata$subtype_noprog),"main_type"]))
colnames(heat_tidy_3)[ncol(heat_tidy_3)] <- "main_type"
gwas = unique(heat_tidy_3$GWAS)
celltype = unique(heat_tidy_3$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_3))
{
  sig_mat[heat_tidy_3$GWAS[i], heat_tidy_3$Celltype_id[i]] = heat_tidy_3$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,  .001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))


heat_tidy_3$Celltype_id[which(heat_tidy_3$Celltype_id=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"

heat_tidy_3$Celltype <- heat_tidy_3$Celltype_id

annot_colors <-  hue_pal()(length(unique(heat_tidy_3$main_type)))
annot_colors <- list(main_type = annot_colors)
names(annot_colors$main_type) <- levels(heat_tidy_3$main_type)

p3 <- tidyheatmap(df = heat_tidy_3,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -2.5,
                  color_legend_max = 0,
                  # cutree_rows = 10,
                  # cutree_cols = 4,
                  fontsize_number = 12,
                  border_color = "grey",
                  display_numbers = sig_mat,
                  # angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  colors = c("darkblue", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_subtype-noprog.pdf", h = 10, w = 20)
p3
dev.off()

p3 <- tidyheatmap(df = heat_tidy_3[-which(heat_tidy_3$GWAS%in%c("CHRONLARGE","SLE_FG")),],
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -2.5,
                  color_legend_max = 0,
                  # cutree_rows = 10,
                  # cutree_cols = 4,
                  fontsize_number = 12,
                  border_color = "grey",
                  display_numbers = sig_mat[-which(gwas%in%c("CHRONLARGE","SLE_FG")),],
                  # angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  colors = c("darkblue", "white"),
                  annotation_colors = annot_colors)



pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_subtype-noprog-NoChron.pdf", h = 10, w = 20)
p3
dev.off()

GWAS_subset <- unique(heat_tidy_3$GWAS[which(heat_tidy_3$P<0.05)])
heat_tidy_3_subset <- heat_tidy_3[which(heat_tidy_3$GWAS%in%GWAS_subset),]

gwas = unique(heat_tidy_3_subset$GWAS)
celltype = unique(heat_tidy_3_subset$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_3_subset))
{
  sig_mat[heat_tidy_3_subset$GWAS[i], heat_tidy_3_subset$Celltype_id[i]] = heat_tidy_3_subset$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,  .001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))


p4 <- tidyheatmap(df = heat_tidy_3_subset[-which(heat_tidy_3_subset$GWAS%in%c("CHRONLARGE","SLE_FG")),],
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -2.5,
                  color_legend_max = 0,
                  # cutree_rows = 10,
                  # cutree_cols = 4,
                  fontsize_number = 12,
                  border_color = "grey",
                  display_numbers = sig_mat[-which(gwas%in%c("CHRONLARGE","SLE_FG")),],
                  # angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  colors = c("darkblue", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_finngen_subtype-noprog-sigonly-NoChron.pdf", h = 8, w = 22)
p4
dev.off()

## 


MAGMA_results_variation_combined <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs =   magma_dirs_all_varition,
  ctd = face_subtype_ctd,
  ctd_species = "human", 
  ctd_name = "face_subtype_ctd",
  upstream_kb = 100,
  downstream_kb = 100,
  nThread = 8,
  force_new = TRUE,
  run_linear = TRUE,
  run_conditional = FALSE)

# names(MAGMA_results_finngen_subtypes) <- gsub("finngen_R11_","", names(MAGMA_results_finngen_subtypes))
# names(MAGMA_results_finngen_subtypes) <- gsub("Q17_","", names(MAGMA_results_finngen_subtypes))
# names(MAGMA_results_variation_combined) <- gsub(".formatted.tsv.bgz","", names(MAGMA_results_variation_combined))
saveRDS(MAGMA_results_variation_combined, "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/MAGMA_results_variation_combined.rds")

MAGMA_results_variation_combined <- readRDS(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/MAGMA_results_variation_combined.rds")



merged_res <- MAGMA.Celltyping::merge_results(MAGMA_results_variation_combined, level = 1, filetype = "ctAssocsLinear")

heat <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res,
                                          fdr_thresh = 1, fill_var = "-log10p")
merged_res_2 <- MAGMA.Celltyping::merge_results(MAGMA_results_variation_combined, level = 2, filetype = "ctAssocsLinear")

heat_2 <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res_2,
                                            fdr_thresh = 1, fill_var = "-log10p")

heat$data$GWAS <- gsub("XiongZ_31763980.","", heat$data$GWAS)
heat$data$GWAS <- gsub(".formatted.tsv.bgz","", heat$data$GWAS)

heat_2$data$GWAS <- gsub("XiongZ_31763980.","", heat_2$data$GWAS)
heat_2$data$GWAS <- gsub(".formatted.tsv.bgz","", heat_2$data$GWAS)

heat_tidy <- as.tibble(cbind(heat$data,cell_type_metadata[match(heat$data$Celltype_id,cell_type_metadata$subtype),"main_type"]))
colnames(heat_tidy)[ncol(heat_tidy)] <- "main_type"

variation_gwas_metadata <- read.table("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/variation_gwas_metadata.txt", header = T, sep = "\t")

heat_tidy <- as.tibble(cbind(heat_tidy,variation_gwas_metadata[match(heat_tidy$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))

heat_tidy_2 <- as.tibble(cbind(heat_2$data,cell_type_metadata[match(heat_2$data$Celltype_id,cell_type_metadata$subtype),c(3)]))
heat_tidy_2 <- as.tibble(cbind(heat_tidy_2,variation_gwas_metadata[match(heat_tidy_2$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))



# heat_tidy$Celltype_id[which(heat_tidy$Celltype_id=="other_blood_cells")] <- "immune"
# heat_tidy$Celltype_id[which(heat_tidy$Celltype_id=="red_blood_cells")] <- "erythrocytes"

p1 <- tidyheatmap(df = heat_tidy,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  clustering_distance_rows = "manhattan",
                  clustering_distance_cols = "manhattan",
                  clustering_method = "complete",
                  colors = c("#950606", "white"),
                  annotation_colors = annot_colors)

heat_tidy_2$Celltype_id[which(heat_tidy_2$Celltype_id=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_2$Celltype_id[which(heat_tidy_2$Celltype_id=="SOX2_SOX10+")] <- "SOX2-_SOX10+"

heat_tidy_2$Celltype <- gsub(" ","",heat_tidy_2$Celltype)
heat_tidy_2$Celltype[which(heat_tidy_2$Celltype=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_2$Celltype[which(heat_tidy_2$Celltype=="SOX2_SOX10+")] <- "SOX2-_SOX10+"

p2 <- tidyheatmap(df = heat_tidy_2,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = Celltype,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  clustering_distance_rows = "manhattan",
                  clustering_distance_cols = "manhattan",
                  clustering_method = "complete",
                  colors = c("#950606", "white"),
                  annotation_colors = annot_colors)


gwas = unique(heat_tidy$GWAS)
celltype = unique(heat_tidy$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy))
{
  sig_mat[heat_tidy$GWAS[i], heat_tidy$Celltype_id[i]] = heat_tidy$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,0.001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))

sig_mat_longer <- as.data.frame(heat_tidy)
sig_mat_longer <- sig_mat_longer[,which(colnames(sig_mat_longer)%in%c("GWAS","Celltype_id","P","main_type"))]

sig_counts_maintype <- sig_mat_longer %>%
  group_by(main_type) %>%
  summarize(count_sig = sum(P < 0.05))%>% 
  arrange(desc(count_sig))

sig_counts_subtype <- sig_mat_longer %>%
  group_by(Celltype_id) %>%
  summarize(count_sig = sum(P < 0.05)) %>% 
  arrange(desc(count_sig))

sig_counts_subtype$main_type <- sig_mat_longer$main_type[match(sig_counts_subtype$Celltype_id,sig_mat_longer$Celltype_id)]

p1 <- tidyheatmap(df = heat_tidy,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  fontsize_number = 12,
                  fontsize = 12,
                  display_numbers = sig_mat,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("#950606", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_combined_gwas_subtype.pdf", w = 28, h = 14)
p1
dev.off()


sig_df <- data.frame(matrix(as.character(sig_mat), nrow = 57, ncol = 89))
rownames(sig_df) <- rownames(sig_mat)
colnames(sig_df) <- unique(heat_tidy$Celltype_id)
sig_df <-melt(as.matrix(sig_df))


GWAS_to_plot <-  unique(heat_tidy$GWAS)
p <- list()
p_top20 <- list()
for (name in GWAS_to_plot){
  data.to.plot <- heat_tidy[which(heat_tidy$GWAS==name),]
  data.to.plot$fill <- data.to.plot$main_type
  data.to.plot <- data.to.plot[order(data.to.plot$log10p),]
  data.to.plot$Celltype_id <- factor(data.to.plot$Celltype_id,levels=data.to.plot$Celltype_id)
  
  sig <- sig_df[which(sig_df$Var1==name),]
  data.to.plot$significance <- ""
  data.to.plot$significance[match(sig$Var2,data.to.plot$Celltype_id)] <- sig$value
  
  p[[name]] <- ggplot(data.to.plot,
                      aes(x=Celltype_id,y=-log10p,fill=fill))+
    geom_bar(stat = "identity") +  # `stat = "identity"` ensures bars reflect values
    geom_text(aes(label = significance, y = -log10p),
              size = 6, vjust = -0.1) +  # Adjust y to place text above bars
    geom_hline(yintercept = -log10(0.05),linetype=2,linewidth = 0.25) +
    labs(title = name, x = "Cell Type", y = "-log10p") +
    theme_minimal() + theme(
      axis.text.x = element_text(angle = 45, hjust = 1,size=8),  # Angles X-axis text
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color="black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  p_top20[[name]] <- ggplot(data.to.plot[1:20,],
                            aes(x=Celltype_id,y=-log10p,fill=fill))+
    geom_bar(stat = "identity") +  # `stat = "identity"` ensures bars reflect values
    geom_text(aes(label = significance, y = -log10p),
              size = 6, vjust = -0.15) +  # Adjust y to place text above bars
    geom_hline(yintercept = -log10(0.05),linetype=2,linewidth = 0.25) +
    labs(title = name, x = "Cell Subtype", y = "-log10p") +
    scale_fill_manual(values=celltype_colors) +
    theme_minimal() + theme(
      axis.text.x = element_text(angle = 45, hjust = 1,size=8),  # Angles X-axis text
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color="black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
}

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_combined_gwas_subtype_barplots.pdf",height=14,width=28)
p
dev.off()

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_combined_gwas_subtype_top20_barplots.pdf",height=6,width=8)
p_top20
dev.off()


#generate matrix of significance symbols for plot overlay
gwas = unique(heat_tidy_2$GWAS)
celltype = unique(heat_tidy_2$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_2))
{
  sig_mat[heat_tidy_2$GWAS[i], heat_tidy_2$Celltype_id[i]] = heat_tidy_2$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,  .001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))

heat_tidy_2$Celltype <- gsub(" ","",heat_tidy_2$Celltype)
p2 <- tidyheatmap(df = heat_tidy_2,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = Celltype,
                  annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  display_numbers = sig_mat,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("#950606", "white"),
                  annotation_colors = annot_colors)


pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_combined_gwas_main.pdf", w = 8.5, h = 11)
p2
dev.off()

sig_df <- data.frame(matrix(as.character(sig_mat), nrow = 57, ncol = 8))
rownames(sig_df) <- rownames(sig_mat)
colnames(sig_df) <- unique(heat_tidy_2$Celltype_id)
sig_df <-melt(as.matrix(sig_df))


GWAS_to_plot <-  unique(heat_tidy_2$GWAS)
p <- list()
for (name in GWAS_to_plot){
  data.to.plot <- heat_tidy_2[which(heat_tidy_2$GWAS==name),]
  data.to.plot$fill <- factor(data.to.plot$Celltype_id,levels=levels(initial.filt_singlet_noneuro$celltype))
  data.to.plot <- data.to.plot[order(data.to.plot$log10p),]
  data.to.plot$Celltype_id <- factor(data.to.plot$Celltype_id,levels=data.to.plot$Celltype_id)
  
  sig <- sig_df[which(sig_df$Var1==name),]
  data.to.plot$significance <- ""
  data.to.plot$significance[match(sig$Var2,data.to.plot$Celltype_id)] <- sig$value
  
  p[[name]] <- ggplot(data.to.plot,
                      aes(x=Celltype_id,y=-log10p,fill=fill))+
    geom_bar(stat = "identity") +  # `stat = "identity"` ensures bars reflect values
    geom_text(aes(label = significance, y = -log10p),
              size = 6, vjust = -0.1) +  # Adjust y to place text above bars
    geom_hline(yintercept = -log10(0.05),linetype=2,linewidth = 0.25) +
    labs(title = name, x = "Cell Type", y = "-log10p") +
    theme_minimal() + theme(
      axis.text.x = element_text(angle = 45, hjust = 1,size=8),  # Angles X-axis text
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color="black"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
}

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_combined_gwas_main_barplots.pdf",height=14,width=28)
p
dev.off()




###
## xiong facial variation only ###

MAGMA_results_xiong <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs =   magma_dirs_xiong,
  ctd = face_subtype_ctd,
  ctd_species = "human", 
  ctd_name = "face_subtype_ctd",
  upstream_kb = 100,
  downstream_kb = 100,
  nThread = 8,
  force_new = TRUE,
  run_linear = TRUE,
  run_conditional = FALSE)

# saveRDS(MAGMA_results_xiong, "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/MAGMA_results_xiong_only.rds")
# 
MAGMA_results_xiong <- readRDS(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/MAGMA_results_xiong_only.rds")


merged_res <- MAGMA.Celltyping::merge_results(MAGMA_results_xiong, level = 1, filetype = "ctAssocsLinear")

heat <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res,
                                          fdr_thresh = 1, fill_var = "-log10p")
merged_res_2 <- MAGMA.Celltyping::merge_results(MAGMA_results_xiong, level = 2, filetype = "ctAssocsLinear")

heat_2 <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res_2,
                                            fdr_thresh = 1, fill_var = "-log10p")

merged_res_3 <- MAGMA.Celltyping::merge_results(MAGMA_results_xiong, level = 3, filetype = "ctAssocsLinear")

heat_3 <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res_3,
                                            fdr_thresh = 1, fill_var = "-log10p")

heat$data$GWAS <- gsub("XiongZ_31763980.","", heat$data$GWAS)
heat$data$GWAS <- gsub(".formatted.tsv.bgz","", heat$data$GWAS)

heat_2$data$GWAS <- gsub("XiongZ_31763980.","", heat_2$data$GWAS)
heat_2$data$GWAS <- gsub(".formatted.tsv.bgz","", heat_2$data$GWAS)

heat_3$data$GWAS <- gsub("XiongZ_31763980.","", heat_3$data$GWAS)
heat_3$data$GWAS <- gsub(".formatted.tsv.bgz","", heat_3$data$GWAS)

heat_tidy <- as.tibble(cbind(heat$data,cell_type_metadata[match(heat$data$Celltype_id,cell_type_metadata$subtype),"main_type"]))
colnames(heat_tidy)[ncol(heat_tidy)] <- "main_type"

variation_gwas_metadata <- read.table("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/variation_gwas_metadata.txt", header = T, sep = "\t")

heat_tidy <- as.tibble(cbind(heat_tidy,variation_gwas_metadata[match(heat_tidy$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))

heat_tidy_2 <- as.tibble(cbind(heat_2$data,cell_type_metadata[match(heat_2$data$Celltype_id,cell_type_metadata$subtype),c(3)]))
heat_tidy_2 <- as.tibble(cbind(heat_tidy_2,variation_gwas_metadata[match(heat_tidy_2$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))

heat_3$data$Celltype_id[which(heat_3$data$Celltype_id=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_3 <- as.tibble(cbind(heat_3$data,cell_type_metadata[match(heat_3$data$Celltype_id,cell_type_metadata$subtype_noprog),"main_type"]))
colnames(heat_tidy_3)[ncol(heat_tidy_3)] <- "main_type"
heat_tidy_3 <- as.tibble(cbind(heat_tidy_3,variation_gwas_metadata[match(heat_tidy_3$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))


gwas = unique(heat_tidy$GWAS)
celltype = unique(heat_tidy$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy))
{
  sig_mat[heat_tidy$GWAS[i], heat_tidy$Celltype_id[i]] = heat_tidy$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,0.001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))

sig_mat_longer <- as.data.frame(heat_tidy)
sig_mat_longer <- sig_mat_longer[,which(colnames(sig_mat_longer)%in%c("GWAS","Celltype_id","P","main_type"))]

sig_counts_maintype <- sig_mat_longer %>%
  group_by(main_type) %>%
  summarize(count_sig = sum(P < 0.05))%>% 
  arrange(desc(count_sig))

sig_counts_subtype <- sig_mat_longer %>%
  group_by(Celltype_id) %>%
  summarize(count_sig = sum(P < 0.05)) %>% 
  arrange(desc(count_sig))

sig_counts_subtype$main_type <- sig_mat_longer$main_type[match(sig_counts_subtype$Celltype_id,sig_mat_longer$Celltype_id)]

p1 <- tidyheatmap(df = heat_tidy,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  # annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  fontsize_number = 12,
                  fontsize = 12,
                  display_numbers = sig_mat,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("darkgreen", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_xiong-only_gwas_subtype.pdf", w = 28, h = 14)
p1
dev.off()


heat_tidy_2$Celltype_id[which(heat_tidy_2$Celltype_id=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_2$Celltype_id[which(heat_tidy_2$Celltype_id=="SOX2_SOX10+")] <- "SOX2-_SOX10+"

heat_tidy_2$Celltype <- gsub(" ","",heat_tidy_2$Celltype)
heat_tidy_2$Celltype[which(heat_tidy_2$Celltype=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_2$Celltype[which(heat_tidy_2$Celltype=="SOX2_SOX10+")] <- "SOX2-_SOX10+"

#generate matrix of significance symbols for plot overlay
gwas = unique(heat_tidy_2$GWAS)
celltype = unique(heat_tidy_2$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_2))
{
  sig_mat[heat_tidy_2$GWAS[i], heat_tidy_2$Celltype_id[i]] = heat_tidy_2$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,  .001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))

heat_tidy_2$Celltype <- gsub(" ","",heat_tidy_2$Celltype)
p2 <- tidyheatmap(df = heat_tidy_2,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = Celltype,
                  annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  display_numbers = sig_mat,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("darkgreen", "white"),
                  annotation_colors = annot_colors)


pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_xiong-only_gwas_main.pdf", w = 8.5, h = 11)
p2
dev.off()

heat_tidy_3 <- heat_tidy_3[-grep(".h",heat_tidy_3$GWAS),]
gwas = unique(heat_tidy_3$GWAS)
celltype = unique(heat_tidy_3$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_3))
{
  sig_mat[heat_tidy_3$GWAS[i], heat_tidy_3$Celltype_id[i]] = heat_tidy_3$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,0.001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))


p3 <- tidyheatmap(df = heat_tidy_3,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  # annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  fontsize_number = 12,
                  fontsize = 12,
                  display_numbers = sig_mat,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("darkgreen", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_xiong-only_gwas_subtype-noprog.pdf", w = 28, h = 14)
p3
dev.off()


GWAS_subset <- unique(heat_tidy_3$GWAS[which(heat_tidy_3$P<0.05)])
heat_tidy_3_subset <- heat_tidy_3[which(heat_tidy_3$GWAS%in%GWAS_subset),]
gwas = unique(heat_tidy_3_subset$GWAS)
celltype = unique(heat_tidy_3_subset$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_3_subset))
{
  sig_mat[heat_tidy_3_subset$GWAS[i], heat_tidy_3_subset$Celltype_id[i]] = heat_tidy_3_subset$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,0.001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))


p4 <- tidyheatmap(df = heat_tidy_3_subset,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  # annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  fontsize_number = 12,
                  fontsize = 12,
                  display_numbers = sig_mat,
                  # angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("darkgreen", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_xiong-only_gwas_subtype-sigonly.pdf", w = 28, h = 12)
p4
dev.off()

##
## bonfontae facial variation only ###

MAGMA_results_bonfontae <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs =   magma_dirs_bonfontae,
  ctd = face_subtype_ctd,
  ctd_species = "human", 
  ctd_name = "face_subtype_ctd",
  upstream_kb = 100,
  downstream_kb = 100,
  nThread = 8,
  force_new = TRUE,
  run_linear = TRUE,
  run_conditional = FALSE)

# saveRDS(MAGMA_results_bonfontae, "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/MAGMA_results_bonfontae_only.rds")
# 
MAGMA_results_bonfontae <- readRDS(file = "/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/MAGMA/MAGMA_results_bonfontae_only.rds")


merged_res <- MAGMA.Celltyping::merge_results(MAGMA_results_bonfontae, level = 1, filetype = "ctAssocsLinear")

heat <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res,
                                          fdr_thresh = 1, fill_var = "-log10p")
merged_res_2 <- MAGMA.Celltyping::merge_results(MAGMA_results_bonfontae, level = 2, filetype = "ctAssocsLinear")

heat_2 <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res_2,
                                            fdr_thresh = 1, fill_var = "-log10p")

merged_res_3 <- MAGMA.Celltyping::merge_results(MAGMA_results_bonfontae, level = 3, filetype = "ctAssocsLinear")

heat_3 <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res_3,
                                            fdr_thresh = 1, fill_var = "-log10p")

heat$data$GWAS <- gsub(".formatted.tsv.bgz","", heat$data$GWAS)
heat_2$data$GWAS <- gsub(".formatted.tsv.bgz","", heat_2$data$GWAS)
heat_3$data$GWAS <- gsub(".formatted.tsv.bgz","", heat_3$data$GWAS)

heat_tidy <- as.tibble(cbind(heat$data,cell_type_metadata[match(heat$data$Celltype_id,cell_type_metadata$subtype),"main_type"]))
colnames(heat_tidy)[ncol(heat_tidy)] <- "main_type"

variation_gwas_metadata <- read.table("/mnt/isilon/cotney_lab_hpc/cotneyj/scRNA-Seq_GWAS/variation_gwas_metadata.txt", header = T, sep = "\t")

heat_tidy <- as.tibble(cbind(heat_tidy,variation_gwas_metadata[match(heat_tidy$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))

heat_tidy_2 <- as.tibble(cbind(heat_2$data,cell_type_metadata[match(heat_2$data$Celltype_id,cell_type_metadata$subtype),c(3)]))
heat_tidy_2 <- as.tibble(cbind(heat_tidy_2,variation_gwas_metadata[match(heat_tidy_2$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))

heat_3$data$Celltype_id[which(heat_3$data$Celltype_id=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_3 <- as.tibble(cbind(heat_3$data,cell_type_metadata[match(heat_3$data$Celltype_id,cell_type_metadata$subtype_noprog),"main_type"]))
colnames(heat_tidy_3)[ncol(heat_tidy_3)] <- "main_type"
heat_tidy_3 <- as.tibble(cbind(heat_tidy_3,variation_gwas_metadata[match(heat_tidy_3$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))



heat_tidy_2$Celltype_id[which(heat_tidy_2$Celltype_id=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_2$Celltype_id[which(heat_tidy_2$Celltype_id=="SOX2_SOX10+")] <- "SOX2-_SOX10+"

heat_tidy_2$Celltype <- gsub(" ","",heat_tidy_2$Celltype)
heat_tidy_2$Celltype[which(heat_tidy_2$Celltype=="SOX2+_SOX10_")] <- "SOX2+_SOX10-"
heat_tidy_2$Celltype[which(heat_tidy_2$Celltype=="SOX2_SOX10+")] <- "SOX2-_SOX10+"


gwas = unique(heat_tidy$GWAS)
celltype = unique(heat_tidy$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy))
{
  sig_mat[heat_tidy$GWAS[i], heat_tidy$Celltype_id[i]] = heat_tidy$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,0.001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))

sig_mat_longer <- as.data.frame(heat_tidy)
sig_mat_longer <- sig_mat_longer[,which(colnames(sig_mat_longer)%in%c("GWAS","Celltype_id","P","main_type"))]

sig_counts_maintype <- sig_mat_longer %>%
  group_by(main_type) %>%
  summarize(count_sig = sum(P < 0.05))%>% 
  arrange(desc(count_sig))

sig_counts_subtype <- sig_mat_longer %>%
  group_by(Celltype_id) %>%
  summarize(count_sig = sum(P < 0.05)) %>% 
  arrange(desc(count_sig))

sig_counts_subtype$main_type <- sig_mat_longer$main_type[match(sig_counts_subtype$Celltype_id,sig_mat_longer$Celltype_id)]

p1 <- tidyheatmap(df = heat_tidy,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  # annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  fontsize_number = 12,
                  fontsize = 12,
                  display_numbers = sig_mat,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("#4B0063", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_bonfontae-only_gwas_subtype.pdf", w = 28, h = 14)
p1
dev.off()


#generate matrix of significance symbols for plot overlay
gwas = unique(heat_tidy_2$GWAS)
celltype = unique(heat_tidy_2$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_2))
{
  sig_mat[heat_tidy_2$GWAS[i], heat_tidy_2$Celltype_id[i]] = heat_tidy_2$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,  .001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))

heat_tidy_2$Celltype <- gsub(" ","",heat_tidy_2$Celltype)
p2 <- tidyheatmap(df = heat_tidy_2,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = Celltype,
                  annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  display_numbers = sig_mat,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("#950606", "white"),
                  annotation_colors = annot_colors)


pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_bonfontae-only_gwas_main.pdf", w = 8.5, h = 11)
p2
dev.off()

gwas = unique(heat_tidy_3$GWAS)
celltype = unique(heat_tidy_3$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_3))
{
  sig_mat[heat_tidy_3$GWAS[i], heat_tidy_3$Celltype_id[i]] = heat_tidy_3$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,0.001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))


p3 <- tidyheatmap(df = heat_tidy_3,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  fontsize_number = 12,
                  fontsize = 12,
                  display_numbers = sig_mat,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("#950606", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_bonfontae-only_gwas_subtype-noprog.pdf", w = 28, h = 14)
p3
dev.off()

GWAS_subset <- unique(heat_tidy_3$GWAS[which(heat_tidy_3$P<0.05)])
heat_tidy_3_subset <- heat_tidy_3[which(heat_tidy_3$GWAS%in%GWAS_subset),]
gwas = unique(heat_tidy_3_subset$GWAS)
celltype = unique(heat_tidy_3_subset$Celltype_id)

sig_mat = matrix(0,nrow = length(gwas), ncol = length(celltype))
rownames(sig_mat) = gwas
colnames(sig_mat) = celltype

for(i in 1:nrow(heat_tidy_3_subset))
{
  sig_mat[heat_tidy_3_subset$GWAS[i], heat_tidy_3_subset$Celltype_id[i]] = heat_tidy_3_subset$P[i]
}

sig_mat <- symnum(sig_mat, corr = FALSE,
                  cutpoints = c(0,0.001,.01,.05, .1, 1),
                  symbols = c("***","**","*","."," "))


p4 <- tidyheatmap(df = heat_tidy_3_subset,
                  rows = GWAS,
                  columns = Celltype_id,
                  values = log10p,
                  scale = "none",
                  annotation_col = main_type,
                  annotation_row = c("Region", "Type"),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  color_legend_n = 9,
                  color_legend_min = -3,
                  color_legend_max = 0,
                  fontsize_number = 12,
                  fontsize = 12,
                  display_numbers = sig_mat,
                  angle_col = 315,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "complete",
                  colors = c("#950606", "white"),
                  annotation_colors = annot_colors)

pdf("/scr1/users/manchela/Data/Craniofacial-scRNAseqManuscript/final-figs/craniofacial_variation_bonfontae-only_gwas_subtype-noprog-sigonly.pdf", w = 28, h = 12)
p4
dev.off()


