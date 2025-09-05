
#neurogenomicslab/magma.celltyping
Sys.setenv(GITHUB_PAT="your token")
install.packages("piggyback")
Sys.setenv('R_MAX_VSIZE'=64000000000)
if(!require("remotes")) install.packages("remotes")
library(MAGMA.Celltyping)
library(qqman)
face_ctd <- EWCE::load_rdata("ctd_face_ctd.rda") 

face_subtype_ctd <- EWCE::load_rdata("ctd_face_subtype_ctd.rda") 

storage_dir <- c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/MAGMA")

genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_ANKYLOGLOSSIA.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_APVR.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_ASD.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 449747,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_AND_CARIES.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_HARD_PALATE_NOS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453074,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_CLEFT_PALATE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_OR_LIP_AND_PALATE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453074,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_PALATE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453539,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_DEFORMITI_FEET.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451335,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_LENS_MALFO.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451518,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_ANTER_SEGMENT_EYE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451621,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_CIRCULATO_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_DEFORMAT_MUSCULOS_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_EAR_CAUSI_IMPAIRM_HEARING.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451480,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GALLB_BILE_DUCTS_LIVER.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452554,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GENITAL_ORGANS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GREAT_ARTERIES.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 449138,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GREAT_VEINS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 448543,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_NERVOUS_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_OVARIES_FALLOP_TUBES_BROAD_LIGAM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 253663,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_POSTERI_SEGMENT_EYE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451518,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_PULMONARY_TRICU_VALVES.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 448693,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_RETI.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451442,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_TRACHEA_BRONCHUS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453484,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_URINARY_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONOTR_DEFEC.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CYSTIC_KIDNEY_DISEA.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452843,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_LVOTO_BROAD.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_LVOTO_NARROW.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OSTEOCHONDROD_W_DEFECTS_GROWTH_TUBULAR_BONES_SPINE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450566,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CHROMOSOME_ABNORMALITI_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453278,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_BRAIN.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453157,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_CIRCULATO_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 449240,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_DIGES_SYSTEM1.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_DIGES_SYSTEM2.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452357,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_EAR.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451803,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_FACE_NECK.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452348,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_FEMALE_GENITALIA.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450456,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_HEART.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 448657,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_INTESTINE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452555,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_KIDNEY.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451936,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_LIMBS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450750,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_MALE_GENITAL_ORGANS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450370,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_NERVOUS_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453326,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450857,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_PERIP_VASCULAR_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 448862,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_RESPI_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453246,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_SKIN.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451959,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_SKULL_FACE_BONES.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450601,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_TONGUE_MOUTH_PHARYNX.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452740,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_UPPER_ALIME_TRACT.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 452503,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_URINARY_SYSTEM.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451965,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_OSTEOCHONDROD.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450779,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_SPECIFE_CONGEN_MALFO_SYNDR_AFFECTING_MULTIPLE_SYSTEMS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 451126,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_PHAKOMATOSES_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450945,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_RVOTO.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_SEPTA_DEFEC.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_UNDES_TESTICLE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 450490,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_RX_CROHN_1STLINE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_RX_CROHN_2NDLINE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 453733,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_SLE_FG.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 339484,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_CHRONLARGE.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 434250,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_CHRONNAS.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 434744,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_CHRONOTH.formatted.tsv.bgz",
  genome_build = "GRCh38",
  N = 434933,
  population = "eur",
  upstream_kb = 100,
  downstream_kb = 100
)
magma_dirs_finngen_immune = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_CHRONLARGE.formatted.tsv.bgz",
                              "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_SLE_FG.formatted.tsv.bgz")
magma_dirs_finngen_other = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_ANKYLOGLOSSIA.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_APVR.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_ASD.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_DEFORMITI_FEET.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_LENS_MALFO.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_ANTER_SEGMENT_EYE.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_CIRCULATO_SYSTEM.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_DEFORMAT_MUSCULOS_SYSTEM.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_EAR_CAUSI_IMPAIRM_HEARING.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GALLB_BILE_DUCTS_LIVER.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GENITAL_ORGANS.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GREAT_ARTERIES.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_GREAT_VEINS.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_NERVOUS_SYSTEM.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_OVARIES_FALLOP_TUBES_BROAD_LIGAM.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_POSTERI_SEGMENT_EYE.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_PULMONARY_TRICU_VALVES.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_RETI.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_TRACHEA_BRONCHUS.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONGEN_MALFO_URINARY_SYSTEM.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CONOTR_DEFEC.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CYSTIC_KIDNEY_DISEA.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_LVOTO_BROAD.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_LVOTO_NARROW.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OSTEOCHONDROD_W_DEFECTS_GROWTH_TUBULAR_BONES_SPINE.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CHROMOSOME_ABNORMALITI_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_DIGES_SYSTEM1.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_EAR.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_FACE_NECK.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_HEART.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_SKIN.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_TONGUE_MOUTH_PHARYNX.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_CONGEN_MALFO_UPPER_ALIME_TRACT.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_OSTEOCHONDROD.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_OTHER_SPECIFE_CONGEN_MALFO_SYNDR_AFFECTING_MULTIPLE_SYSTEMS.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_PHAKOMATOSES_NOT_ELSEW_CLASSIFIED.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_RVOTO.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_SEPTA_DEFEC.formatted.tsv.bgz",
                             "/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_UNDES_TESTICLE.formatted.tsv.bgz")


#spot check manhattan plots
library(data.table)
library(dplyr)
library(manhattanly)

finngen_cleft = fread("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_CLEFT_PALATE.formatted.tsv.bgz")
finngen_cleft.subset <- finngen_cleft %>% 
  filter(-log10(P)>1)

head(finngen_cleft.subset)

manhattan(finngen_cleft.subset, chr="CHR", bp="BP", snp="SNP", p="P", annotatePval = 0.0001)

#interactive version
manhattanly(finngen_cleft.subset, chr="CHR", bp="BP", snp="SNP", p="P", gene = "NEAREST_GENES")

library(locuszoomr)
library(LDlinkR)
#API token aba03b266910
token <- c("aba03b266910")


loc <- locus(data = finngen_cleft, chrom = "CHR", labs = "SNP", pos = "BP", p = "P", gene = 'IRF6', index_snp = "rs642961", flank = 1e5, ens_db = "EnsDb.Hsapiens.v86")
loc <- link_LD(loc, pop = "FIN", genome_build = "grch38", r2d = "r2", method = "proxy", token = token)
locus_plot(loc, label_x = c(4, -5), labels = c("index", "rs570516915"))

magma_dirs_finngen_cleft  = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_AND_CARIES.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_HARD_PALATE_NOS.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_CLEFT_PALATE.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_LIP_OR_LIP_AND_PALATE.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_Q17_CLEFT_PALATE.formatted.tsv.bgz")
magma_dirs_finngen = c(magma_dirs_finngen_cleft, magma_dirs_finngen_other,magma_dirs_finngen_immune)
magma_dirs_xiong = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.AlL-AlR.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.AlL-Ls.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.AlL-Sn.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.AlR-ChL.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.EnL-AlL.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.EnL-Sn.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.EnR-AlR.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.EnR-ChL.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.EnR-ChR.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.EnR-N.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.EnR-Prn.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.EnR-Sn.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.ExL-AlL.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.ExR-ChL.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.ExR-ChR.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.ExR-N.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.N-En.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.N-ExL.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.N-Prn.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.Prn-AlL.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.Prn-AlR.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.Prn-EnL.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.Prn-Ls.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/XiongZ_31763980.Sn-AlR.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/GCST011096.h.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/28067908-GCST004132-EFO_0000384.crohns.h.formatted.tsv.bgz")
magma_dirs_bonfontae = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Brow_Ridge_Protrusion_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Brow_Ridge_Protrusion_2.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Brow_Ridge_Protrusion_3.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Chin_Protrusion_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Chin_Protrusion_2.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Columella_Inclination.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Columella_Size.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Ear_Size.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Eye_Position_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Forehead_Protrusion_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Jaw_Protrusion_2.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Jaw_Protrusion_5.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Jaw_Slope_2.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Lip_Protrusion.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Lip_Thickness_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Lip_Thickness_Ratio_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Lip_Thickness_Ratio_2.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Lower_Face_Flatness.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Lower_Lip_Protrusion.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Lower_Lip_Thickness_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Lower_Lip_Thickness_2.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nasion_Depth_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nasion_Depth_2.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nasion_Position_2.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nose_Height.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nose_Protrusion.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nose_Roundness_1.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nose_Roundness_3.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nose_Size.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Nostril_Size.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Philtrum_Length.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/Upper_Face_Flatness.formatted.tsv.bgz")
magma_dirs_all_varition = c(magma_dirs_xiong, magma_dirs_bonfontae)
magma_dirs_nsclp = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/GMKF_CPseq_allpopulations_TDT_plinkSE_050124.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.meta.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/GCST011096.h.formatted.tsv.bgz")
#magma_dirs_nsclp = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/GMKF_CPseq_allpopulations_TDT_plinkSE_050124.com.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/GMKF_CPseq_allpopulations_TDT_plinkSE_050124.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/GMKF_CPseq_allpopulations_TDT_plinkSE_050124.par.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.beaty.tdt.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.marazita.cc.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.marazita.tdt.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.meta.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/pofc.geneva.meta.allpops.CLCLP.meta.update2018.meta.marazita.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/GCST011096.h.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/28067908-GCST004132-EFO_0000384.crohns.h.formatted.tsv.bgz")
magma_dirs_microsomia = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/craniofacial_microsomia.formatted.tsv.bgz")
magma_dirs_immune = c("/Users/cotneyj/Desktop/scRNA-Seq_GWAS/GCST011096.h.formatted.tsv.bgz","/Users/cotneyj/Desktop/scRNA-Seq_GWAS/28067908-GCST004132-EFO_0000384.crohns.h.formatted.tsv.bgz")
magma_dirs_craniofacial_disease =c(magma_dirs_nsclp,magma_dirs_microsomia,magma_dirs_immune)
magma_dirs_craniofacial_disease_variation = c(magma_dirs_nsclp,magma_dirs_microsomia,magma_dirs_xiong)
magma_dirs_craniofacial_disease_variation2 = c(magma_dirs_nsclp,magma_dirs_microsomia,magma_dirs_immune,magma_dirs_bonfontae)
test = c(magma_dirs_craniofacial_disease, magma_dirs_finngen_subset)


MAGMA_results_finngen_subtypes <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs =   magma_dirs_finngen,
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
saveRDS(MAGMA_results_finngen_subtypes, file = "MAGMA_results_finngen_subtypes.rds")

MAGMA_results_finngen_subtypes <- readRDS("MAGMA_results_finngen_subtypes.rds")

library(ggplot2)
library(tidyverse)
library(tidyheatmaps)

#subtype finngen analysis
merged_res <- MAGMA.Celltyping::merge_results(MAGMA_results_finngen_subtypes, level = 1, filetype = "ctAssocsLinear")

heat <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res,
                                          fdr_thresh = 1, fill_var = "-log10p")
merged_res_2 <- MAGMA.Celltyping::merge_results(MAGMA_results_finngen_subtypes, level = 2, filetype = "ctAssocsLinear")

heat_2 <- MAGMA.Celltyping::results_heatmap(merged_results = merged_res_2,
                                            fdr_thresh = 1, fill_var = "-log10p")

cell_type_metadata <- readxl::read_xlsx("annotation_expanded.xlsx", sheet = 1)
cell_type_metadata <- cell_type_metadata[,1:3]



heat_tidy <- as_tibble(cbind(heat$data,cell_type_metadata[match(heat$data$Celltype_id,cell_type_metadata$subtype),c(3)]))

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

heat_tidy_2 <- as_tibble(heat_2$data)

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
                  cutree_rows = 12,
                  cutree_cols = 5,
                  angle_col = 315,
                  fontsize_number = 12,
                  fontsize_col = 12,
                  border_color = "grey",
                  display_numbers = sig_mat,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  clustering_method = "ward.D2",
                  colors = c("#950606", "white"))

pdf("plots/craniofacial_finngen_subtype.pdf", h = 22, w = 17)
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
            cutree_rows = 10,
            cutree_cols = 4,
            fontsize_number = 12,
            border_color = "grey",
            display_numbers = sig_mat,
            angle_col = 315,
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean",
            clustering_method = "ward.D2",
            colors = c("#950606", "white"))



pdf("plots/craniofacial_finngen_main.pdf", h = 8.5, w = 11)
p2
dev.off()


MAGMA_results_variation_combined <- readRDS(file = "MAGMA_results_variation_combined.rds")



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

heat_tidy <- as.tibble(cbind(heat$data,cell_type_metadata[match(heat$data$Celltype_id,cell_type_metadata$subtype),c(3)]))

variation_gwas_metadata <- read.table("variation_gwas_metadata.txt", header = T, sep = "\t")

heat_tidy <- as.tibble(cbind(heat_tidy,variation_gwas_metadata[match(heat_tidy$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))

heat_tidy_2 <- as.tibble(cbind(heat_2$data,cell_type_metadata[match(heat_2$data$Celltype_id,cell_type_metadata$subtype),c(3)]))
heat_tidy_2 <- as.tibble(cbind(heat_tidy_2,variation_gwas_metadata[match(heat_tidy_2$GWAS,variation_gwas_metadata$Trait),c(2,3,4)]))





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
                  colors = c("#950606", "white"))

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
                  colors = c("#950606", "white"))


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
                  colors = c("#950606", "white"))

pdf("plots/craniofacial_variation_combined_gwas_subtype.pdf", w = 8.5, h = 11)
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
                  colors = c("#950606", "white"))


pdf("plots/craniofacial_variation_combined_gwas_main.pdf", w = 8.5, h = 11)
p2
dev.off()
