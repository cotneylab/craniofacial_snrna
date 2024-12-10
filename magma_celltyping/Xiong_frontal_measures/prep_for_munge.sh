while read p; do
 echo "reformatted_xiongz <- MungeSumstats::format_sumstats(
  path=\"XiongZ_31763980."$p".csv.gz\",
  save_path = \"XiongZ_31763980."$p".formatted.tsv.gz\",
  force_new = TRUE,
  ref_genome=\"GRCh37\",
  convert_ref_genome=\"GRCh38\",
  dbSNP = 155,
  allele_flip_check = TRUE,
  bi_allelic_filter = FALSE,
  write_vcf=FALSE,
  tabix_index=TRUE,
  nThread = 4,
  log_mungesumstats_msgs=TRUE,
  log_folder= \"XiongZ_31763980."$p".mungesumstats_log\",
  mapping_file = xiongz_gwas_colheader
)" >> generate_munge_for_R.txt
done < XiongZ_study_traits.txt