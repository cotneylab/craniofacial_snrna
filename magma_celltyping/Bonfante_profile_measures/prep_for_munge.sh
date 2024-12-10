while read x y; do
 echo "reformatted_ruiz <- MungeSumstats::format_sumstats(
  path=\"RUIZ/"$x".tsv\",
  save_path = \""$y".formatted.tsv.bgz\",
  force_new = TRUE,
  ref_genome=\"GRCh37\",
  convert_ref_genome=\"GRCh38\",
  local_chain = \"GRCh37_to_GRCh38.chain\",
  dbSNP = 155,
  allele_flip_check = TRUE,
  bi_allelic_filter = TRUE,
  write_vcf=FALSE,
  tabix_index=TRUE,
  nThread = 8,
  log_mungesumstats_msgs=TRUE,
  log_folder= \""$y".mungesumstats_log\"
  )
  " >> generate_munge_for_R.txt
done < trait_names.txt
