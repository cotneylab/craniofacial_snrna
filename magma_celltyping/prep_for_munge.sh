while read x; do
 echo "reformatted_finngen <- MungeSumstats::format_sumstats(
  path=\"FINNGEN/finngen_R11_"$x".gz\",
  save_path = \"finngen_R11_"$x".formatted.tsv.bgz\",
  force_new = TRUE,
  ref_genome=\"GRCh38\",
  convert_ref_genome=\"GRCh38\",
  local_chain = \"GRCh37_to_GRCh38.chain\",
  dbSNP = 155,
  allele_flip_check = TRUE,
  bi_allelic_filter = TRUE,
  write_vcf=FALSE,
  tabix_index=TRUE,
  nThread = 8,
  log_mungesumstats_msgs=TRUE,
  log_folder= \"finngen_R11_"$x".mungesumstats_log\"
  )
  " >> generate_munge_for_R.txt
done < congenital_studies.txt

while read x; do
 echo "reformatted_finngen <- MungeSumstats::format_sumstats(
  path=\"FINNGEN/finngen_R11_"$x".gz\",
  save_path = \"finngen_R11_"$x".formatted.tsv.bgz\",
  force_new = TRUE,
  ref_genome=\"GRCh38\",
  convert_ref_genome=\"GRCh38\",
  local_chain = \"GRCh37_to_GRCh38.chain\",
  dbSNP = 155,
  allele_flip_check = TRUE,
  bi_allelic_filter = TRUE,
  write_vcf=FALSE,
  tabix_index=TRUE,
  nThread = 8,
  log_mungesumstats_msgs=TRUE,
  log_folder= \"finngen_R11_"$x".mungesumstats_log\"
  )
  " >> generate_munge_for_R.txt
done < immune_studies.txt