#filter studies for those with at least one significant association
#grep "genome-wide significant variants (P<5e-8)" finngen_*.mungesumstats_log/finngen_*_log_msg.txt | grep -v "0 genome-wide significant variants (P<5e-8)" | cut -f1 -d " " | uniq | sed -e 's/\./\t/g' -e 's/finngen_R11_//g' | cut -f1 > FINNGEN/significant_studies.txt

grep -f significant_studies.txt finngen_R11_manifest.tsv | cut -f1,4,5 | while read a b c; do
 echo "genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = \"/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_"$a".formatted.tsv.bgz\",
  genome_build = \"GRCh38\",
  N = "$(($b+$c))",
  population = \"eur\",
  upstream_kb = 100,
  downstream_kb = 100
  )" >> generate_finngen_celltype_associations_for_R.txt
  done

while read p; do
echo "\"/Users/cotneyj/Desktop/scRNA-Seq_GWAS/finngen_R11_"$p".formatted.tsv.bgz\"" >> magma_dirs.txt
done < significant_studies.txt

while read p; do
echo "ctAssocs"$p" <- MAGMA_results$\`finngen_R11_"$p".formatted.tsv.bgz\`\$ctAssocsLinear" >> ctAssocsLinear.txt
done < significant_studies.txt



