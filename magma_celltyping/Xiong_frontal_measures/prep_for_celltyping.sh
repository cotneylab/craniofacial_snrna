while read p; do
 echo "genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = \"/Users/justincotney/Desktop/scRNA-Seq_GWAS/XiongZ_31763980."$p".formatted.tsv.bgz\",
  genome_build = \"GRCh38\",
  N = 10115,
  population = \"eur\",
  upstream_kb = 100,
  downstream_kb = 100
  )" >> generate_celltype_associations_for_R.txt
  done < significant_studies.txt

while read p; do
echo "\"/Users/justincotney/Desktop/scRNA-Seq_GWAS/XiongZ_31763980."$p".formatted.tsv.bgz\"" >> magma_dirs.txt
done < significant_studies.txt

while read p; do
echo "ctAssocs"$p" <- MAGMA_results$\`XiongZ_31763980."$p".formatted.tsv.bgz\`\$ctAssocsLinear" >> ctAssocsLinear.txt
done < significant_studies.txt