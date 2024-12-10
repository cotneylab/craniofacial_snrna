while read p x; do
 echo "genesOut <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = \"/Users/cotneyj/Desktop/scRNA-Seq_GWAS/"$x".formatted.tsv.bgz\",
  genome_build = \"GRCh38\",
  N = 6169,
  population = \"amr\",
  upstream_kb = 100,
  downstream_kb = 100
  )" >> generate_ruiz_celltype_associations_for_R.txt
  done < trait_names.txt

while read p x; do
echo "\"/Users/cotneyj/Desktop/scRNA-Seq_GWAS/"$p".formatted.tsv.bgz\"" >> magma_dirs.txt
done < significant_studies.txt

while read p x; do
echo "ctAssocs"$p" <- MAGMA_results$\`"$p".formatted.tsv.bgz\`\$ctAssocsLinear" >> ctAssocsLinear.txt
done < significant_studies.txt
