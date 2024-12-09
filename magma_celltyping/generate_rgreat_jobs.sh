while read x; do
echo ""$x" <- \"hoskens_"$x".bed\"
bed."$x" <- import("$x", genome = \"hg19\")
job_"$x" = submitGreatJob(bed."$x",
                      genome = \"hg19\")
"$x".tbl = getEnrichmentTables(job_"$x")
"$x"_gene_regions = getRegionGeneAssociations(job_"$x", verbose = great_opt$verbose)
"$x"_genes <- unique(as.data.frame("$x"_gene_regions$annotated_genes)$value)
done < trait_list.txt
