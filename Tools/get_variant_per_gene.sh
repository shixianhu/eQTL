#!\bin\bash
while IFS='' read line; do
	filename="${line}_eqtls.txt"
	grep $line 165test_clumpstrict_hwe_nomiss_maf5percent_allcor_FDR5percent_multi_ciseQTLs.txt | awk '{print $1}' > "./variant_ids/$filename"
done < geneids_toclump.txt
