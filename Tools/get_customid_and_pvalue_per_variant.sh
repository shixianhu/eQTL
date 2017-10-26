#!\bin\bash
while IFS='' read line; do
	filename_part=$(echo ${line} | cut -d '_' -f 1)
	filename1="${filename_part}_doubleids_complete.txt"
	filename2="${filename_part}_doubleids.txt"
	filename="${filename_part}_doubleids.txt"
	grep -wf $line ../165test_passonly.noXnoY.hwe.nomiss.maf5percent_doubleids_cor.txt | awk '{print $1"\t"$2}' > "$filename"
	grep -f $line ../165test_clumpstrict_hwe_nomiss_maf5percent_allcor_FDR5percent_multi_ciseQTLs.txt | awk '{print $4}' > temp_file.txt
	cut -f 1 temp_file.txt | paste "$filename2" - > "$filename1"
done < filelist_eqtls.txt
