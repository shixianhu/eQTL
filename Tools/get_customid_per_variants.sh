#!\bin\bash
while IFS='' read line; do
	filename_part=$(echo ${line} | cut -d '_' -f 1)
	filename="${filename_part}_doubleids.txt"
	grep -wf $line ../165test_passonly.noXnoY.hwe.nomiss.maf5percent_doubleids_cor.txt | awk '{print $1"\t"$2}' > "$filename"
done < filelist_eqtls.txt
