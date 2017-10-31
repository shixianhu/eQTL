#!\bin\bash
while IFS='' read line; do
	filename_part=$(echo ${line} | cut -d '_' -f 1)
	filename1="${filename_part}_doubleids_complete.txt"
	filename2="${filename_part}_doubleids.txt"
	grep -w $filename_part ../165test_clumpstrict_hwe_nomiss_maf5percent_modinf_FDR5percent_multi_ciseQTLs.txt > temp_qtls.txt
	while IFS='' read variant; do 
		grep -w $variant temp_qtls.txt | awk '{print $4}' >> temp_file.txt 
	done < $line
	cut -f 1 temp_file.txt | paste "$filename2" - > "$filename1"
	rm temp_file.txt
done < filelist_doubleids_list.txt 
