#!\bin\bash
while IFS='' read line; do
	filename_part=$(echo ${line} | cut -d '_' -f 1)
	filename="${filename_part}_doubleids_list.txt"
	awk '{print $1}' $line > "$filename"
done < filelist_doubleids.txt
