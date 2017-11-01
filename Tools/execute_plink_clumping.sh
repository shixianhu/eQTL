#!/bin/bash
#SBATCH --job-name=clumpingpergene
#SBATCH --output=clumpingpergene_165test_maf5_nomiss_allcor.out
#SBATCH --error=clumpingpergene_165test_maf5_nomiss_allcor.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --nodes=1

module load plink
while IFS='' read line; do
	filename_part=$(echo ${line} | cut -d '_' -f 1)
	filename="${filename_part}_clumpedstrict"
	plink --bfile /groups/umcg-weersma/tmp04/Ruggero/SN0090243/clumping_test/clumping_wes_test/westot --clump "${line}" --clump-snp-field customID --clump-field P --r2 dprime with-freqs --clump-verbose --clump-p1 1 --clump-p2 1 --clump-kb 500 --clump-r2 0.2 --out "${filename}"
done < filelist_doubleids_complete.txt
