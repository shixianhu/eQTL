#!/bin/bash
#BATCH --job-name=cis
#SBATCH --error=cis.err
#SBATCH --output=cis.out
#SBATCH --mem=8gb
#SBATCH --time=64:00:00
#SBATCH --cpus-per-task=6

ml plink

#awk '{print $2}' All_pairs.txt | sort | uniq > Probe.txt

export  LD_LIBRARY_PATH=/home/umcg-hushixian/gemma/gcc-5.4.0-lib-3x53yv4v144c9xp0/lib

cat xaa.txt | while read line

do

grep -w $line ../All_pairs.txt | awk '{print $1}' > tmp.snp.txt
plink --bfile /groups/umcg-gastrocol/tmp04/Inflamed_eQTL/Process/WholeAnalysis/Genotype/plink/whole.sample --extract tmp.snp.txt --make-bed --out tmp.analysis
awk -v col=$line 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' ../Matched_table/Reordered.phenotype.txt > tmp.expression.txt
sed -i '1d' tmp.expression.txt
~/gemma-0.98-linux-static -bfile tmp.analysis -p tmp.expression.txt -km 1 -k ../Kinship/IBS.mibs -lmm 4 -o $line.outcome -miss 0.99
rm tmp*

done
