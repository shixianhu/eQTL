#!/bin/bash
# $1 indicates the path of raw samples. In the input folder, one sample has one independent folder with two pair-end fastq files. The folder name should be the sample name. the fastq file should be sample_1.fastq and sample_2.fastq 

if [[ $# != 1 ]]
then
  break
  echo "Please enter the correct sample input path"
fi

input=$1

for file in $input/*

do
sample=$(basename $file)
echo '#!/bin/bash'  > rnaseq.${sample}.sh
echo "#SBATCH --job-name=RNAseq.${sample}" >> rnaseq.${sample}.sh
echo "#SBATCH --error=RNAseq.${sample}.err" >> rnaseq.${sample}.sh
echo "#SBATCH --output=RNAseq.${sample}.out" >> rnaseq.${sample}.sh
echo "#SBATCH --mem=15gb" >> rnaseq.${sample}.sh
echo "#SBATCH --time=6:00:00" >> rnaseq.${sample}.sh
echo "#SBATCH --cpus-per-task=6" >> rnaseq.${sample}.sh

echo "bash /groups/umcg-weersma/tmp01/Shixian/RNA-seq/pipeline/Processing.RNAseq.sh ${sample}" >> rnaseq.${sample}.sh

done
