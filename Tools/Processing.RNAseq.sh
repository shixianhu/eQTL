#!/bin/bash

# $1=sample name

if [[ $# != 1 ]]
then
  break
  echo "No sample found!!!"
fi

ml FastQC
ml R
ml STAR
ml Java
ml SAMtools
ml sambamba

sample=$1
echo "===================> $sample now is processing"

out="/groups/umcg-weersma/tmp01/Shixian/RNA-seq/Output"
mkdir -p $out/fastqc_pre
mkdir -p $out/fastqc_post
mkdir -p $out/samtools
mkdir -p $out/star
mkdir -p $out/htseq
mkdir -p $out/trimmomatic

input="/groups/umcg-weersma/tmp01/Shixian/RNA-seq/test_samples/${sample}/"

# FasQC check
fastqc -q -t 4 $input/${sample}\_1.fq.gz $input/${sample}\_2.fq.gz --outdir $out/fastqc_pre
echo -e "#{sample} pre fastqc is done"

# Trimmomatic
java -jar /home/umcg-hushixian/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -phred33 /$input/${sample}\_1.fq.gz /$input/${sample}\_2.fq.gz \
   $out/trimmomatic/${sample}\_1_paired.fq $out/trimmomatic/${sample}\_1_unpaired.fq \
   $out/trimmomatic/${sample}\_2_paired.fq $out/trimmomatic/${sample}\_2_unpaired.fq \
   ILLUMINACLIP:/home/umcg-hushixian/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 HEADCROP:8 MINLEN:50
echo -e "${sample} Trimmomatic is done"

# FastQC again
fastqc -q -t 4 $out/trimmomatic/${sample}\_1_paired.fq $out/trimmomatic/${sample}\_2_paired.fq --outdir $out/fastqc_post
echo -e "${sample} post fasqc is done"

# STAR alignment
cd $out/star/
STAR --outFileNamePrefix ${sample}.star --readFilesIn $out/trimmomatic/${sample}\_1_paired.fq $out/trimmomatic/${sample}\_2_paired.fq \
 --genomeDir "/groups/umcg-weersma/tmp01/Shixian/RNA-seq/pipeline/database/StarIndex/"  --genomeLoad NoSharedMemory --runThreadN 6  --outFilterMultimapNmax 5 --outFilterMismatchNmax 8 --outReadsUnmapped Fastx --outSAMunmapped Within --outSAMstrandField intronMotif
echo -e "${sample} star is done"

# SAMtools and Picard 
samtools view -bS -o $out/samtools/${sample}.bam $out/star/${sample}.starAligned.out.sam
cd $out/samtools/
sambamba sort -o ${sample}.sorted.bam -p -t 50  ${sample}.bam
sambamba index -p -t 50 ${sample}.sorted.bam
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${sample}.bam O=${sample}.marked_duplicates.bam METRICS_FILE=${sample}.marked_dup_metrics.txt  ASSUME_SORT_ORDER=queryname
echo -e "${sample} mark duplicates is done"
java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics I=${sample}.bam O=${sample}.insert_size_metrics.txt M=0.5 H=${sample}.insert_size_histogram.pdf
echo -e "${sample} insert size check is done"

# HTseq gene count
cd $out/htseq
python /apps/software/HTSeq/0.9.1-foss-2015b-Python-2.7.11/bin/htseq-count -s no -m union -t exon \
       -i gene_id $out/star/${sample}.starAligned.out.sam \
        //groups/umcg-gastrocol/tmp04/Inflamed_eQTL/RNAseq_data/test/database/gencode.v19.annotation.patched_contigs.gtf > ${sample}.count.txt
echo -e "${sample} htseq is done"

# clean processed files
echo -e "${sample} bam files compressing"
zip -r $out/samtools/${sample}.BAM.files.zip $out/samtools/${sample}*bam
echo -e "${sample} sam files compressing"
zip -r $out/star/${sample}.SAM.files.zip $out/star/${sample}*sam
echo -e "${sample} trimmed fastq compressing"
zip -r $out/trimmomatic/${sample}*FQ.files.zip $out/trimmomatic/${sample}*fq

echo -e "===================> ${sample} processing is done"

