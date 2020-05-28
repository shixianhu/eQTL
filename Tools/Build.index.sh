#!/bin/bash
#SBATCH --job-name=GenomeIndex
#SBATCH --error=GenomeIndex.err
#SBATCH --output=GenomeIndex.out
#SBATCH --mem=60gb
#SBATCH --time=15:59:00
#SBATCH --cpus-per-task=6

ml STAR

STAR --runThreadN 5 \
--runMode genomeGenerate \
--genomeDir /groups/umcg-weersma/tmp01/Shixian/RNA-seq/pipeline/database/StarIndex \
--genomeFastaFiles /groups/umcg-weersma/tmp01/Shixian/RNA-seq/pipeline/database/Homo_sapiens_assembly19.fasta \
--sjdbGTFfile /groups/umcg-weersma/tmp01/Shixian/RNA-seq/pipeline/database/gencode.v19.annotation.patched_contigs.gtf \
--sjdbOverhang 150
