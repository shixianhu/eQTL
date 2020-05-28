
# RNAseq raw data processing

*Samples input folder check*
---
- Create a folder contains all samples
- One sample has one folder with two fastq files inside. sample_1.fastq and sample_2.fastq


*Build reference index for STAR mapping*
---
- Check STAR genome index.
- Reference genome is from GTEx V7 (https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq)

```
sbatch Build.index.sh
```


*Create job file for each sample*
---
- the pipeline in *Generate.process.sh* contains the following:

1. fasqc for pre-processed raw sequencing reads
2. remove adapters, primers (illumina true-seq) and low-quality reads (average score 25, minlen 50)
3. fasqc for post-processed sequecing reads
4. star alingemnt to reference genome 
5. samtools and picard to generate marked duplicates and insert size metrix
6. htseq to genenrate gene-level counts
7. compress and clean processed sam and bam files

```
sbatch Generate.process.sh $SAMPLE_PATH/
```


*Submit all jobs*
---
- sbatch all jobs.
