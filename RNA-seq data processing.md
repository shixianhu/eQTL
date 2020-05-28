
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
```
sbatch Generate.process.sh $SAMPLE_PATH/
```


*Create job file for each sample*
---
- sbatch all jobs.
