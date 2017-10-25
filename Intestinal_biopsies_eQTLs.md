Intestinal Biopsies eQTL
=========================

Creator: Ruggero Barbieri

Year: 2017

0.Filter the VCF file
--------------------------------------------
#Locally or on the cluster do:
```
module load VCFtools; # cluster only
vcftools --vcf input.vcf --keep wes_165_ids.txt --recode --out input_165; #extracting the 165 samples
vcftools --vcf input_165.vcf --recode --remove-filtered-all --out input_165_passonly; #keeping only the variants that pass the QC filters
vcftools --vcf input_165_passonly.vcf --recode --not-chr X --not-chr Y --out input_165_passonly_noXnoY; #removing X and Y chromosomes
vcftools --vcf input_165_passonly_noXnoY.vcf --recode --hwe 0.000001 --out input_165_passonly_noXnoY_hwe; #filter by HW equilibrium
vcftools --vcf input_165_passonly_noXnoY_hwe.vcf --recode --maf 0.05 --out input_165_passonly_noXnoY_hwe_maf5percent; #filter by MAF
vcftools --vcf input_165_passonly_noXnoY_hwe_maf3percent.vcf --recode --max-missing 1 --out input_165_passonly_noXnoY_hwe_maf3percent_nomiss; #remove variants with missing genotypes
#Use the command "scp" to copy the resulting VCF to your machine if necessary
```

1.Prepare the data necessary for Matrix-eQTL
--------------------------------------------

#Files necessary: input VCF file, gene expression file, Identifier conversion files, Phenotypes files
#NOTE: The '#' on the sample names line in the VCF should be removed for the script to work
Execute meQTL_data_preparation.R; if on the cluster do:
```
module load R
Rscript meQTL_data_preparation

```
#Add 'ID' to the header for the first column

2.Run Matrix-eQTL and extract variants to clump
-----------------------------------------------

#Files necessary: Genotype file, Expression file, Gene positions file, SNP position files, Covariants files
Execute MatrixeQTL_pre-clumping_analysis.R; if on the cluster do:
```
module load R
Rscript MatrixeQTL_pre-clumping_analysis.R

```
#The "multi" file needs to be uploaded on the cluster for best results. Use the "scp" command to do it

3.Convert IDs for clumping
--------------------------
#Locally or on the cluster do:
```
grep -v "#" input.vcf | awk '{print $3"\t"$1":"$2":""$4":"$5}' > input_doubleids.txt
#Add the header to the resulting file: "rsid" for the first column, "customid" for the second
#In R do:
temp_doubleids <- read.table(file="input_doubleids.txt",header=T,stringsAsFactors=F)

temp_noids <- temp_doubleids[grep("\\.",temp_doubleids$rsid),]
temp_noids$rsid <- gsub(":","_",temp_noids$customid)
temp_rsids <- temp_doubleids[grep("\\.",temp_doubleids$rsid,invert = T),]
temp_doubleids_cor <- rbind(temp_rsids,temp_noids)

dup_check <- duplicated(temp_doubleids_cor$rsid)
temp_doubleids_cor$rsid[which(dup_check== TRUE)] <- paste(temp_doubleids_cor$rsid[which(dup_check== TRUE)],"2",sep="_")
dup_check <- duplicated(temp_doubleids_cor$rsid)
temp_doubleids_cor$rsid[which(dup_check== TRUE)] <- paste(temp_doubleids_cor$rsid[which(dup_check== TRUE)],"3",sep="_")
dup_check <- duplicated(temp_doubleids_cor$rsid)
temp_doubleids_cor$rsid[which(dup_check== TRUE)] <- paste(temp_doubleids_cor$rsid[which(dup_check== TRUE)],"4",sep="_")
dup_check <- duplicated(temp_doubleids_cor$rsid)
temp_doubleids_cor$rsid[which(dup_check== TRUE)] <- paste(temp_doubleids_cor$rsid[which(dup_check== TRUE)],"5",sep="_")
dup_check <- duplicated(temp_doubleids_cor$rsid)
temp_doubleids_cor$rsid[which(dup_check== TRUE)] <- paste(temp_doubleids_cor$rsid[which(dup_check== TRUE)],"6",sep="_")

write.table(temp_doubleids_cor, file="input_doubleids_cor.txt", quote=F, row.names=F, sep="\t")

#Use the "scp" command to upload this file to the cluster, in the same folder as the "multi" eQTL file  
```
4.Per gene clumping
-------------------
###Work in Progress

Plink command to execute:
```
#Do this on the cluster
module load Plink
plink --bfile /groups/umcg-weersma/tmp04/Ruggero/SN0090243/clumping_test/clumping_wes_test/westot --clump gene_file.txt --clump-snp-field customid --clump-field P --r2 dprime with-freqs --clump-verbose --clump-p1 1 --clump-p2 1 --clump-kb 500 --clump-r2 0.2 --out gene_clumped_file.clumped
```

5.Plot eQTLs boxplots
---------------------
#Locally, execute eQTL_plots.R (coming soon)
