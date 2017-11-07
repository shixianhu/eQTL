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
vcftools --vcf input_165_passonly_noXnoY.vcf --recode --hwe 0.001 --out input_165_passonly_noXnoY_hwe; #filter by HW equilibrium
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
#Files necessary: input VCF file (eg: input.VCF)
#Locally or on the cluster do:
```
grep -v "#" input.vcf | awk '{print $3"\t"$1":"$2":"$4":"$5}' > input_doubleids.txt
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
#Files necessary: 'multi variant per gene' eQTL file (eg: input_multi_eQTLs.txt), ID conversion file (eg: input_doubleids_cor.txt) 

#Do this on the cluster:
```
awk '{print $2}' input_multi_ciseQTLs.txt | uniq | grep "ENS" > geneids_toclump.txt
bash get_variant_per_gene.sh
ls -lh | awk '{print $9}' | grep 'ENSG[0-9]*_eqtls\.txt' > filelist_eqtls.txt
bash get_customid_per_variant.sh
ls -lh | awk '{print $9}' | grep 'ENSG[0-9]*_doubleids_complete.txt' > filelist_doubleids_complete.txt
bash get_variant_from_doubleids.sh
ls -lh | awk '{print $9}' | grep "[0-9]_doubleids_list\.txt" > filelist_doubleids_list.txt
for i in *_doubleids_complete.txt; do sed -i '1i rsID\tcustomID\tP' $i; done
sbatch execute_plink_clumping.sh
for i in *_clumpedstrict.clumped; do sed -i 's/  */ /g' $i; done
grep "^ [0-9]* 1" *_clumpedstrict.clumped > summary_clumpedstrict.clumped.txt
awk 'BEGIN{FS=OFS=" "}{$1="";}1' summary_clumpedstrict.clumped.txt > input_summary.clumped.txt
```

5.Plot eQTLs boxplots
---------------------
#Files necessary: Phenotype files, clumped eQTLs file, Genotype file
#Locally, execute eQTL_plots_ggplot.R;if on the cluster do:
```
module load R
Rscript eQTL_plots_ggplot.R

```
