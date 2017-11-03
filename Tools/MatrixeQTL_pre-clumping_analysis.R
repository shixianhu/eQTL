###Matrix-eQTL sample code with and added portion to extract non-independent eQTLs to be clumped

#Creator: Ruggero Barbieri

#Date:2017


# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)
base.dir = getwd();

## Location of the package with the data files.
#base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
#useModel = modelLINEAR_CROSS; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR;

# Genotype file name
snps_location_file_name = paste(base.dir,"test_whole_exome/MeQTL_165test_clumpstrict_maf5percent_hwe_nomiss_snp_pos.txt", sep="");
SNP_file_name = paste(base.dir,"MeQTL_165test_clumpstrict_maf5percent_hwe_nomiss_genotype.txt", sep="");

# Gene expression file name
expression_file_name = paste(base.dir,"MeQTL_gavin_test_counts.txt", sep="");
gene_location_file_name = paste(base.dir, "MeQTL_ready_gene_pos.txt", sep="");

# Covariates file name
# Set to character() for no covariates
#covariates_file_name =  character();
covariates_file_name = paste(base.dir,"MeQTL_liu_test_modinf_allcov_cov.txt", sep="");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Distance for local gene-SNP pairs, left and right from the border of the gene. since we're not looking at trans only 5e5 is fine
cisDist = 5e5;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 5000;      # read file in slices of 5,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  #pvOutputThreshold     = pvOutputThreshold_tra,
  pvOutputThreshold     = 0,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

# collect cis eQTLs in 'all'
all <- me$cis$eqtls
all.fdr <- all[which(all$FDR<0.05),] #subset snp-gene pairs with fdr less than 0.05

###Extract variant per gene from eQTLs
geneid_dup <- unique(as.character(all.fdr$gene[duplicated(all.fdr$gene)]))
all.multi <- NULL
for (x in geneid_dup){
  geneid <- x
  eqtls <- all.fdr[grep(x,all.fdr$gene),]
  all.multi <- rbind(all.multi,eqtls)
}

geneid_single <- setdiff(all.fdr$gene,geneid_dup)
all.single <- all.fdr[all.fdr$gene %in% geneid_single,]

write.table(all.fdr,file="165test_clumpstrict_hwe_nomiss_maf5percent_modinf_FDR5percent_ciseQTLs.txt", quote=F, row.names=F, sep="\t")
write.table(all.single,file="165test_clumpstrict_hwe_nomiss_maf5percent_modinf_FDR5percent_single_ciseQTLs.txt", quote=F, row.names=F, sep="\t")
write.table(all.multi,file="165test_clumpstrict_hwe_nomiss_maf5percent_modinf_FDR5percent_multi_ciseQTLs.txt", quote=F, row.names=F, sep="\t")
