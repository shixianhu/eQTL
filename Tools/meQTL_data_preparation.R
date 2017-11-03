###Data preparation for a Matrix-eQTL run

#Creator: Ruggero Barbieri

#Date:2016

# Read in counts expression and ID conversion
setwd("Desktop")
vsd_test <- read.table(file="test_run_QTLs/complete_QC_counts_vsd_table_test.txt", header = T,check.names=F )
umcg_to_biopsy <- read.table("Biopsies_docs/Biopsies_documents/UMCGtoBiopsy.txt", header=T, stringsAsFactors = F)
row.names(umcg_to_biopsy) <- umcg_to_biopsy$Biopsy_ID
umcg_to_biopsy$Biopsy_ID <- NULL
# Select right pts from counts with ID conversion
selected_patients <- umcg_to_biopsy[match(colnames(vsd_test)[-1],row.names(umcg_to_biopsy)),]

# Read in VCF and link to UMCG research ID
input_vcf <- read.table("test_run_QTLs/Liu_hwe_regions_250KB_unique.165test.passonly.noXnoY.hwe.maf3percent.recode.vcf", sep="\t", stringsAsFactors=F, header = T, check.names = F)
hosptowes <- read.table("WES_documents/UMCGtoWESconv.txt", sep="\t", stringsAsFactors=F, header = T)

temp_vec <- NULL 
for(x in selected_patients){
  temp_vec <- c(temp_vec,grep(unlist(strsplit(x,"_"))[1],hosptowes$UMCG_ID))
}

hosptowes2 <- hosptowes[temp_vec,]
selection <- input_vcf[,as.character(hosptowes2$WES_ID)]
colnames(selection) <- hosptowes2$UMCG_ID

# Create a custom ID for variants without RS ID: chr_pos_refA_altA
temp_id_list <- vector(length = length(rownames(input_vcf)))
for(y in 1:length(rownames(input_vcf))){
  if(liu_vcf$ID[y] == ".") {
          temp_id_list[y] <- paste(input_vcf$CHROM[y],input_vcf$POS[y],input_vcf$REF[y],input_vcf$ALT[y],sep="_")}
  else {temp_id_list[y] <- input_vcf$ID[y]}
}

# Make variant names unique (indels vary in length)
temp_id_list2 <- temp_id_list
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"2",sep="_")
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"3",sep="_")
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"4",sep="_")
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"5",sep="_")
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"6",sep="_")

# Creat new table with genotypes
sel_sim <- matrix(ncol=dim(selection)[2],nrow=dim(selection)[1])
for (z in 1:dim(selection)[1]) {
  a <- selection[z,]
  a[grepl("^0/0",a)] <- 0
  a[grepl("^0/1",a)] <- 1
  a[grepl("^1/1",a)] <- 2
  sel_sim[z,] <- as.numeric(a)
}

#sel_sim[is.na(sel_sim)] <- 0 # deprecated
sel_sim <- cbind(temp_id_list2,sel_sim)
rownames(sel_sim) <- sel_sim[,1]
sel_sim <- sel_sim[,-1]
colnames(sel_sim) <- colnames(selection)

# create SNP position table with identifiers and chromosomes in position
snp_pos <- cbind(temp_id_list2,input_vcf$CHROM,input_vcf$POS)
colnames(snp_pos) <- c("snp","chr","pos")

# write genotype and SNP position files
write.table(sel_sim,"test_run_QTLs/MeQTL_Liu_165test_clumpstrict_maf3percent_genotype.txt",sep="\t", quote=F)
write.table(snp_pos,"test_run_QTLs/MeQTL_Liu_165test_clumpstrict_maf3percent_snp_pos.txt",sep="\t", quote=F,row.names = F)

##Preparing covariates files

# loading and pasting phenotype files, per batch
pheno_p <- read.csv("Biopsies_docs/Biopsies_documents/Pilot_proposal.csv", header=T, sep=";", row.names = 2)
pheno_p$Research_ID <- NULL
pheno_p$Sequence_run <- c(rep("Pilot", 20))

values <- read.table("Biopsies_docs/Biopsies_documents/1603_IBD_RNAseq_80.expression.genelevel.v75.htseq.txt.table.txt", header=T, check.names=F, row.names = 1 )
pheno <- read.csv("Biopsies_docs/Biopsies_documents/Plate_80.csv", header=T, sep=";",row.names = 1)
pheno1 <- pheno[colnames(values),]
pheno1$Sequence_run <- c(rep("Second", 79))
colnames(pheno1) <- colnames(pheno_p)

pheno_100 <- read.csv("Biopsies_docs/Biopsies_documents/Plate_100.csv", header=T, sep=";",row.names = 1)
pheno_100$Sequence_run <- c(rep("Third", 100))

pheno_100_2 <- read.csv("Biopsies_docs/Biopsies_documents/Plate_100_2.csv", header=T, sep=";",row.names = 1)
pheno_100_2$Sequence_run <- c(rep("Fourth", 100))

# paste all phenofiles
pheno_tot <- rbind(pheno_p,pheno1)
pheno_tot$Date_category <- NULL
colnames(pheno_100) <- colnames(pheno_tot)
pheno_all <- rbind(pheno_tot,pheno_100)
rownames(pheno_100_2)[1] <- "807_6Atris"
colnames(pheno_100_2) <- colnames(pheno_all)
pheno_all2 <- rbind(pheno_all,pheno_100_2)

# load extra files for sex and age of biopsy and add to big pheno table
sex_280 <- read.table("Biopsies_docs/Biopsies_documents/UMCGtoBiopsy+Sex.txt",header=T, row.names = 2,check.names = F)
ageatbiopsy_280 <- read.table("Biopsies_docs/Biopsies_documents/AgeAtBiopsy.txt",header=T, row.names = 1,check.names = F)

# paste sex and age to pheno table
pheno_all2 <- cbind(pheno_all2,sex_280[rownames(pheno_all2),])
pheno_all2 <- cbind(pheno_all2,ageatbiopsy_280[rownames(pheno_all2),])
colnames(pheno_all2)[11] <- "age_at_biopsy"

# load the counts file again
vsd_test <- read.table(file="test_run_QTLs/complete_QC_counts_vsd_table_test.txt", header = T,check.names=F )
# select from pheno file the patients of interest, using identifiers form the counts table (for example the 165-pt set)
pheno_test <- pheno_all2[colnames(vsd_test),]
pheno_test2 <- pheno_test[-1,]

# create numeric variable for 'batch'
batch <- pheno_test2$Sequence_run
batch[grepl("Pilot",batch)] <- 1
batch[grepl("Second",batch)] <- 2
batch[grepl("Third",batch)] <- 3
batch[grepl("Fourth",batch)] <- 4
# create numeric variable for 'sex'
sex <- as.character(pheno_test2$Sex)
sex[grepl("Male",sex)] <- 0
sex[grepl("Female",sex)] <- 1
# load 'age'
age <- pheno_test2$age_at_biopsy

# conversion to numeric of inplammation, diagnosis and tissue
inf <- as.character(pheno_test2$Inflammation)
inf[grepl("light_I",inf)] <- "NI"
inf[grepl("NI",inf)] <- 0
inf[grepl("I",inf)] <- 1

diag <- as.character(pheno_test2$Diagnosis)
diag[grepl("IBDU",diag)] <- "Ulcerative colitis"
diag[grepl("Crohn's disease",diag)] <- 0
diag[grepl("Ulcerative colitis",diag)] <- 1

loc <- as.character(pheno_test2$Location)
loc[grep("ileum",loc,invert=T)] <- "colon"
loc[grepl("colon",loc)] <- 0
loc[grepl("ileum",loc)] <- 1

# insert column 'ID'
batch_id <- cbind(colnames(sel_sim),batch)
colnames(batch_id)[1] <- "ID"
batch_id2 <- t(batch_id) # transpose

cov_id <- cbind(colnames(sel_sim),batch,sex,age)
colnames(cov_id)[1] <- "ID"
cov_id2 <- t(cov_id) # transpose

labeled_allcov <- t(cbind(colnames(sel_sim),batch,sex,age,inf,diag,loc)) # transpose

write.table(batch_id2,"MeQTL_gavin_test_cov_batch.txt",sep="\t", quote=F, col.names = F)
write.table(cov_id2,"MeQTL_liu_test_cov.txt",sep="\t", quote=F, col.names = F)
write.table(labeled_allcov,"MeQTL_liu_test_allcov.txt",sep="\t", quote=F, col.names = F)

##Preparing Inflammation based data
vsd_test_inf <- read.table(file="test_run_QTLs/infOnly_QC_dedup_values_R.txt", header = T,check.names=F )
umcg_to_biopsy <- read.table("Biopsies_docs/Biopsies_documents/UMCGtoBiopsy.txt", header=T, stringsAsFactors = F)
row.names(umcg_to_biopsy) <- umcg_to_biopsy$Biopsy_ID
umcg_to_biopsy$Biopsy_ID <- NULL
selected_patients <- umcg_to_biopsy[match(colnames(vsd_test_inf),row.names(umcg_to_biopsy)),]
input_vcf <- read.table("test_run_QTLs/Liu_hwe_regions_250KB.recode.unique.vcf", sep="\t", stringsAsFactors=F, header = T, check.names = F)
hosptowes <- read.table("WES_documents/UMCGtoWESconv.txt", sep="\t", stringsAsFactors=F, header = T)

temp_vec <- NULL 
for(x in selected_patients){
  temp_vec <- c(temp_vec,grep(unlist(strsplit(x,"_"))[1],hosptowes$UMCG_ID))
}

hosptowes2 <- hosptowes[temp_vec,]

temp_id_list <- vector(length = length(rownames(input_vcf)))
for(y in 1:length(rownames(input_vcf))){
  if(input_vcf$ID[y] == ".") {
    temp_id_list[y] <- paste(input_vcf$CHROM[y],input_vcf$POS[y],input_vcf$REF[y],input_vcf$ALT[y],sep="_")}
  else {temp_id_list[y] <- input_vcf$ID[y]}
}

temp_id_list2 <- temp_id_list
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"2",sep="_")
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"3",sep="_")
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"4",sep="_")
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"5",sep="_")
dup_check <- duplicated(temp_id_list2)
temp_id_list2[which(dup_check== TRUE)] <- paste(temp_id_list2[which(dup_check== TRUE)],"6",sep="_")

colnames(vsd_test_inf) <- hosptowes2$UMCG_ID
selection <- input_vcf[,as.character(hosptowes2$WES_ID)]
colnames(selection) <- hosptowes2$UMCG_ID

sel_sim <- matrix(ncol=dim(selection)[2],nrow=dim(selection)[1])
for (z in 1:dim(selection)[1]) {
  a <- selection[z,]
  a[grepl("^0/0",a)] <- 0
  a[grepl("^0/1",a)] <- 1
  a[grepl("^1/1",a)] <- 2
  #print(a)
  sel_sim[z,] <- as.numeric(a)
}
sel_sim[is.na(sel_sim)] <- 0
sel_sim <- cbind(temp_id_list2,sel_sim)
rownames(sel_sim) <- sel_sim[,1]
sel_sim <- sel_sim[,-1]
colnames(sel_sim) <- colnames(selection)
write.table(vsd_test_inf,"MeQTL_inf_test_counts.txt", sep="\t", quote=F)
write.table(sel_sim,"test_run_QTLs/MeQTL_liu_inf_250KB_genotype.txt",sep="\t", quote=F)

pheno_inf <- read.table(file="test_run_QTLs/infOnly_QC_dedup_pheno_R.txt", header = T,check.names=F,stringsAsFactors=F )

batch <- pheno_inf$Sequence_run
batch[grepl("Pilot",batch)] <- 1
batch[grepl("Second",batch)] <- 2
batch[grepl("Third",batch)] <- 3
batch[grepl("Fourth",batch)] <- 4

sex <- as.character(pheno_inf$Sex)
sex[grepl("Male",sex)] <- 0
sex[grepl("Female",sex)] <- 1

age <- pheno_inf$age_at_biopsy

diag <- as.character(pheno_inf$Diagnosis)
diag[grepl("IBDU",diag)] <- "Ulcerative colitis"
diag[grepl("Crohn's disease",diag)] <- 0
diag[grepl("Ulcerative colitis",diag)] <- 1

loc <- as.character(pheno_inf$Location)
loc[grep("ileum",loc,invert=T)] <- "colon"
loc[grepl("colon",loc)] <- 0
loc[grepl("ileum",loc)] <- 1


cov_id <- cbind(colnames(sel_sim),batch,sex,age)
colnames(cov_id)[1] <- "ID"
cov_id2 <- t(cov_id)

labeled_allcov_inf <- t(cbind(colnames(sel_sim),batch,sex,age,inf,diag,loc))

write.table(cov_id2,"MeQTL_liu_inf_cov.txt",sep="\t", quote=F, col.names = F)
write.table(labeled_allcov_inf,"MeQTL_liu_inf_allcov.txt",sep="\t", quote=F, col.names = F)
