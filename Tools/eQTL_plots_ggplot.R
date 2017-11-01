setwd("Desktop")

ploteqtl <- function(x) {
  ggplot(x,aes_string(x="Genotype",y="Expression",color="Inflammation")) +
    geom_boxplot(color="black") +
    geom_jitter(position=position_jitter(0.2))+
    labs(title=paste("SNP",mut,"Gene",gene,"Pv",pval,"FDR",fdr),x="Genotype",y="Normalized Counts")
}

ploteqtlcompare <- function(x) {
  ggplot() +
    geom_boxplot(data=x,aes(x=Genotype,y=Expression,factor(Inflammation),color=Inflammation)) +
    geom_jitter(position=position_jitter(width=.1, height=0))+
    labs(title=paste("SNP",mut,"Gene",gene,"Pv",pval,"FDR",fdr),x="Genotype",y="Normalized Counts")
}

vsd_test <- read.table(file="test_run_QTLs/complete_QC_counts_vsd_table_test.txt", header = T,check.names=F )
vsd_test <- read.table(file="test_run_QTLs/complete_QC_counts_vsd_table_test.txt", header = T,check.names=F )
umcg_to_biopsy <- read.table("Biopsies_docs/Biopsies_documents/UMCGtoBiopsy.txt", header=T, stringsAsFactors = F)
row.names(umcg_to_biopsy) <- umcg_to_biopsy$Biopsy_ID
umcg_to_biopsy$Biopsy_ID <- NULL
genotypes_test_we <- read.table("test_run_QTLs/MeQTL_165test_clumpstrict_maf5percent_hwe_nomiss_genotype.txt",header=T,stringsAsFactors = F,fill=T)
hosptowes <- read.table("WES_documents/UMCGtoWESconv.txt", sep="\t", stringsAsFactors=F, header = T)

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

pheno_tot <- rbind(pheno_p,pheno1)
pheno_tot$Date_category <- NULL
colnames(pheno_100) <- colnames(pheno_tot)
pheno_all <- rbind(pheno_tot,pheno_100)
rownames(pheno_100_2)[1] <- "807_6Atris"
colnames(pheno_100_2) <- colnames(pheno_all)
pheno_all2 <- rbind(pheno_all,pheno_100_2)

pheno_vsd <-pheno_all2[colnames(vsd),]
pheno_vsd2 <-pheno_vsd[-1,]
pheno_vsd2$Inflammation <- gsub("light_I","NI",pheno_vsd2$Inflammation)

pheno_test <- pheno_all2[colnames(vsd_test),]
pheno_test2 <- pheno_test[-1,]
pheno_test2$Inflammation <- gsub("light_I","NI",pheno_test2$Inflammation)

colnames(vsd_test) <- colnames(genotypes_test_we)
rownames(pheno_test2) <- colnames(genotypes_test_we)[-1]

ensidtogene <- read.table("Ensembl_GRCh37_ID_comparison.txt",header=T,stringsAsFactors=F, fill=T)

eQTLs_table <- read.table("test_run_QTLs/165test_clumpstrict_hwe_nomiss_maf5percent_modinf_FDR5percent_all_clumped_ciseQTLs.txt",header=T,stringsAsFactors=F)
eQTLs_table <- eQTLs_table[order(as.numeric(eQTLs_table$FDR)),]
mutlist <- eQTLs_table$snps
mutlist_clean <- mutlist[startsWith(mutlist,"rs")]
mutlist_clean_pass <- intersect(mutlist_clean,genotypes_test_we$ID)
eQTLs_table_pass <- eQTLs_table[eQTLs_table$snps %in% mutlist_clean_pass,]
genelist <- eQTLs_table_pass$gene


for (i in 1:100) {
  
  mut <- eQTLs_table_pass$snps[i]
  gene <- ensidtogene[ensidtogene$Gene_ID == genelist[i],2]
  pval <- signif(eQTLs_table_pass$pvalue[i],digits=2)
  fdr <- signif(eQTLs_table$FDR[i],digits=2)
  
  pheno_sim <- pheno_test2[,1:4]
  pheno_sim$Location[grep("ileum",pheno_sim$Location,invert=T)] <- "colon"
  pheno_sim <- cbind(pheno_sim,unlist(genotypes_test_we[which(genotypes_test_we$ID == mut),-1]),
                     unlist(vsd_test[which(vsd_test$ID == genelist[i]),-1]))
  colnames(pheno_sim)[5:6] <- c("Genotype","Expression")
  
  homref <- "AA"
  het <- "AB"
  homalt <- "BB"
  
  pheno_sim$Genotype <- gsub("0",homref,pheno_sim$Genotype)
  pheno_sim$Genotype <- gsub("1",het,pheno_sim$Genotype)
  pheno_sim$Genotype <- gsub("2",homalt,pheno_sim$Genotype)
  geno_labels <- c(homref,het,homalt)
  
  pheno_sim_inf <- pheno_sim[which(pheno_sim$Inflammation == "I"),]
  pheno_sim_nin <- pheno_sim[which(pheno_sim$Inflammation == "NI"),]
  
  eqtl_plot<- ploteqtl(pheno_sim)
  eqtl_plot_compare<- ploteqtlcompare(pheno_sim)
  
  #pdf(eqtl_plot,file=paste0("test_run_QTLs/eQTL_plots/eQTL_plots_165test_maf5percent_modinf/",gene,"-",mut,".pdf")) #crea tabella per tutti i pdf
  #print(eqtl_plot)
  pdf(eqtl_plot,file=paste0("test_run_QTLs/eQTL_plots/eQTL_plots_165test_maf5percent_modinf/",gene,"-",mut,"_inf_compare.pdf"))
  print(eqtl_plot_compare)
  dev.off()
  pheno_sim <- NULL
}
