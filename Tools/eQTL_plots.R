setwd("Desktop")
vsd_test <- read.table(file="test_run_QTLs/complete_QC_counts_vsd_table_test.txt", header = T,check.names=F )
vsd <- read.table(file="Biopsies_docs/Biopsies_documents/complete_QC_counts_vsd_table_pil280.txt",header=T,check.names=F)
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

plot_list <- list()

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

  p <- ggplot(pheno_sim,aes(x=Genotype,y=Expression,color=Inflammation)) +
      geom_boxplot(color="black") +
      geom_jitter(position=position_jitter(0.2))+
      labs(title=paste("SNP",mut,"Gene",gene,"Pv",pval,"FDR",fdr),x="Genotype",y="Normalized Counts")

  plot_list[[i]] <- p
}

pdf(file="test_run_QTLs/eQTL_plots/165test_clumpstrict_maf5percent_modinf_nomiss_FDR5percent_clumped_top100_topsnps.pdf")
for (j in 1:100){
  print(plot_list[[j]])
}
dev.off()  


##Inflammation related dataset (use with relevant dataset)

#infvsnin_samples <- read.table("test_run_QTLs/infvsnin_RNA_samples.txt",header=T,stringsAsFactors=F)
#hosptowes <- read.table("WES_documents/UMCGtoWESconv.txt", sep="\t", stringsAsFactors=F, header = T)

#infsamples <- infvsnin_samples[infvsnin_samples$Inflammation=="I",]
#rownames(infsamples) <- infsamples$ID
#ninsamples <- infvsnin_samples[infvsnin_samples$Inflammation=="NI",]
#rownames(ninsamples) <- ninsamples$ID

#pheno_test_inf <- pheno_vsd2[rownames(infsamples),]
#vsd_test_inf <- vsd[,rownames(pheno_test_inf)]
#vsd_test_inf <- cbind(vsd$ID,vsd_test_inf)
#colnames(vsd_test_inf)[1] <- "ID"
#selected_patients_inf <- umcg_to_biopsy[match(rownames(pheno_test_inf),row.names(umcg_to_biopsy)),]
#temp_vec_inf <- NULL 
#for(x in selected_patients_inf){
#  temp_vec_inf <- c(temp_vec_inf,grep(unlist(strsplit(x,"_"))[1],hosptowes$UMCG_ID))
#}
#hosptowes2 <- hosptowes[temp_vec_inf,]
#selection <- genotypes_test[,as.character(hosptowes2$UMCG_ID)]
#selection <- cbind(genotypes_test$ID,selection)
#colnames(selection)[1] <- "ID"

#pheno_test_nin <- pheno_vsd2[rownames(ninsamples),]
#selected_patients_nin <- umcg_to_biopsy[match(rownames(pheno_test_nin),row.names(umcg_to_biopsy)),]
#vsd_test_nin <- vsd[,rownames(pheno_test_nin)]
#vsd_test_nin <- cbind(vsd$ID,vsd_test_nin)
#colnames(vsd_test_nin)[1] <- "ID"

#nin75list <- read.table(file="test_run_QTLs/nin75_gene_snps_list.txt", header=T, stringsAsFactors=F )
#nin75list_clean <- nin75list[which(startsWith(nin75list$ID,"rs")),]
#genelist <- nin75list_clean$Ensid

#nin75_eQTLs <- read.table(file="test_run_QTLs/liu_250KB_nin75_FDR_noHLA_ciseQTLs.txt", header=T, stringsAsFactors=F )
#nin75only_eQTLs <- nin75_eQTLs[match(nin75list_clean$ID,nin75_eQTLs$snps),]

#pdf(file="test_run_QTLs/eQTL_plots/liu_test_nin75_infvsnin.pdf")

#for (j in 1:length(nin75list_clean$ID))
#for (j in 1:length(inf75list_clean$ID))  {
  #mutation <- nin75list_clean$ID[j]
  #mutation <- inf75list_clean$ID[j]
  #gene <- ensidtogene[ensidtogene$Gene_ID == genelist[j],2]
  #pvalue <- signif(nin75only_eQTLs$pvalue[j],digits=3)
  #FDR <- signif(nin75only_eQTLs$FDR[j],digits=3)
  #pvalue <- signif(inf75only_eQTLs$pvalue[j],digits=3)
  #FDR <- signif(inf75only_eQTLs$FDR[j],digits=3)
  
  
  #pheno_sim_inf <- pheno_test_inf[,1:4]
  #pheno_sim_inf$Location[grep("ileum",pheno_sim_inf$Location,invert=T)] <- "colon"
  #pheno_sim_nin <- pheno_test_nin[,1:4]
  #pheno_sim_nin$Location[grep("ileum",pheno_sim_nin$Location,invert=T)] <- "colon"

  #pheno_sim_inf <- cbind(pheno_sim_inf,unlist(selection[which(selection$ID == mutation),-1]),
  #                     unlist(vsd_test_inf[which(vsd_test_inf$ID == genelist[j]),-1]))
  #colnames(pheno_sim_inf)[5:6] <- c("Genotype","Expression")
  #pheno_sim_nin <- cbind(pheno_sim_nin,unlist(selection[which(selection$ID == mutation),-1]),
  #                     unlist(vsd_test_nin[which(vsd_test_nin$ID == genelist[j]),-1]))
  #colnames(pheno_sim_nin)[5:6] <- c("Genotype","Expression")

##Temporary solution for alleles

  #homref <- "AA"
  #het <- "AB"
  #homalt <- "BB"
#####

  #pheno_sim_inf$Genotype <- gsub("0",homref,pheno_sim_inf$Genotype)
  #pheno_sim_inf$Genotype <- gsub("1",het,pheno_sim_inf$Genotype)
  #pheno_sim_inf$Genotype <- gsub("2",homalt,pheno_sim_inf$Genotype)

  #pheno_sim_nin$Genotype <- gsub("0",homref,pheno_sim_nin$Genotype)
  #pheno_sim_nin$Genotype <- gsub("1",het,pheno_sim_nin$Genotype)
  #pheno_sim_nin$Genotype <- gsub("2",homalt,pheno_sim_nin$Genotype)

  #geno_labels <- c(homref,het,homalt)

  #boxplot(pheno_sim_inf$Expression[pheno_sim_inf$Genotype == homref],
  #      pheno_sim_inf$Expression[pheno_sim_inf$Genotype == het],
  #      pheno_sim_inf$Expression[pheno_sim_inf$Genotype == homalt],
  #      pheno_sim_nin$Expression[pheno_sim_nin$Genotype == homref],
  #      pheno_sim_nin$Expression[pheno_sim_nin$Genotype == het],
  #      pheno_sim_nin$Expression[pheno_sim_nin$Genotype == homalt],
  #      names=c(rep(geno_labels,2)),col=c(rep("red",3),rep("blue",3)),
  #      main=paste("SNP",mutation,"Gene",gene,"Pv",pvalue,"FDR",FDR))
  #stripchart(list(pheno_sim_inf$Expression[pheno_sim_inf$Genotype == homref],
  #              pheno_sim_inf$Expression[pheno_sim_inf$Genotype == het],
  #              pheno_sim_inf$Expression[pheno_sim_inf$Genotype == homalt],
  #              pheno_sim_nin$Expression[pheno_sim_nin$Genotype == homref],
  #              pheno_sim_nin$Expression[pheno_sim_nin$Genotype == het],
  #              pheno_sim_nin$Expression[pheno_sim_nin$Genotype == homalt]),
  #         vertical=T,method="jitter", add=T, pch=20)
#}

#dev.off()
