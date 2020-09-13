args<-commandArgs(T)
file = args[1]

library(TwoSampleMR)
library(MRInstruments)
library(coloc)

eQTL.hits=read.table(file,sep = "\t",header = T,stringsAsFactors = F)
eQTL.hits=eQTL.hits[,c("ExpressionGene","Chr","rsID","Pos","Allele1","Allele0","AllelFre","Beta","SE","p_wald")]

IBD.data <- extract_outcome_data(
  snps = eQTL.hits$rsID,
  outcomes = 'ebi-a-GCST004131',access_token=NULL
)
CD.data <- extract_outcome_data(
  snps = eQTL.hits$rsID,
  outcomes = 'ebi-a-GCST004132',access_token=NULL
)
UC.data <- extract_outcome_data(
  snps = eQTL.hits$rsID,
  outcomes = 'ebi-a-GCST004133',access_token=NULL
)
coeliac.data<- extract_outcome_data(
  snps = eQTL.hits$rsID,
  outcomes = 'ukb-b-8631',access_token=NULL
)
diverticulitis.data<- extract_outcome_data(
  snps = eQTL.hits$rsID,
  outcomes = 'ukb-b-14796',access_token=NULL
)
colon.cancer.data<- extract_outcome_data(
  snps = eQTL.hits$rsID,
  outcomes = 'ukb-b-20145',access_token=NULL
)

IBD.data=IBD.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome","samplesize.outcome")]
CD.data=CD.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome","samplesize.outcome")]
UC.data=UC.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome","samplesize.outcome")]
coeliac.data=coeliac.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome","samplesize.outcome")]
diverticulitis.data=diverticulitis.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome","samplesize.outcome")]
colon.cancer.data=colon.cancer.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome","samplesize.outcome")]

tmp=eQTL.hits[eQTL.hits$rsID %in% IBD.data$SNP,]
IBD.data=IBD.data[IBD.data$SNP %in% tmp$rsID,]
IBD.data=merge(IBD.data,tmp[,c("AllelFre","rsID"),drop=F],by.x = "SNP",by.y="rsID")
tmp=tmp[order(tmp$rsID),]
IBD.data=IBD.data[order(IBD.data$SNP),]

IBD.coloc <- coloc.abf(dataset1=list(pvalues=tmp$p_wald,N=280,type="quant"),
                    dataset2=list(pvalues=IBD.data$pval.outcome,N=IBD.data$samplesize.outcome[1],type="quant"),
                    MAF=tmp$AllelFre)

tmp=eQTL.hits[eQTL.hits$rsID %in% CD.data$SNP,]
CD.data=CD.data[CD.data$SNP %in% tmp$rsID,]
CD.data=merge(CD.data,tmp[,c("AllelFre","rsID"),drop=F],by.x = "SNP",by.y="rsID")
tmp=tmp[order(tmp$rsID),]
CD.data=CD.data[order(CD.data$SNP),]

CD.coloc <- coloc.abf(dataset1=list(pvalues=tmp$p_wald,N=280,type="quant"),
                       dataset2=list(pvalues=CD.data$pval.outcome,N=CD.data$samplesize.outcome[1],type="quant"),
                       MAF=tmp$AllelFre)

tmp=eQTL.hits[eQTL.hits$rsID %in% UC.data$SNP,]
UC.data=UC.data[UC.data$SNP %in% tmp$rsID,]
UC.data=merge(UC.data,tmp[,c("AllelFre","rsID"),drop=F],by.x = "SNP",by.y="rsID")
tmp=tmp[order(tmp$rsID),]
UC.data=UC.data[order(UC.data$SNP),]

UC.coloc <- coloc.abf(dataset1=list(pvalues=tmp$p_wald,N=280,type="quant"),
                       dataset2=list(pvalues=UC.data$pval.outcome,N=UC.data$samplesize.outcome[1],type="quant"),
                       MAF=tmp$AllelFre)

tmp=eQTL.hits[eQTL.hits$rsID %in% coeliac.data$SNP,]
coeliac.data=coeliac.data[coeliac.data$SNP %in% tmp$rsID,]
coeliac.data=merge(coeliac.data,tmp[,c("AllelFre","rsID"),drop=F],by.x = "SNP",by.y="rsID")
tmp=tmp[order(tmp$rsID),]
coeliac.data=coeliac.data[order(coeliac.data$SNP),]

coeliac.coloc <- coloc.abf(dataset1=list(pvalues=tmp$p_wald,N=280,type="quant"),
                       dataset2=list(pvalues=coeliac.data$pval.outcome,N=coeliac.data$samplesize.outcome[1],type="quant"),
                       MAF=tmp$AllelFre)

tmp=eQTL.hits[eQTL.hits$rsID %in% diverticulitis.data$SNP,]
diverticulitis.data=diverticulitis.data[diverticulitis.data$SNP %in% tmp$rsID,]
diverticulitis.data=merge(diverticulitis.data,tmp[,c("AllelFre","rsID"),drop=F],by.x = "SNP",by.y="rsID")
tmp=tmp[order(tmp$rsID),]
diverticulitis.data=diverticulitis.data[order(diverticulitis.data$SNP),]

diverticulitis.coloc <- coloc.abf(dataset1=list(pvalues=tmp$p_wald,N=280,type="quant"),
                       dataset2=list(pvalues=diverticulitis.data$pval.outcome,N=diverticulitis.data$samplesize.outcome[1],type="quant"),
                       MAF=tmp$AllelFre)

tmp=eQTL.hits[eQTL.hits$rsID %in% colon.cancer.data$SNP,]
colon.cancer.data=colon.cancer.data[colon.cancer.data$SNP %in% tmp$rsID,]
colon.cancer.data=merge(colon.cancer.data,tmp[,c("AllelFre","rsID"),drop=F],by.x = "SNP",by.y="rsID")
tmp=tmp[order(tmp$rsID),]
colon.cancer.data=colon.cancer.data[order(colon.cancer.data$SNP),]

colon.cancer.coloc <- coloc.abf(dataset1=list(pvalues=tmp$p_wald,N=280,type="quant"),
                       dataset2=list(pvalues=colon.cancer.data$pval.outcome,N=colon.cancer.data$samplesize.outcome[1],type="quant"),
                       MAF=tmp$AllelFre)



IBD.result=as.data.frame(t(IBD.coloc$summary))
IBD.result$disease="IBD"

CD.result=as.data.frame(t(CD.coloc$summary))
CD.result$disease="CD"

UC.result=as.data.frame(t(UC.coloc$summary))
UC.result$disease="UC"

coeliac.result=as.data.frame(t(coeliac.coloc$summary))
coeliac.result$disease="coeliac"

diverticulitis.result=as.data.frame(t(diverticulitis.coloc$summary))
diverticulitis.result$disease="diverticulitis"

colon.cancer.result=as.data.frame(t(colon.cancer.coloc$summary))
colon.cancer.result$disease="colon.cancer"

results=rbind(IBD.result,CD.result,UC.result,coeliac.result,diverticulitis.result,colon.cancer.result)
results$Probe=eQTL.hits$ExpressionGene[1]

write.table(results,paste0(results$Prob[1],".coloc.txt"),sep = "\t",row.names = F,quote = F)
