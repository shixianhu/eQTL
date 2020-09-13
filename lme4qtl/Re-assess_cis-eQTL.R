# we used GEMMA for an interaction model (gene expression = SNP + inflammation + SNP * inflammation) to calculate inflammation-dependent cis-eQTL
# however, we only got summary statistics of interaction term from GEMMA by default
# to get full summary statistics of term 'SNP' and 'inflammation', we used lme4qtl R package to re-calculate 190 significant (FDR interaction <0.05) from GEMMA
# in addition, we also compared the results derived from GEMMA and lme4qtl, and compared three models to adjust for random effect

library(ggsci)
library(GMMAT)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(grid)
library(lme4)
library(lme4qtl)
library(foreach)
library(crayon)
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

# import gene expression table
expression=read.table("exlInflammationPCs_corrected.txt",
                      sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
expression=t(expression)
expression=as.data.frame(expression)

# import phenotype file
metadata=read.table("complete_pheno_table_pil280.txt",sep = " ",header = T,stringsAsFactors = F,check.names = F)
rownames(metadata)=metadata$ID
metadata=metadata[rownames(metadata) %in% rownames(expression),]
coupling=read.table("coupling_file.txt",sep = "\t",header = F,stringsAsFactors = F)
coupling=coupling[coupling$V2 %in% rownames(expression),]

# set covariates
metadata$Location[metadata$Location!="ileum"]="colon"
metadata$Inflammation[metadata$Inflammation=="NI"]=0
metadata$Inflammation[metadata$Inflammation=="I"]=1
metadata$Location[metadata$Location=="ileum"]=0
metadata$Location[metadata$Location=="colon"]=1
metadata$Diagnosis[metadata$Diagnosis=="CD"]=0
metadata$Diagnosis[metadata$Diagnosis=="UC"]=1
metadata$Inflammation=as.numeric(metadata$Inflammation)
metadata=merge(metadata,coupling,by.x="row.names",by.y = "V2",all=F)

# import genotype file
genotype=read.table("SNP.sig.dosage.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
genotype$`-`=NULL
genotype=as.data.frame(t(genotype))

# import 190 significant inflammation-dependent cis-eQTLs
signals=read.table("Merged.FDR.sig.independent.txt",header = T,sep = "\t",stringsAsFactors = F)

# import IBS matrix using as a random effect
matrx=read.table("IBS.mibs",header = F,stringsAsFactors = F)
matrx.id=read.table("IBS.mibs.id.txt",header = F,stringsAsFactors = F)
rownames(matrx)=matrx.id$V2
colnames(matrx)=matrx.id$V2

# ========================== compare results from GEMMA and lme4qtl ===========================

# we compared Z score of the interaction term in the interaction model from GEMMA and lme4qtl
lme4qtl.results=read.table("Models.compare.txt",sep = "\t",header = T,stringsAsFactors = F)
gemma.results=read.table("Merged.FDR.sig.independent.txt",sep = "\t",header = T,stringsAsFactors = F)
gemma.results$z=gemma.results$Beta/bb$SE

mm=merge(lme4qtl.results[,c("Probe","SNP","IBS.interaction.t")],gemma.results[,c("ExpressionGene","z")],
         by.x="Probe",by.y="ExpressionGene",all=F)
for(i in 1:nrow(mm)){
  if(mm$z[i]<0){
    mm$IBS.interaction.t[i]=-abs(mm$IBS.interaction.t[i])
  }else if(mm$z[i]>0){
    mm$IBS.interaction.t[i]=abs(mm$IBS.interaction.t[i])
  }
}

ggplot(mm, aes(x=z, y=IBS.interaction.t)) + 
  geom_point()+xlab("Z score of interaction term from GEMMA")+ylab("Z score of interaction term from lme4qtl")
ggsave("Plot/GEMMA.compare.pdf",width = 6,height = 5)

# ========================== compare three models to adjust for random effect ===========================

# three models compare
# 1)	Interaction model without adjusting for repeat measurements:
# 'lm' function: gene = intercept + PCs (1,3~18) + SNP + Inflammation + SNP*Inflammation 
# 2)	Interaction model using patients ID vector (1|ID) for repeat measurements:
# 'lmer' function: gene = intercept + PCs (1,3~18) + SNP + Inflammation + SNP*Inflammation + (1|ID)
# 3)	Interaction model using IBS matrix for repeat measurements:
# 'relmatLmer' function: gene = intercept + PCs (1,3~18) + SNP + Inflammation + SNP*Inflammation + (IBS matrix)

resut_matrix = foreach(i=1:nrow(signals),.combine = rbind) %do%  {
  tryCatch({
    cat(green(i,"calculating +++","\n"))
    tmp.gene=signals$ExpressionGene[i]
    tmp.snp=signals$rsID[i]
    tmp.genotype=genotype[,tmp.snp,drop=F]
    tmp.genotype[tmp.genotype==-1]=NA
    tmp.genotype=na.omit(tmp.genotype)
    tmp.table=merge(tmp.genotype,metadata[,c("ID","Inflammation","V1","UMCG_ID")],by.x="row.names",by.y="V1",all=F)
    tmp.table=merge(tmp.table,expression[,tmp.gene,drop=F],by.x="ID",by.y = "row.names",all=F)
    tmp.matrix=matrx[rownames(matrx) %in% tmp.table$Row.names,]
    tmp.matrix=tmp.matrix[,colnames(matrx) %in% tmp.table$Row.names]
    tmp.table=tmp.table[tmp.table$Row.names %in% rownames(tmp.matrix),]
    tmp.table=tmp.table[order(tmp.table$Row.names),]
    tmp.matrix=tmp.matrix[order(rownames(tmp.matrix)),]
    tmp.matrix=tmp.matrix[,order(colnames(tmp.matrix))]
    tmp.matrix=as.matrix(tmp.matrix)
    
    fit0<-relmatLmer(tmp.table[,6]~tmp.table[,3]+tmp.table[,4]+tmp.table[,3]*tmp.table[,4]+(1|Row.names),data=tmp.table,relmat = list(Row.names = tmp.matrix))
    fit1<-lmer(tmp.table[,6]~tmp.table[,3]+tmp.table[,4]+tmp.table[,3]*tmp.table[,4]+(1|UMCG_ID),data=tmp.table)
    fit2<-lm(tmp.table[,6]~tmp.table[,3]+tmp.table[,4]+tmp.table[,3]*tmp.table[,4],data=tmp.table)
    coefs.fit0 <- data.frame(coef(summary(fit0)))
    coefs.fit1 <- data.frame(coef(summary(fit1)))
    coefs.fit2 <- data.frame(coef(summary(fit2)))
    coefs.fit0$p.z <- 2 * (1 - pnorm(abs(coefs.fit0$t.value)))
    coefs.fit1$p.z <- 2 * (1 - pnorm(abs(coefs.fit1$t.value)))
    return.string = data.frame(SNP=tmp.snp,Probe=tmp.gene,
                               lm.snp.Pvalue=coefs.fit2$Pr...t..[2],
                               lm.snp.beta=coefs.fit2$Estimate[2],
                               lm.snp.t=coefs.fit2$t.value[2],
                               lm.inflammation.Pvalue=coefs.fit2$Pr...t..[3],
                               lm.inflammation.beta=coefs.fit2$Estimate[3],
                               lm.inflammation.t=coefs.fit2$t.value[3],
                               lm.interaction.Pvalue=coefs.fit2$Pr...t..[4],
                               lm.interaction.beta=coefs.fit2$Estimate[4],
                               lm.interaction.t=coefs.fit2$t.value[4],
                               IBS.snp.Pvalue=coefs.fit0$p.z[2],
                               IBS.snp.beta=coefs.fit0$Estimate[2],
                               IBS.snp.t=coefs.fit0$t.value[2],
                               IBS.inflammation.Pvalue=coefs.fit0$p.z[3],
                               IBS.inflammation.beta=coefs.fit0$Estimate[3],
                               IBS.inflammation.t=coefs.fit0$t.value[3],
                               IBS.interaction.Pvalue=coefs.fit0$p.z[4],
                               IBS.interaction.beta=coefs.fit0$Estimate[4],
                               IBS.interaction.t=coefs.fit0$t.value[4],
                               ID.snp.Pvalue=coefs.fit1$p.z[2],
                               ID.snp.beta=coefs.fit1$Estimate[2],
                               ID.snp.t=coefs.fit1$t.value[2],
                               ID.inflammation.Pvalue=coefs.fit1$p.z[3],
                               ID.inflammation.beta=coefs.fit1$Estimate[3],
                               ID.inflammation.t=coefs.fit1$t.value[3],
                               ID.interaction.Pvalue=coefs.fit1$p.z[4],
                               ID.interaction.beta=coefs.fit1$Estimate[4],
                               ID.interaction.t=coefs.fit1$t.value[4]
)},error=function(e){})
  
}

