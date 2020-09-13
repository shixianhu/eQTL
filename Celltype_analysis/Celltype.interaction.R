# we used xCell to predict cell type enrichment using RNA-seq data
# we assessed SNP * celltype enrichment effect for the 190 inflammation-dependent cis-eQLs

library(xCell)
library(foreach)
library(lme4qtl)
library(nlme)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
library(randomcoloR)
library(crayon)
library(ggforce)

# ================================ cell type deconvolution ================================ 

# import gene expression table 
count=read.table("ExpressionTable.txt",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
dgeFull <- calcNormFactors(dgeFull, method="TMM")    
expression_norm=cpm(dgeFull,log = TRUE, prior.count = 1)  
expression_norm=as.data.frame(t(expression_norm))

# import phenotype data
metadata=read.table("complete_pheno_table_pil280.txt",sep = " ",header = T,stringsAsFactors = F,check.names = F)
rownames(metadata)=metadata$ID
metadata$Location[metadata$Location!="ileum"]="colon"
metadata$Inflammation[metadata$Inflammation=="NI"]=0
metadata$Inflammation[metadata$Inflammation=="I"]=1
metadata$Location[metadata$Location=="ileum"]=0
metadata$Location[metadata$Location=="colon"]=1
metadata$Diagnosis[metadata$Diagnosis=="CD"]=0
metadata$Diagnosis[metadata$Diagnosis=="UC"]=1
inflammation=metadata[metadata$Inflammation==1,]
control=metadata[metadata$Inflammation==0,]

# deconvolution analysis on all 64 cell types in XCell
deconvolution=xCellAnalysis(expression_norm)
deconvolution=as.data.frame(deconvolution)
deconvolution=as.data.frame(t(deconvolution))
deconvolution=deconvolution[,c("cDC",	"Macrophages M1",	"NK cells",	"pDC",	"Macrophages M2",
                               "CD4+ naive T-cells",	"CD4+ Tcm",	"CD8+ naive T-cells",	"CD8+ Tcm",	
                               "Tgd cells",	"Th2 cells",	"Tregs",	"Th1 cells",
                               "NKT",	"CD8+ Tem",	"CD4+ Tem",	"Class-switched memory B-cells",
                               "Plasma cells",	"naive B-cells",	"Memory B-cells",	"Basophils",
                               "Mast cells",	"Neutrophils",	"Eosinophils",	"Endothelial cells",
                               "Epithelial cells",	"Fibroblasts",	"Smooth muscle"),drop=F]

# ================================ cell type enrichment compare ================================ 

# compare cell type enrichemnt between inflamed and non-inflamed biopsies
compare=matrix(nrow = ncol(deconvolution),ncol = 2)
compare=as.data.frame(compare)
colnames(compare)=c("CellType","Pvalue")
for(i in 1:ncol(deconvolution)){
  cell=colnames(deconvolution)[i]
  inflamed=deconvolution[rownames(deconvolution) %in% inflammation$ID,]
  noninf=deconvolution[rownames(deconvolution) %in% control$ID,]
  mm=wilcox.test(inflamed[,cell],noninf[,cell])
  Pvalue=mm$p.value
  compare$CellType[i]=cell
  compare$Pvalue[i]=Pvalue
}
compare$FDR=p.adjust(compare$Pvalue)
compare=compare[order(compare$FDR),]

compare_all=data.frame(Score=NA,Group=NA,Cell=NA,Sample=NA)
for(i in 1:nrow(compare)){
  cell=compare$CellType[i]
  sub.in=inflamed[,cell,drop=F]
  sub.non=noninf[,cell,drop=F]
  sub.in$Group="Inflammation"
  sub.non$Group="Non_inflammation"
  sub.in$Cell=cell
  sub.non$Cell=cell
  colnames(sub.in)[1]="Score"
  colnames(sub.non)[1]="Score"
  sub.data=rbind(sub.in,sub.non)
  sub.data$Sample=rownames(sub.data)
  rownames(sub.data)=NULL
  compare_all=rbind(compare_all,sub.data)
}
compare_all=na.omit(compare_all)
compare_all$Cell=factor(compare_all$Cell,levels = (unique(compare_all$Cell)))

ggplot(compare_all, aes(compare_all$Cell,compare_all$Score)) + 
  geom_boxplot(aes(fill=Group),position=position_dodge(1),outlier.shape = NA)+
  theme_bw()+
  scale_fill_manual(values=c("#339900","#FF9900"))+
  scale_color_manual(values=c("#339900","#FF9900"))+
  guides(color=FALSE)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
#  theme(axis.text.x = element_text(size = 6,hjust = 0))+
  xlab("CellTypes")+ylab("EnrichmentScore")
  #+coord_flip()
ggsave("CellType.In_NonIn.pdf",width = 8,height = 5)

# ================================ cell type * SNP interaction analysis ================================ 

# import gene expression table after confounder adjusting
expression=read.table("ExpressionTable.18PCs.txt",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
all_cell=deconvolution
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
all_cell = apply(all_cell,2,invrank)
all_cell=as.data.frame(all_cell,stringsAsFactors = F)

# import 190 significant (FDR <0.05) inflammation-dependent cis-eQTLs
eInteraction=read.table("Merged.FDR.sig.independent.txt",
                        sep = "\t",header = T,stringsAsFactors = F)
eSNP=read.table("SNP.sig.dosage.txt",
                sep = "\t",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
eSNP=as.data.frame(t(eSNP))
coup=read.table("coupling_file.txt")
coup=coup[coup$V1 %in% rownames(eSNP),]
eSNP=eSNP[rownames(eSNP) %in% coup$V1,]
coup=coup[order(coup$V1),]
eSNP=eSNP[order(rownames(eSNP)),]
rownames(eSNP)=coup$V2
eSNP=merge(eSNP,metadata,by="row.names")
rownames(eSNP)=eSNP$Row.names
eSNP$Row.names=NULL

expression=as.data.frame(t(expression))
expression=expression[rownames(expression) %in% coup$V2,]
expression=expression[order(rownames(expression)),]
expression=merge(expression,metadata[,"UMCG_ID",drop=F],by="row.names")
rownames(expression)=expression$Row.names
expression$Row.names=NULL

# linear mixed interaction model: gene expression = SNP + celltype + SNP * celltype + IBS
eInteraction=merge(eInteraction,annotation,by.x="ExpressionGene",by.y = "Gene",all=F)
eInteraction=eInteraction[,c("rsID","id")]

matrx=read.table("IBS.mibs",header = F,stringsAsFactors = F)
matrx.id=read.table("IBS.mibs.id.txt",header = F,stringsAsFactors = F)
rownames(matrx)=matrx.id$V2
colnames(matrx)=matrx.id$V2

coupling=read.table("coupling_file.txt",sep = "\t",stringsAsFactors = F,header = F)
celltype=data.frame(SNP=NA,Gene=NA,Cell=cell,
                    IBS.snp.Pvalue=NA,
                    IBS.snp.beta=NA,
                    IBS.snp.t=NA,
                    IBS.enrichscore.Pvalue=NA,
                    IBS.enrichscore.beta=NA,
                    IBS.enrichscore.t=NA,
                    IBS.interaction.Pvalue=NA,
                    IBS.interaction.beta=NA,
                    IBS.interaction.t=NA)

for(n in 1:28){

  print(n)
  cell=colnames(all_cell)[n]
  celleffect=foreach(i=1:nrow(eInteraction),.combine = rbind) %do%  {
    cat(green(i,"+++++++++++++++++++","\n"))
  snp=eInteraction$rsID[i]
  probe=eInteraction$id[i]
  genotype=eSNP[,c(snp,"UMCG_ID"),drop=F]
  genotype[genotype==-1]=NA
  genotype=na.omit(genotype)
  phenotype=expression[,c(probe,"UMCG_ID"),drop=F]
  phenotype=merge(phenotype,genotype,by="row.names",all=F)
  rownames(phenotype)=phenotype$Row.names
  phenotype$UMCG_ID.x=as.factor(phenotype$UMCG_ID.x)
  aa=phenotype[,c(probe,snp,"UMCG_ID.x")]
  aa=merge(aa,all_cell[,cell,drop=F],by="row.names")
  aa=merge(aa,coupling,by.x="Row.names",by.y="V2")
  
  tmp.matrix=matrx[rownames(matrx) %in% aa$V1,]
  tmp.matrix=tmp.matrix[,colnames(matrx) %in% aa$V1]
  aa=aa[aa$V1 %in% rownames(tmp.matrix),]
  aa=aa[order(aa$V1),]
  tmp.matrix=tmp.matrix[order(rownames(tmp.matrix)),]
  tmp.matrix=tmp.matrix[,order(colnames(tmp.matrix))]
  tmp.matrix=as.matrix(tmp.matrix)
  aa[,3]=as.numeric(as.character(aa[,3]))
  aa$gene=aa[,2]

  fit0<-relmatLmer(gene ~ aa[,3]+aa[,5]+aa[,3]*aa[,5]+(1|Row.names),relmat = list(Row.names = tmp.matrix),data=aa)

  coefs.fit0 <- data.frame(coef(summary(fit0)))
  coefs.fit0$p.z <- 2 * (1 - pnorm(abs(coefs.fit0$t.value)))
  return.string = data.frame(SNP = snp, Gene = probe,Cell=cell,
                             IBS.snp.Pvalue=coefs.fit0$p.z[2],
                             IBS.snp.beta=coefs.fit0$Estimate[2],
                             IBS.snp.t=coefs.fit0$t.value[2],
                             IBS.enrichscore.Pvalue=coefs.fit0$p.z[3],
                             IBS.enrichscore.beta=coefs.fit0$Estimate[3],
                             IBS.enrichscore.t=coefs.fit0$t.value[3],
                             IBS.interaction.Pvalue=coefs.fit0$p.z[4],
                             IBS.interaction.beta=coefs.fit0$Estimate[4],
                             IBS.interaction.t=coefs.fit0$t.value[4])
}
  celltype=rbind(celltype,celleffect)
  
}

celltype.clean=na.omit(celltype)
celltype.clean$IBS.interaction.FDR=p.adjust(celltype.clean$IBS.interaction.Pvalue)
celltype.clean$IBS.snp.FDR=p.adjust(celltype.clean$IBS.snp.Pvalue)
celltype.clean$IBS.enrichscore.FDR=p.adjust(celltype.clean$IBS.enrichscore.Pvalue)

