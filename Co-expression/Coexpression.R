library(foreach)
args=commandArgs(T)
i=args[1]
expression=read.table("TMM_expression.table.Log2Transformed.ProbesCentered.SamplesZTransformed.18PCAsOverSamplesRemoved.txt",
                      row.names = 1,header = T,stringsAsFactors = F,sep = "\t")
metadata=read.table("complete_pheno_table_pil280.txt",sep = " ",header = T,stringsAsFactors = F,check.names = F)
rownames(metadata)=metadata$ID
metadata=metadata[rownames(metadata) %in% rownames(expression),]
coupling=read.table("coupling_file.txt",sep = "\t",header = F,stringsAsFactors = F)
coupling=coupling[coupling$V2 %in% rownames(expression),]

metadata$Location[metadata$Location!="ileum"]="colon"
metadata$Inflammation[metadata$Inflammation=="NI"]=0
metadata$Inflammation[metadata$Inflammation=="I"]=1
metadata$Location[metadata$Location=="ileum"]=0
metadata$Location[metadata$Location=="colon"]=1
metadata$Diagnosis[metadata$Diagnosis=="CD"]=0
metadata$Diagnosis[metadata$Diagnosis=="UC"]=1
metadata$Inflammation=as.numeric(metadata$Inflammation)
metadata=merge(metadata,coupling,by.x="row.names",by.y = "V2",all=F)

matrx=read.table("IBS.mibs",header = F,stringsAsFactors = F)
matrx.id=read.table("IBS.mibs.id.txt",header = F,stringsAsFactors = F)
rownames(matrx)=matrx.id$V2
colnames(matrx)=matrx.id$V2

CoExpression = foreach(n=1:nrow(expression),.combine = rbind) %do%  {
  fit0<-relmatLmer(tmp.table[,6]~tmp.table[,3]+(1|Row.names),data=tmp.table,relmat = list(Row.names = tmp.matrix)
  coefs.fit0$p.z <- 2 * (1 - pnorm(abs(coefs.fit0$t.value)))
  return.string = data.frame(IBS.gene.Pvalue=coefs.fit0$p.z[2],
                             IBS.gene.beta=coefs.fit0$Estimate[2],
                             IBS.gene.t=coefs.fit0$t.value[2])

}

CoExpression$FDR=p.adjust(CoExpression$IBS.gene.Pvalue)
write.table(CoExpression, file = paste(i,"coexpression.txt",sep = "."),
            sep = '\t',row.names = F,quote = F)
