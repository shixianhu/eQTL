library(foreach)
args=commandArgs(T)
i=args[1]
expression=read.table("TMM_expression.table.Log2Transformed.ProbesCentered.SamplesZTransformed.18PCAsOverSamplesRemoved.txt",
                      row.names = 1,header = T,stringsAsFactors = F,sep = "\t")

coupling=read.table("coupling_file.txt",sep = "\t",header = F,stringsAsFactors = F)
coupling=coupling[coupling$V2 %in% rownames(expression),]

matrx=read.table("IBS.mibs",header = F,stringsAsFactors = F)
matrx.id=read.table("IBS.mibs.id.txt",header = F,stringsAsFactors = F)
rownames(matrx)=matrx.id$V2
colnames(matrx)=matrx.id$V2
tmp.table=as.data.frame(t(expression[,i,drop=F]))

CoExpression = foreach(n=1:nrow(expression),.combine = rbind) %do%  {
  tmp.gene=rownames(expression)[n]
  tmp.gene=as.data.frame(t(tmp.gene))
  tmp.table=merge(tmp.table,expression[,tmp.gene,drop=F],by = "row.names",all=F)
  fit0<-relmatLmer(tmp.table[,1]~tmp.table[,2]+(1|Row.names),data=tmp.table,relmat = list(Row.names = tmp.matrix)
  coefs.fit0$p.z <- 2 * (1 - pnorm(abs(coefs.fit0$t.value)))
  return.string = data.frame(IBS.gene.Pvalue=coefs.fit0$p.z[2],
                             IBS.gene.beta=coefs.fit0$Estimate[2],
                             IBS.gene.t=coefs.fit0$t.value[2])

}

CoExpression$FDR=p.adjust(CoExpression$IBS.gene.Pvalue)
write.table(CoExpression, file = paste(i,"coexpression.txt",sep = "."),
            sep = '\t',row.names = F,quote = F)
