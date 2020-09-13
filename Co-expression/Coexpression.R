library(foreach)
args=commandArgs(T)
i=args[1]
expression=read.table("TMM_expression.table.Log2Transformed.ProbesCentered.SamplesZTransformed.18PCAsOverSamplesRemoved.txt",
                      row.names = 1,header = T,stringsAsFactors = F,sep = "\t")

CoExpression = foreach(n=1:nrow(expression),.combine = rbind) %do%  {
  mm=cor.test(as.matrix(expression[rownames(expression)==i,]),as.matrix(expression[n,]),method = "spearman")
  p=mm$p.value
  coef=mm$estimate
  print(n)

  return.string = data.frame(TargetGene = i, CoexpressGene=rownames(expression)[n],
                             Coeffecient = coef,Pvalue = p)
}

CoExpression$FDR=p.adjust(CoExpression$Pvalue)
write.table(CoExpression, file = paste(i,"coexpression.txt",sep = "."),
            sep = '\t',row.names = F,quote = F)
