---
creator: "Shixian"
date: "03/14/2019"
RNA-seq Data: "171 individuals; 299 biopsy"
Genomic Data: "171 individuals; WES+GSA"
Sample Excluded: "8CD"
Sample Included: "185CD + 106UC+IBDU"
---

# eQTL analysis based on mucosal biopsy RNA-seq in IBD

This project is to identify the eQTL effect in context of inflammation and non-inflammation in mucosal biopsy in IBD



*Models used (generalized linear mixed model):*
---

 - Differentially gene expression (DGE) analysis
```
DGE of inflammation:
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 = 𝛼 + 𝛽𝑆𝑁𝑃 + (1+3~18)𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠(IBS matrix) + 𝜀

DGE of disease diagnosis:
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 = 𝛼 + 𝛽𝑆𝑁𝑃 + (2~18)𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠(IBS matrix) + 𝜀

DGE of biopsy location:
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 = 𝛼 + 𝛽𝑆𝑁𝑃 + (2~18)𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠(IBS matrix) + 𝜀
```

 - cis-eQTL analysis 
```
General cis-eQTLs
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 = 𝛼 + 𝛽𝑆𝑁𝑃 + 18𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠(IBS matrix) + 𝜀

Inflammation-dependent cis-eQTLs
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 = 𝛼 + 𝛽𝑆𝑁𝑃 + (1+3~18)𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠(IBS matrix) + 𝛽𝑖𝑛𝑓𝑙𝑎𝑚𝑚𝑎𝑡𝑖𝑜𝑛 + 𝛽𝑆𝑁𝑃 × 𝑖𝑛𝑓𝑙𝑎𝑚𝑚𝑎𝑡𝑖𝑜𝑛 + 𝜀
```

 - Pair-wise gene co-expression analysis
 ```
 𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 = 𝛼 + 𝛽gene + 18𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠(IBS matrix) + 𝜀
 ```

 - Cell-type specific analysis 
 ```
 𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛 = 𝛼 + 𝛽𝑆𝑁𝑃 + 18𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠(IBS matrix) +𝛽celltype_enrich_socre + 𝛽𝑆𝑁𝑃 × celltype_enrich_socre + 𝜀
 ```
 
 
*RNA-seq data QC*
---
```
1. Reads alignment percentage < 90%; mapped reads < 30 million.     ---> 4 samples are removed
2. Duplicate samples check                                          ---> 2 samples are removed
3. Outliers from expression data (PCA check).                       ---> 2 samples are removed
4. Phenotype mismatch and lightly-inflamed samples                  ---> 11 samples are removed
```

# Part 1. differential gene expression (DGE) analysis

*step 1. Generate group (CD/UC, colon/ileum, inflamed/non-inflamed) file*
---

- GEMMA is used for eQTL anlysis instead of DGE analysis, however, we can replace the genotype bimbam file with phenotype file in GEMMA.

```
For example, we code sample inflammation as 0 (inflamed samples) and 1 (non-inflamed samples) in Inflammation.bimbam

rs00000 A T 1 1 0 1 1 1 1 0 1 1 0 1 0 1 1 1 

where rs00000, A and T are random but nessacery for generating bimbam file to run in GEMMA. 
```

*step 2. Run GEMMA for DGE annlysis*
---

```
gemma-0.98-linux-static -g Inflammation.bimbam -p gene.expression.txt -lmm 4 -km 1 -k IBS.mibs –o DGE.out

```


# Part 2. cis-eQTL analysis


*step 1. Normalization and log transformation*
---

- Use the expression matrix to run the TMM normalization (edge R).

```
library(edgeR)
library(limma)

count=read.table("ExpressionTable.txt",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
dgeFull <- DGEList(count, remove.zeros = TRUE)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
timmed=cpm(dgeFull,log = TRUE, prior.count = 1) 
timmed=as.data.frame(timmed)
timmed=data.frame(rownames(timmed),timmed,check.names = F)
colnames(timmed)[1]="probe"
```

*step 2. Remove PCs*
---

- By adjusting for a set of PCs, we try to remove confouders effects.

```
# get PCs of gene expression

expression=timmed[,2:ncol(timmed)]
rownames(expression)=rownames(timmed)
pca=prcomp(timmed,scale = TRUE)
eigenvalue=get_eig(pca)
ind <- get_pca_ind(pca)
pca_matrix=as.data.frame(ind$coord)
pca_matrix=pca_matrix[1:18,]

# corrected for targeted PCs (for example, the first 18PCs)

pca_matrix=pca_matrix[1:18,]
corrected_data = apply(expression,2,function(x){
  x.resid = resid(lm(x ~ ., data=pca_matrix))
  return(x.resid)
})
```


*step 3.1. eQTL analysis - Match expression data to genotype data*
---

 - All_pairs.txt (gene-SNP pairs)
 - plink files (genotype file)
 - Normalized.txt (gene expression data after removing PCs)
 - coupling file (connect biopsy ID to WES+genotype ID)

```
# match the samples orders of phenotype file and genotype file

Rscript Phenotype.Prepare.R Normalized.txt Plink.fam

vim Reordered.phenotype.txt and add "-"
```


*step 3.2. eQTL analysis - Generate relatedness file*
---

- This step is to consider kinship as a random effect in mixed linear model.

```
ml plink

plink --bfile genoytpe.plink --distance ibs Kinship/IBS

```

*step 3.3. eQTL analysis - Run GEMMA*
---

- Intestinal cis-eQTL model in GEMMA
```
gemma-0.98-linux-static -bfile genotype.plink -p gene.expression.txt -km 1 -k IBS.mibs -lmm 4 -o $line.outcome 

```

- Add 𝑆𝑁𝑃 × 𝑖𝑛𝑓𝑙𝑎𝑚𝑚𝑎𝑡𝑖𝑜𝑛 (gxe) in an interaction model in GEMMA
```
gemma-0.98-linux-static -bfile genotype.plink -p tmp.expression.txt -gxe covariate.txt -km 1 -k IBS.mibs -lmm 4 -o $line.outcome 

```

*step 3.4. eQTL analysis - Merging results*
---

- To merge all eQTL results of each expression gene (*eg. ENSG00000072135*) in output folder.

```
echo -e "ExpressionGene\tChr\trsID\tPos\tMissingSample\tAllele1\tAllele0\tAllelFre\tBeta\tSE\tlogl_H1\tl_remle\tl_mle\tp_wald\tp_lrt\tp_score" > Merge.assoc.txt

for i in /groups/umcg-gastrocol/tmp04/Inflamed_eQTL/Previous_process/GEMMA_mixed_model/eQTL_CD/output/*.outcome.assoc.txt

do

name=$(basename $i)
expression=${name%.outcome.assoc.txt}
awk -v var="$expression" '{OFS="\t"}{if (NR!=1) print var,$0}' $i >> Merge.assoc.txt

done
