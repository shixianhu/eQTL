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

# Part 1. cis-eQTL analysis


*step 1. Normalization and log transformation*
---
- Use the expression matrix with included samples to run the TMM normalization.
```
An assumption of TMM is the majority of the genes are not differentially expressed. 

The main aim in TMM normalization is to account for library size variation between samples of interest, accounting for the fact that some extremely differentially expressed genes would impact negatively the normalization procedure.

A trimmed mean is the average after removing the upper and lower x% of the data.
```
We use edgeR to run TMM normalization.
```
library(edgeR)
library(limma)
library(RColorBrewer)
library(mixOmics)
library(VennDiagram)
library(HTSFilter)

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

---> output: Pheno.txt Reordered.phenotype.txt

vim Reordered.phenotype.txt and add "-"
```


*step 3.2. eQTL analysis - Generate relatedness file*
---

- This step is to consider kinship as a random effect in mixed linear model.

```
ml plink

plink --bfile genoytpe.plink --distance ibs Kinship/IBS

---> output: IBS.cluster; IBS.log; IBD.mibs; IBD.mibs.ID;, IBS.nosex

```


*step 3.3. eQTL analysis - Loop for each expression probe using GEMMA*
---

- Note for all *gcc* dependencies. 

```
ml plink

awk '{print $2}' All_pairs.txt | sort | uniq > Probe.txt

export  LD_LIBRARY_PATH=/home/umcg-hushixian/gemma/gcc-5.4.0-lib-3x53yv4v144c9xp0/lib

cat Probe.txt | while read line

do

grep -w $line All_pairs.txt | awk '{print $1}' > tmp.snp.txt
plink --bfile genotype.plink --extract tmp.snp.txt --make-bed --out tmp.analysis
awk -v col=$line 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' Reordered.phenotype.txt > tmp.expression.txt
sed -i '1d' tmp.expression.txt 

~/gemma/bin/gemma -bfile tmp.analysis \
-p tmp.expression.txt \
-km 1 -k Kinship/IBS.mibs \
-lmm 4 -o $line.outcome \
-miss 0.99

rm tmp* 
# this removing is very important TAKE CARE !!!!!!

done
```
- Add 𝑆𝑁𝑃 × 𝑖𝑛𝑓𝑙𝑎𝑚𝑚𝑎𝑡𝑖𝑜𝑛 interaction term in GEMMA
```
~/gemma/bin/gemma -bfile tmp.analysis \
-p tmp.expression.txt \
-gxe covariate.txt \
-km 1 -k Kinship/IBS.mibs \
-lmm 4 -o $line.outcome \
-miss 0.99
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
