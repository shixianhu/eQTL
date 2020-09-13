- we used GEMMA for an interaction model (gene expression = SNP + inflammation + SNP * inflammation) to calculate inflammation-dependent cis-eQTL

- however, we only got summary statistics of interaction term from GEMMA by default. "In this case, for each SNP in turn, GEMMA will fit a linearmixed model that controls both the SNP main effect and environmental main effect, while testing for the interaction effect." (https://www.xzlab.org/software/GEMMAmanual.pdf)

- to get full summary statistics of term 'SNP' and 'inflammation', we used lme4qtl R package to re-calculate 190 significant (FDR interaction <0.05) from GEMMA

- we compared the results derived from GEMMA and lme4qtl, the Z score of the interaction term (SNP * inflammation). See GEMMA.lme4qtl.interaction.model.results.compare.pdf

- in addition, we also compared three models to adjust for random effect. See Two_random_effect.compare.Z.pdf and No_random_effect.sompare.Z.pdf

```
 1)	Interaction model without adjusting for repeat measurements:
 
 'lm' function: gene = intercept + PCs (1,3~18) + SNP + Inflammation + SNP*Inflammation
 
 2)	Interaction model using patients ID vector (1|ID) for repeat measurements:
 
 'lmer' function: gene = intercept + PCs (1,3~18) + SNP + Inflammation + SNP*Inflammation + (1|ID)
 
 3)	Interaction model using IBS matrix for repeat measurements:
 
 'relmatLmer' function: gene = intercept + PCs (1,3~18) + SNP + Inflammation + SNP*Inflammation + (IBS matrix)
 
```
