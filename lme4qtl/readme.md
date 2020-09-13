- we used GEMMA for an interaction model (gene expression = SNP + inflammation + SNP * inflammation) to calculate inflammation-dependent cis-eQTL

- however, we only got summary statistics of interaction term from GEMMA by default

- to get full summary statistics of term 'SNP' and 'inflammation', we used lme4qtl R package to re-calculate 190 significant (FDR interaction <0.05) from GEMMA

- in addition, we also compared the results derived from GEMMA and lme4qtl, and compared three models to adjust for random effect
