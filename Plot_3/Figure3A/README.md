# Spearman's for Metastatic Burden
The genomic features FGA and TMB were considered significantly correlated with metastatic burden if:
* Spearman’s correlation between the FGA/TMB and the number of metastatic sites is significant with a q-value < 0.05
  - Where q-value is p-value adjusted for false discovery rate
* The feature is signficant with p-value < 0.05 when used as a predictive variable in a linear regression model
  - Before running regression, power transformations were performed on FGA and TMB with Tukey’s ladder of powers to harmonize their distributions
  - Then they were normalised between 0 - 1
  - Multivariable linear regression was run to model metastatic site count using the feature and sample type as coefficients, e.g. `MET_SITE_COUNT ~ SAMPLE_TYPE + FGA`
  - This condition is required so that the ratio of primary to metastatic samples doesn't act as a confounding factor as it has an association with metastatic burden
