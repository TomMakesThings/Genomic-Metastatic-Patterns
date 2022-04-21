# Significance
Significance between primary and metastasis samples for FGA, WGD, TMB and TMB was quantified using the Mann-Whitney U test (`significant_subtypes.csv`). Q-values were calculated by adjusting p-values for false discovery rate and a q-value < 0.05 was considered statistically significant.

# TMB High and WGD
The fraction of samples per subtype for primary and metastatic samples with a tumor mutational burden >= 10 mutations per megabase was calculated. This was combined with WGD and number of patients in file `tmb_high_and_wgd.csv`.

# Arm-level CNAs
The frequency of arm-level copy number alterations between primary and metastatic tumors was calculated as the fraction of AMP / DEL per chromsome arm, e.g. 12p, predicted by ASCETS. This is saved in `subtype_arm_level_fraction_alteration.csv` and is required for the heatmap in plot 2B.

# FGA and TMB
Our calculated values are found [here](https://github.com/TomMakesThings/Genomics-II-Group/blob/main/Plot_2/FGA/calculated_TMB_and_FGA.csv).
