library(tibble)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Set q-value threshold
q_value_threshold <- 0.05

# Open calculated TMB and FGA data
samples_data <- read.csv(file = "Plot_2/FGA/calculated_TMB_and_FGA.csv",
                         header = TRUE, fill = TRUE)
cna_data <- read.csv(file = "Plot_2/aSCNAs/sample_arm_level_cna.csv",
                         header = TRUE, fill = TRUE)

# Add calculated WGD
samples_data$Our_WGD <- cna_data$ARM_WGD
# Get subtypes
all_subtypes <- unique(samples_data$SUBTYPE)

# Record TMB high for primary and metastasis samples
tmb_high_data <- data.frame()

# Record p-values from Mann-Whitney U test of all subtypes
fga_pvalues <- c()
wgd_pvalues <- c()
tmb_pvalues <- c()
tmb_high_pvalues <- c()

for (subtype in all_subtypes) {
  # Get samples for the subtype
  subtype_rows <- samples_data[which(samples_data$SUBTYPE == subtype), ]
  
  # Find primary and metastasis samples
  primary_subtype <- subtype_rows[which(subtype_rows$SAMPLE_TYPE == "Primary"), ]
  metastasis_subtype <- subtype_rows[which(subtype_rows$SAMPLE_TYPE == "Metastasis"), ]
  
  # Get samples with TMB >= 10 mut/Mb
  primary_tmb_high <- primary_subtype[which(primary_subtype$Our_TMB >= 10), ]
  metastasis_tmb_high <- metastasis_subtype[which(metastasis_subtype$Our_TMB >= 10), ]
  
  # Calculate fraction of samples with high metastatic burden
  primary_tmb_frac <- nrow(primary_tmb_high) / nrow(primary_subtype)
  metastasis_tmb_frac <- nrow(metastasis_tmb_high) / nrow(metastasis_subtype)
  
  # Record the TMB high fraction
  tmb_high_data <- rbind(tmb_high_data,
                         data.frame(SUBTYPE = subtype,
                                    SAMPLE_TYPE = "Primary",
                                    TMB_HIGH_FRACTION = primary_tmb_frac))
  tmb_high_data <- rbind(tmb_high_data,
                         data.frame(SUBTYPE = subtype,
                                    SAMPLE_TYPE = "Metastasis",
                                    TMB_HIGH_FRACTION = metastasis_tmb_frac))
  
  # Mann-Whitney U test to compare primary and metastasis samples
  man_whitney_fga <- wilcox.test(primary_subtype$Our_FGA, metastasis_subtype$Our_FGA)
  man_whitney_wgd <- wilcox.test(primary_subtype$Our_WGD, metastasis_subtype$Our_WGD)
  man_whitney_tmb <- wilcox.test(primary_subtype$Our_TMB, metastasis_subtype$Our_TMB)
  fga_pvalues <- c(fga_pvalues, man_whitney_fga$p.value)
  wgd_pvalues <- c(wgd_pvalues, man_whitney_fga$p.value)
  tmb_pvalues <- c(tmb_pvalues, man_whitney_tmb$p.value)
  
  if (primary_tmb_frac == 0 | metastasis_tmb_frac == 0) {
    # Set highest p-value if no high TMB samples found
    tmb_high_pvalues <- c(tmb_high_pvalues, 1)
  } else {
    # Otherwise perform Mann-Whitney U test
    man_whitney_tmb_high <- wilcox.test(primary_tmb_high$Our_TMB, metastasis_tmb_high$Our_TMB)
    tmb_high_pvalues <- c(tmb_high_pvalues, man_whitney_tmb_high$p.value)
  }
}

# Convert from p-values to q-values and return which subtypes are significant
findSigniciant <- function(pvalues, subtype_list, threshold) {
  qvalues <- p.adjust(tmb_high_pvalues, method = "fdr")
  
  return(qvalues < threshold)
}

# Find significant subtypes for FGA, WGD, TMB and TMB high
fga_significant <- findSigniciant(fga_pvalues, all_subtypes, q_value_threshold)
wgd_significant <- findSigniciant(wgd_pvalues, all_subtypes, q_value_threshold)
tmb_significant <- findSigniciant(tmb_pvalues, all_subtypes, q_value_threshold)
tmb_high_significant <- findSigniciant(tmb_high_pvalues, all_subtypes, q_value_threshold)

significant_df <- data.frame(SUBTYPE = all_subtypes,
                             FGA_SIGNIFICANT = fga_significant,
                             WGD_SIGNIFICANT = wgd_significant,
                             TMB_SIGNIFICANT = tmb_significant,
                             TMB_HIGH_SIGNIFICANT = tmb_high_significant)

# Save to file
write.csv(significant_df, "Plot_2/Figure2B/significant_subtypes.csv", row.names = FALSE)
