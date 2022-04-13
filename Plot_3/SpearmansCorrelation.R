library(DescTools)
library(ggplot2)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Open clinical sample data
samples_data <- read.table(file = "Data/data_clinical_sample.txt", sep = '\t', 
                      quote = "", header = TRUE, fill = TRUE)

# Get columns of interest
samples_data <- samples_data[c("SAMPLE_ID", "SUBTYPE", "PRIMARY_SITE",
                               "METASTATIC_SITE", "MET_SITE_COUNT", "FGA", 
                               "TMB_NONSYNONYMOUS", "TUMOR_PURITY")]

message(paste("Found", nrow(samples_data), "samples"))


# Calculate Spearman's correlation coefficient between FGA/TMB and mutational burden
calculateSpearmans <- function(data) {
  fga_correlation <- SpearmanRho(x = data$FGA,
                                 y = data$MET_SITE_COUNT,
                                 conf.level = 0.95)
  tmb_correlation <- SpearmanRho(x = data$TMB_NONSYNONYMOUS,
                                 y = data$MET_SITE_COUNT,
                                 conf.level = 0.95)
  # Extract p-values
  fga_pvalue <- cor.test(data$FGA,
                         data$MET_SITE_COUNT,
                         method = "spearman",
                         exact = FALSE)$p.value
  tmb_pvalue <- cor.test(data$TMB_NONSYNONYMOUS,
                         data$MET_SITE_COUNT,
                         method = "spearman",
                         exact = FALSE)$p.value

  # Add results to dataframe
  spearmans_data <- data.frame(SUBTYPE = data$SUBTYPE[1], # Subtype name
                               N.SAMPLES = nrow(data), # Subtype size
                               FGA.RHO = unname(fga_correlation["rho"]), # FGA spearmans
                               FGA.CI.LOW = unname(fga_correlation["rho"]), # FGA confidence interval lower
                               FGA.CI.HIGH = unname(fga_correlation["rho"]), # FGA confidence interval upper
                               FGA.PVAL = fga_pvalue, # FGA p-value
                               TMB.RHO = unname(fga_correlation["rho"]), # TMB spearmans
                               TMB.CI.LOW = unname(fga_correlation["rho"]), # TMB confidence interval lower
                               TMB.CI.HIGH = unname(fga_correlation["rho"]), # TMB confidence interval upper
                               TMB.PVAL = tmb_pvalue) # TMP p-value
  
  return(spearmans_data)
}

# Find the 50 subtype names
cancer_subtypes <- unique(samples_data$SUBTYPE)

# Calculate Spearman's correlation for all cancer samples
spearmans_df <- calculateSpearmans(samples_data)
spearmans_df$SUBTYPE <- "PanCan"

# Add Spearman's correlation for each subtype
for (subtype in cancer_subtypes) {
  # Find samples for the subtype
  subtype_samples <- samples_data[which(samples_data$SUBTYPE == subtype), ]
  # Filter to find those with metastasis
  subtype_samples <- subtype_samples[which(subtype_samples$MET_SITE_COUNT > 0), ]
  # Add row for spearman's correlation with FGA and TMB
  spearmans_df <- rbind(spearmans_df, calculateSpearmans(subtype_samples))
}

# Adjust p-values for false discovery rate
spearmans_df$FGA.QVAL <- p.adjust(spearmans_df$FGA.PVAL, method = "fdr")
spearmans_df$TMB.QVAL <- p.adjust(spearmans_df$TMB.PVAL, method = "fdr")

# Mark as significant if q-value < 0.05
spearmans_df$FGA.SIGNIFICANT <- spearmans_df$FGA.QVAL < 0.05
spearmans_df$TMB.SIGNIFICANT <- spearmans_df$FGA.QVAL < 0.05


# library(reshape)
# mdata <- melt(spearmans_df, id=c("id","time"))
# 
# 
# ggplot(spearmans_df, aes(x = dose, y = len, group = supp, color = supp)) +
#   coord_flip()
# 
# # Use geom_pointrange
# ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
#   geom_pointrange(aes(ymin=len-sd, ymax=len+sd))