library(readxl)
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

# Open data with plot colours
table_s1a <- read_excel("Tables_S1-4/Table_S1.xlsx", sheet = 1, skip = 2)

# Calculate Spearman's correlation coefficient between FGA/TMB and mutational burden
calculateSpearmans <- function(data, colour = "black") {
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
                               FGA.CI.LOW = unname(fga_correlation["lwr.ci"]), # FGA confidence interval lower
                               FGA.CI.HIGH = unname(fga_correlation["upr.ci"]), # FGA confidence interval upper
                               FGA.PVAL = fga_pvalue, # FGA p-value
                               TMB.RHO = unname(tmb_correlation["rho"]), # TMB spearmans
                               TMB.CI.LOW = unname(tmb_correlation["lwr.ci"]), # TMB confidence interval lower
                               TMB.CI.HIGH = unname(tmb_correlation["upr.ci"]), # TMB confidence interval upper
                               TMB.PVAL = tmb_pvalue,
                               COLOUR = colour) # TMP p-value
  
  return(spearmans_data)
}

# Find the 50 subtype names
cancer_subtypes <- unique(samples_data$SUBTYPE)
metastatic_samples <- samples_data[which(samples_data$MET_SITE_COUNT > 0), ]

# Calculate Spearman's correlation across all samples
spearmans_df <- calculateSpearmans(metastatic_samples)
spearmans_df$SUBTYPE <- "PanCan"

# Add Spearman's correlation for each subtype
for (subtype in cancer_subtypes) {
  # Find samples for the subtype
  subtype_samples <- metastatic_samples[which(metastatic_samples$SUBTYPE == subtype), ]
  # Filter to find those with metastasis
  subtype_samples <- subtype_samples[which(subtype_samples$MET_SITE_COUNT > 0), ]
  
  subtype_colour <- table_s1a[which(table_s1a$curated_subtype_display == subtype), ]$color_subtype
  
  # Add row for spearman's correlation with FGA and TMB
  spearmans_df <- rbind(spearmans_df,
                        calculateSpearmans(subtype_samples, colour = subtype_colour))
}

# Adjust p-values for false discovery rate
spearmans_df$FGA.QVAL <- p.adjust(spearmans_df$FGA.PVAL, method = "fdr")
spearmans_df$TMB.QVAL <- p.adjust(spearmans_df$TMB.PVAL, method = "fdr")

# Mark as significant if q-value < 0.05
spearmans_df$FGA.SIGNIFICANT <- spearmans_df$FGA.QVAL < 0.05
spearmans_df$TMB.SIGNIFICANT <- spearmans_df$TMB.QVAL < 0.05

# Find subtypes that are significant
significant_subtypes_df <- spearmans_df[which(spearmans_df$FGA.SIGNIFICANT |
                                                spearmans_df$TMB.SIGNIFICANT), ]

# Set colours for significant values
spearmans_df$FGA.COLOUR <- "grey"
spearmans_df$TMB.COLOUR <- "grey"

for (r in row.names(significant_subtypes_df)) {
  subtype_row <- spearmans_df[r,]
  
  if (subtype_row$FGA.SIGNIFICANT) {
    spearmans_df[r,]$FGA.COLOUR <- subtype_row$COLOUR
  }
  if (subtype_row$TMB.SIGNIFICANT) {
    spearmans_df[r,]$TMB.COLOUR <- subtype_row$COLOUR
  }
}

# Rearrange data for plotting
spearman_plot_data <- data.frame(subtype = rep(paste(spearmans_df$SUBTYPE,
                                                     " (", spearmans_df$N.SAMPLES,
                                                     ")", sep = ""), 2),
                                 type = c(rep("FGA", nrow(spearmans_df)),
                                          rep("TMB", nrow(spearmans_df))),
                                 rho = c(spearmans_df$FGA.RHO,
                                         spearmans_df$TMB.RHO),
                                 qval = c(spearmans_df$FGA.QVAL,
                                          spearmans_df$TMB.QVAL),
                                 conf_lower = c(spearmans_df$FGA.CI.LOW,
                                                      spearmans_df$TMB.CI.LOW),
                                 conf_upper = c(spearmans_df$FGA.CI.HIGH,
                                                spearmans_df$TMB.CI.HIGH))

# Set subtype as a factor so order remains the same as in the dataframe
spearman_plot_data$subtype <- factor(spearman_plot_data$subtype,
                                     levels = rev(paste(spearmans_df$SUBTYPE,
                                                    " (", spearmans_df$N.SAMPLES,
                                                    ")", sep = "")))

# Plot point range with spearman's rho
ggplot(spearman_plot_data, aes(x = subtype, y = rho, group = type, 
                               ymin = conf_lower, ymax = conf_upper, # Confidence interval
                               color = factor(row.names(spearman_plot_data), # Set as factor so order is the same
                                              levels = row.names(spearman_plot_data)),
                               shape = type)) +
  geom_pointrange(position = position_dodge(width = -1)) + # Add space between grouped points
  coord_flip() +
  labs(title = "Spearman's correlation", x = "", y = "") +
  scale_y_continuous(position = "right", limits = c(-0.6, 0.6)) +
  theme(axis.text.y = element_text(hjust = 1, colour = rev(spearmans_df$COLOUR)), # Colour labels
        legend.position = "bottom") + 
  scale_color_manual(values = c(spearmans_df$FGA.COLOUR, spearmans_df$TMB.COLOUR),
                     guide ='none') + # Colour points
  scale_shape_manual(name = "", values = c(16, 18))
