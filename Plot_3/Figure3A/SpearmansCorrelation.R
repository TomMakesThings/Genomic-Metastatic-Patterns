library(readxl)
library(DescTools)
library(ggplot2)
library(rcompanion)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Specify the column names for FGA and TMB in the sample data
measure_columns <- list(fga = "Our_FGA", tmb = "Our_TMB")

# Open calculated TMB and FGA data
samples_data <- read.csv(file = "Plot_2/FGA/calculated_TMB_and_FGA.csv",
                         header = TRUE, fill = TRUE)

# Get columns of interest
samples_data <- samples_data[c("SAMPLE_ID", "SUBTYPE", "SAMPLE_TYPE", "PRIMARY_SITE",
                               "METASTATIC_SITE", "MET_SITE_COUNT",
                               measure_columns[["fga"]], measure_columns[["tmb"]])]

# Open data with plot colours
table_s1a <- read_excel("Tables_S1-4/Table_S1.xlsx", sheet = 1, skip = 2)

# Calculate Spearman's correlation coefficient between FGA/TMB and mutational burden
calculateSpearmans <- function(data, colour = "black") {
  fga_correlation <- SpearmanRho(x = data[[measure_columns[["fga"]]]],
                                 y = data$MET_SITE_COUNT,
                                 conf.level = 0.95)
  tmb_correlation <- SpearmanRho(x = data[[measure_columns[["tmb"]]]],
                                 y = data$MET_SITE_COUNT,
                                 conf.level = 0.95)
  
  # Extract p-values
  fga_pvalue <- cor.test(data[[measure_columns[["fga"]]]],
                         data$MET_SITE_COUNT,
                         method = "spearman",
                         exact = FALSE)$p.value
  tmb_pvalue <- cor.test(data[[measure_columns[["tmb"]]]],
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
  # Set the plot colour
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

# Set default colours for non-significant values
spearmans_df$FGA.COLOUR <- "grey"
spearmans_df$TMB.COLOUR <- "grey"

# Iterate through significant subtypes
for (r in row.names(significant_subtypes_df)) {
  # Set the name of the subtype
  subtype <- spearmans_df[r,]$SUBTYPE
  
  if (subtype != "PanCan") {
    # Find samples for the subtype
    subtype_samples <- samples_data[which(samples_data$SUBTYPE == subtype), ]

    # Repeat for FGA and TMB
    for (genomic_feature in unlist(measure_columns)) {

      # Perform power transformations of FGA / TMB with Tukey's ladder of powers
      tukey_features <- transformTukey(subtype_samples[[genomic_feature]],
                                       plotit = FALSE, quiet = TRUE)
      # Normalise between 0 to 1 by subtracting the minimum and dividing by the maximum
      normalised_features <- tukey_features - min(tukey_features)
      normalised_features <- normalised_features / max(normalised_features)

      # Add the normalised feature to the subtype data
      if (genomic_feature == measure_columns[["fga"]]) {
        norm_column <- "NORM_FGA"
      } else {
        norm_column <- "NORM_TMB"
      }
      subtype_samples[norm_column] <- normalised_features
    
      # Run linear regression of metastatic burden with FGA/TMB as the predictive variable
      # and adjust for sample type
      regression_results <- lm(as.formula(paste("MET_SITE_COUNT ~ SAMPLE_TYPE +", norm_column)),
                               data = subtype_samples)
      # Get the p-value for the FGA/TMB coefficient and check if significant
      regression_coefficients <- data.frame(summary(regression_results)$coefficients)
      coefficient_p_value <- regression_coefficients[norm_column,]$Pr...t..
      regression_significant <- coefficient_p_value < 0.05

      if (genomic_feature == measure_columns[["fga"]]) {
        spearmans_df[r,]$FGA.SIGNIFICANT <- spearmans_df[r,]$FGA.SIGNIFICANT & regression_significant
      } else {
        spearmans_df[r,]$TMB.SIGNIFICANT <- spearmans_df[r,]$TMB.SIGNIFICANT & regression_significant
      }
    }
  }
  
  # Set the row for the subtype
  subtype_row <- spearmans_df[r,]
  
  # Set colours for significant values
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
pdf(file = "Plot_3/Figure3A/Figure_3A.pdf", width = 7, height = 10)
ggplot(spearman_plot_data, aes(x = subtype, y = rho, group = type, 
                               ymin = conf_lower, ymax = conf_upper, # Confidence interval
                               color = factor(row.names(spearman_plot_data), # Set as factor so order is the same
                                              levels = row.names(spearman_plot_data)),
                               shape = type)) +
  geom_hline(yintercept = 0, color = "white", size = 1) + # Set line through y-axis center
  geom_pointrange(position = position_dodge(width = -1)) + # Add space between grouped points
  coord_flip() + # Show spearman's rho on x-axis and subtype on y-axis
  labs(title = "Spearman's correlation", x = "", y = "",
       caption = expression("" %<-% "Lower with increasing metastatic burden      Higher with increasing metastatic burden" %->% "")) +
  theme(plot.title = element_text(size = 14, hjust = 0.5, vjust = -1), # Set title position
        plot.caption.position = "plot", # Set caption position and size
        plot.caption = element_text(hjust = 1, size = 8),
        axis.text.y = element_text(hjust = 1, colour = rev(spearmans_df$COLOUR)), # Colour labels
        legend.position = "bottom",
        panel.grid.major.x = element_blank(), # Hide x-axis background grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dotted", size = 0.8), # Set style of y-axis background grid lines
        panel.grid.minor.y = element_line(linetype = "dotted", size = 0.8)) + 
  scale_y_continuous(position = "right", limits = c(-0.6, 0.6)) + # Set y-axis limits
  scale_color_manual(values = c(spearmans_df$FGA.COLOUR, spearmans_df$TMB.COLOUR), # Colour points
                     guide ='none') +
  scale_shape_manual(name = "", values = c(16, 18)) # Set FGA and TMB point shapes
dev.off()