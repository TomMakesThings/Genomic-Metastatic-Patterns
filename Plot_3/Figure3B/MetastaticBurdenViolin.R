library(ggplot2)
library(vioplot)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Open calculated TMB and FGA data
samples_data <- read.csv(file = "Plot_2/FGA/calculated_TMB_and_FGA.csv",
                          header = TRUE, fill = TRUE)

# Specify the subtypes to plot and their colours
subtypes <- c(pancancer = "PanCan", prostate = "Prostate Adenocarcinoma", 
              uterine = "Uterine Hypermutated", colorectal = "Colorectal MSS")
subtype_colours <- list(pancancer = "#383838", prostate = "#be1e2d",
                        uterine = "#946d2c", colorectal = "#007eb5")

# Find TMB and FGA for samples with 1 to 6+ metastatic sites
metastaticStats <- function(data, metastastic_counts = 1:6) {
  
  for (count in metastastic_counts) {
    if (count >= metastastic_counts[length(metastastic_counts)]) {
      # Find samples with metastatic burden >= 6
      meta_count_samples <- data[which(data$MET_SITE_COUNT >= count), ]
    } else {
      # Find samples with metastatic burden < 6
      meta_count_samples <- data[which(data$MET_SITE_COUNT == count), ]
    }
    
    # Split into primary and metastatic samples
    primary_meta_samples <- meta_count_samples[which(meta_count_samples$SAMPLE_TYPE == "Primary"),]
    metastasis_meta_samples <- meta_count_samples[which(meta_count_samples$SAMPLE_TYPE == "Metastasis"),]
    
    if (count == metastastic_counts[1]) {
      # Initialise data for first count
      tmb_primary <- data.frame(primary_meta_samples$Our_TMB)
      tmb_metastasis <- data.frame(metastasis_meta_samples$Our_TMB)
      fga_primary <- data.frame(primary_meta_samples$Our_FGA)
      fga_metastasis <- data.frame(metastasis_meta_samples$Our_FGA)
      
    } else {
      # Otherwise add columns to data
      tmb_primary <- merge(tmb_primary, data.frame(primary_meta_samples$Our_TMB),
                           by = "row.names", all = TRUE)[-1]
      tmb_metastasis <- merge(tmb_metastasis, data.frame(metastasis_meta_samples$Our_TMB),
                              by = "row.names", all = TRUE)[-1]
      fga_primary <- merge(fga_primary, data.frame(primary_meta_samples$Our_FGA),
                           by = "row.names", all = TRUE)[-1]
      fga_metastasis <- merge(fga_metastasis, data.frame(metastasis_meta_samples$Our_FGA),
                              by = "row.names", all = TRUE)[-1]
    }
    
    # Reset column names to prevent duplicates caused by merge
    colnames(tmb_primary) <- 1:count
    colnames(tmb_metastasis) <- 1:count
    colnames(fga_primary) <- 1:count
    colnames(fga_metastasis) <- 1:count
  }
  
  return(list(tmb_primary = tmb_primary,
              tmb_metastasis = tmb_metastasis,
              fga_primary = fga_primary,
              fga_metastasis = fga_metastasis))
}

# Get TMB and FGA per subtype, split by metastatic burden and primary/metastasis
subtype_stats <- list(pancancer = metastaticStats(samples_data),
                      prostate = metastaticStats(samples_data[which(samples_data$SUBTYPE == subtypes[2]), ]),
                      uterine = metastaticStats(samples_data[which(samples_data$SUBTYPE == subtypes[3]), ]),
                      colorectal = metastaticStats(samples_data[which(samples_data$SUBTYPE == subtypes[4]), ]))

# Create 3 by two subplots
par(mfrow = c(4,2), mar = c(3, 4, 4, 2), mgp = c(2, 1, 0))

for (subtype in names(subtype_stats)) {

  # Set the plot colours based on subtype
  colour_primary <- "white"
  colour_metastasis <- as.character(subtype_colours[subtype])
  plot_title = as.character(subtypes[subtype])
  
  if (subtype == "colorectal") {
    x_label = "Metastatic Burden"
  } else {
    x_label = ""
  }
  
  # Plot both FGA or TMB
  for (type in c("FGA", "TMB")) {
    
    # Select the type of data
    if (type == "TMB") {
      primary_data <- subtype_stats[[subtype]]$tmb_primary
      metastasis_data <- subtype_stats[[subtype]]$tmb_metastasis
      title_postfix <- "- TMB (mut/Mb)"
      add_legend <- TRUE
    } else {
      primary_data <- subtype_stats[[subtype]]$fga_primary * 100
      metastasis_data <- subtype_stats[[subtype]]$fga_metastasis * 100
      title_postfix <- "- FGA (%)"
      add_legend <- FALSE
    }
    
    # Plot split violin plots
    vioplot(primary_data,
            colMed = "white",
            side = "left",
            col = colour_primary,
            xlab = x_label,
            main = paste(plot_title, title_postfix),
            pchMed = 20)
    vioplot(metastasis_data,
            colMed = "white",
            side = "right",
            col = colour_metastasis,
            pchMed = 20,
            add = TRUE)
    
    # Add legend only to right plots
    if (add_legend) {
      legend("topright",
             legend = c("Primary", "Metastasis"),
             fill = c(colour_primary, colour_metastasis),
             cex = 0.8)
    }
  }
}
