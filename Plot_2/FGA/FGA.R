library(dplyr)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Open CNA segmentation file
cna_hg19 <- read.table(file = "Data/data_cna_hg19.seg", header = TRUE)
# Open calculated TMB data
samples_data <- read.csv(file = "Plot_2/TMB/calculated_TMB.csv",
                         header = TRUE, fill = TRUE)

# Calculate fga for each tumour sample
fga <- c()

for (id in unique(cna_hg19$ID)) {
  # Find the CNAs for the sample
  sample_cna <- cna_hg19 %>% filter(ID == id)
  # Calculate the total size of the sequenced genome
  genome_size <- sum(sample_cna$loc.end - sample_cna$loc.start)
  # Determine which rows have absolute log2 copy ratios > 0.2
  log2_thresholded <- as.integer(abs(sample_cna$seg.mean) > 0.2)
  # Scale by the size of each region of the genome
  log2_scaled <- log2_thresholded * (sample_cna$loc.end - sample_cna$loc.start)
  # Divide by genome size to calculate fraction of genome altered (FGA)
  fga <- c(fga, sum(log2_scaled) / genome_size)
}

# Add IDs to calculated FGAs
fga_df <- data.frame(SAMPLE_ID = unique(cna_hg19$ID), Our_FGA = fga)

# Combine with calculated TMB data
fga_samples_data <- merge(samples_data[ , !(names(samples_data) %in% "FGA")],
                          fga_df, by = "SAMPLE_ID")
fga_samples_data$Sup_FGA <- samples_data$FGA

# Save to file
write.csv(fga_samples_data, "Plot_2/FGA/calculated_TMB_and_FGA.csv", row.names = FALSE)
