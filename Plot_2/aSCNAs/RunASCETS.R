# Run ASCETS, the arm-level copy number events caller for targeted sequencing data

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

run_ascets <- F

if (run_ascets) {
  # Import ASCETS functions from script
  source("ASCETS/ascets_resources.R")
  
  # Open CNA segmentation file
  cna_hg19 <- read.table(file = "Data/data_cna_hg19.seg", header = TRUE)
  # Open chromosome arm genomic coordinates for reference hg19
  cytoband_hg19 <- read.table(file = "ASCETS/genomic_arm_coordinates_hg19.txt", header = TRUE)
  
  # Run ASCETS
  ascets_output <- ascets(cna = cna_hg19, cytoband = cytoband_hg19, name = "ascets_results")
  write_outputs_to_file(ascets_output, location = "Plot_2/aSCNAs")
}

# Open results of weighted average segment mean values for each arm in each sample
average_segmeans <- read.table(file = "Plot_2/aSCNAs/ascets_results_arm_weighted_average_segmeans.txt",
                               header = TRUE, check.names = FALSE)
average_segmeans$sample

# Open clinical sample data
samples_data <- read.table(file = "Data/data_clinical_sample.txt", sep = '\t', 
                           quote = "", header = TRUE, fill = TRUE)

# Get columns of interest
samples_data <- samples_data[c("SAMPLE_ID", "SUBTYPE", "SAMPLE_TYPE", 
                               "PRIMARY_SITE", "METASTATIC_SITE", "MET_SITE_COUNT")]

# Combine data to add labels to arm-level CNA results
arm_level_df <- merge(samples_data, average_segmeans,
                      by.x = "SAMPLE_ID", by.y = "sample")

# Save to file
write.csv(arm_level_df, "Plot_2/aSCNAs/arm_level_cna.csv", row.names = FALSE)
