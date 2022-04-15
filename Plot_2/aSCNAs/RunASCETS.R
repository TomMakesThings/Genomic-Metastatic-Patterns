setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Run ASCETS, the arm-level copy number events caller for targeted sequencing data
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
# Results of anuploidy scores
aneuploidy_score <- read.table(file = "Plot_2/aSCNAs/ascets_results_aneuploidy_scores.txt",
                               header = TRUE, check.names = FALSE)
# Results of arm-level calls
arm_level_calls <- read.table(file = "Plot_2/aSCNAs/ascets_results_arm_level_calls.txt",
                              header = TRUE, check.names = FALSE)
# Clinical sample data
samples_data <- read.table(file = "Data/data_clinical_sample.txt", sep = '\t', 
                           quote = "", header = TRUE, fill = TRUE)

# Get columns of interest from clinical data
samples_data <- samples_data[c("SAMPLE_ID", "SUBTYPE", "SAMPLE_TYPE", 
                               "PRIMARY_SITE", "METASTATIC_SITE", "MET_SITE_COUNT")]

# Select autosomal chromosomes from arm-level calls
autosomal_calls <- arm_level_calls[, !names(arm_level_calls) %in% c("Xp", "Xq")]

# Calculate frequency of chromosome amplification per sample
amp_frequency <- c()

for (sample in autosomal_calls$sample) {
  # Get the row for the sample
  row_calls <- autosomal_calls[which(autosomal_calls$sample == sample), ]
  # Extract calls as a vector
  sample_calls <- as.character(unname(as.vector(row_calls))[2:ncol(row_calls)])
  
  # Get the frequency of each type of call, i.e. AMP, DEL, NC, NEUTRAL
  call_frequency <- table(sample_calls)
  sample_amp <- unname(call_frequency["AMP"] / sum(call_frequency))
  
  # Check if frequency is undefined
  if (is.na(sample_amp)) {
    if ("NC" %in% names(call_frequency)) {
      if (call_frequency["NC"] > sum(call_frequency)/2) {
        # Set as NA if over half the calls are unknown
        amp_frequency <- c(amp_frequency, NA)
      } else {
        # Otherwise set as zero as no amplification found
        amp_frequency <- c(amp_frequency, 0)
      }
    } else {
      amp_frequency <- c(amp_frequency, 0)
    }
  } else {
    # Record the calculated amplification frequency
    amp_frequency <- c(amp_frequency, sample_amp)
  }
}

# Combine data to add labels to CNA results
arm_level_df <- samples_data
arm_level_df$ANEUPLOIDY_SCORE <- aneuploidy_score$aneuploidy_score
# Record estimated WGD
arm_level_df$ARM_WGD <- amp_frequency
arm_level_df <- merge(arm_level_df, average_segmeans,
                      by.x = "SAMPLE_ID", by.y = "sample")

# Get list of subtypes
all_subtypes <- unique(arm_level_df$SUBTYPE)
subtype_cna_df <- data.frame()

for (subtype in all_subtypes) {
  # Find primary and metastasis samples for each subtype
  subtype_rows <- arm_level_df[which(arm_level_df$SUBTYPE == subtype), ]
  primary_subtype <- subtype_rows[which(subtype_rows$SAMPLE_TYPE == "Primary"), ]
  metastasis_subtype <- subtype_rows[which(subtype_rows$SAMPLE_TYPE == "Metastasis"), ]
  
  # Calculate average for all numeric columns
  primary_means <- colMeans(primary_subtype[sapply(primary_subtype, is.numeric)])
  metastasis_means <- colMeans(metastasis_subtype[sapply(metastasis_subtype, is.numeric)])
  
  # Add subtype and sample type to beginning
  primary_means_df <- cbind(data.frame(SUBTYPE = subtype, SAMPLE_TYPE = "Primary"),
                            data.frame(t(primary_means)))
  metastasis_means_df <- cbind(data.frame(SUBTYPE = subtype, SAMPLE_TYPE = "Metastasis"),
                            data.frame(t(metastasis_means)))
  
  # Add averages to subtype CNA dataframe
  subtype_cna_df <- rbind(subtype_cna_df, primary_means_df)
  subtype_cna_df <- rbind(subtype_cna_df, metastasis_means_df)
}

# Save to file
write.csv(arm_level_df, "Plot_2/aSCNAs/sample_arm_level_cna.csv", row.names = FALSE)
write.csv(subtype_cna_df, "Plot_2/aSCNAs/subtype_arm_level_cna.csv", row.names = FALSE)
