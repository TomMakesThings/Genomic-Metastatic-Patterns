# Set the mutations MAF file to a format suitable for OncoKB annotation by adding OncoTree codes

library(readr)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Remove NA columns from dataframe
removeEmptyCols <- function(df) {
  return(df[ ,colSums(is.na(df)) < nrow(df)])
}

# Open samples data
samples_data <- read.table(file = "Data/data_clinical_sample.txt", sep = '\t', 
                           quote = "", header = TRUE, fill = TRUE)
# Open mutation data
mutations <- removeEmptyCols(read.delim(file = "Data/data_mutations.txt", 
                                        header = TRUE, fill = TRUE))

# Get the sample IDs and number of mutations
mutation_samples <- mutations["Tumor_Sample_Barcode"]
n_mutations <- nrow(mutations)

# Placeholders to record OncoTree code, sample ID and subtype per mutation
oncotree_codes <- rep(NA, n_mutations)
sample_ids <- rep(NA, n_mutations)
subtypes <- rep(NA, n_mutations)
sample_types <- rep(NA, n_mutations)

# Add OncoTree codes for each mutation
for (r in 1:n_mutations) {
  # Find the sample ID for the mutation
  sample_id <- mutation_samples[r,]
  
  if (!is.na(sample_id)) {
    # Find and record the OncoTree code and sample ID
    sample <- samples_data[which(samples_data$SAMPLE_ID == sample_id), ]
    oncotree_codes[r] <- sample$ONCOTREE_CODE
    sample_ids[r] <- sample_id
    subtypes[r] <- sample$SUBTYPE
    sample_types[r] <- sample$SAMPLE_TYPE
  }
}

# Set the columns required for MafAnnotator
oncokb_format <- mutations[c("Hugo_Symbol", "Variant_Classification",
                             "Tumor_Sample_Barcode", "HGVSp_Short",
                             "Start_Position", "End_Position")]
oncokb_format$Sample_ID <- sample_ids
oncokb_format$OncoTree_Code <- oncotree_codes
oncokb_format$Subtype <- subtypes
oncokb_format$Sample_Type <- sample_types

# Save to file
write.table(oncokb_format, file = "Plot_2/OncoKB/mutations_with_oncotree_codes.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# View results annotated with OncoKb
oncokb_results <- read.delim("Plot_2/OncoKB/oncokb_annotated_mutations.txt",
                             header = TRUE, fill = TRUE, sep = "\t")
