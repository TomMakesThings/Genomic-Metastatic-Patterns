library(readr)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Remove NA columns from dataframe
removeEmptyCols <- function(df) {
  return(df[ ,colSums(is.na(df)) < nrow(df)])
}

# Whether to annotate the mutations MAF file with OncoTree codes
add_oncotree_codes <- F

if (add_oncotree_codes) {
  # Open samples data
  samples_data <- read.table(file = "Data/data_clinical_sample.txt", sep = '\t', 
                             quote = "", header = TRUE, fill = TRUE)
  # Open mutation data
  mutations <- removeEmptyCols(read.table(file = "Data/data_mutations.txt", 
                                          header = TRUE, fill = TRUE))
  
  mutation_samples <- mutations[c("dbSNP_RS", "dbSNP_Val_Status")]
  
  n_mutations <- nrow(mutations)
  
  # Placeholders to record OncoTree code, sample ID and subtype per mutation
  oncotree_codes <- rep(NA, n_mutations)
  sample_ids <- rep(NA, n_mutations)
  subtypes <- rep(NA, n_mutations)
  sample_types <- rep(NA, n_mutations)
  
  # Add OncoTree codes for each mutation
  for (r in 1:n_mutations) {
    mutation <- mutation_samples[r,]
    
    # Find the sample ID for the mutation
    if (substr(mutation$dbSNP_RS[[1]], 1, 1) == "P") {
      sample_id <- mutation$dbSNP_RS[[1]]
    } else if (substr(mutation$dbSNP_Val_Status[[1]], 1, 1) == "P") {
      sample_id <- mutation$dbSNP_Val_Status[[1]]
    } else {
      # One mutation does not have a sample ID
      sample_id <- NA
    }
    
    if (!is.na(sample_id)) {
      # Find and record the OncoTree code and sample ID
      sample <- samples_data[which(samples_data$SAMPLE_ID == sample_id), ]
      oncotree_codes[r] <- sample$ONCOTREE_CODE
      sample_ids[r] <- sample_id
      subtypes[r] <- sample$SUBTYPE
      sample_types[r] <- sample$SAMPLE_TYPE
    }
  }
  
  mutations$OncoTree_Code <- oncotree_codes
  mutations$Sample_ID <- sample_ids
  mutations$Subtype <- subtypes
  mutations$Sample_Type <- sample_types
  
  write.table(mutations, file = "Plot_2/OncoKB/mutations_with_oncotree_codes.txt",
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# oncokb_results <- read.table("Plot_2/OncoKB/oncokb_annotated_mutations.txt",
#                              header = TRUE, fill = TRUE, sep = "\t")

# Open OncoKB annotated mutation data from zip file
# oncokb_results <- removeEmptyCols(data.frame(read.table(unzip("Plot_2/OncoKB/oncokb_annotated_mutations.zip",
#                                                               "oncokb_annotated_mutations.txt"),
#                                                         header = TRUE, fill = TRUE, sep = "\t")))
