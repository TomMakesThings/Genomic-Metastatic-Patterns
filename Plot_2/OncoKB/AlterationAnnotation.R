library(readr)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Remove NA columns from dataframe
removeEmptyCols <- function(df) {
  return(df[ ,colSums(is.na(df)) < nrow(df)])
}

add_oncotree_codes <- T

if (add_oncotree_codes) {
  # Open samples data
  samples_data <- read.table(file = "Data/data_clinical_sample.txt", sep = '\t', 
                             quote = "", header = TRUE, fill = TRUE)
  # Open mutation data
  mutations <- removeEmptyCols(read.table(file = "Data/data_mutations.txt", 
                                          header = TRUE, fill = TRUE))
  
  # Add placeholder columns for OncoTree code and sample ID
  mutations$OncoTree_Code <- NA
  mutations$Sample_ID <- NA
  
  # Add OncoTree codes for each mutation
  for (r in 1:nrow(mutations)) {
    mutation <- mutations[r,]
    
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
      mutations[r,]$OncoTree_Code <- sample$ONCOTREE_CODE
      mutations[r,]$Sample_ID <- sample_id
    }
  }
  
  write.csv(mutations, file = "Plot_2/OncoKB/mutations_with_oncotree_codes.csv")
}



# Open OncoKB annotated mutation data from zip file
oncokb_results <- removeEmptyCols(data.frame(read_csv(unzip("oncokb_results.zip",
                                                            "oncokb_results.csv"))))

oncokb_results$HIGHEST_LEVEL
