setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Open mutation data
mutations <- read.table(file = "Data/data_mutations.txt", 
                        header = TRUE, fill = TRUE)

# Remove NA columns from dataframe
removeEmptyCols <- function(df) {
  return(df[ ,colSums(is.na(df)) < nrow(df)])
}

mutations <- removeEmptyCols(mutations)

unique(mutations$Tumor_Sample_Barcode)
