setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Import ASCETS functions from script
source("ASCETS/ascets_resources.R")

# Open CNA segmentation file
cna_hg19 <- read.table(file = "Data/data_cna_hg19.seg", header = TRUE)

# 
ascets_output <- ascets(cna = cna_hg19, cytoband = , 
                        min_cov = 0.5)
write_outputs_to_file(ascets_output, location = "output_folder/")