# Run ASCETS, the arm-level copy number events caller for targeted sequencing data

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Import ASCETS functions from script
source("ASCETS/ascets_resources.R")

# Open CNA segmentation file
cna_hg19 <- read.table(file = "Data/data_cna_hg19.seg", header = TRUE)
# Open chromosome arm genomic coordinates for reference hg19
cytoband_hg19 <- read.table(file = "ASCETS/genomic_arm_coordinates_hg19.txt", header = TRUE)

# Run ASCETS
ascets_output <- ascets(cna = cna_hg19, cytoband = cytoband_hg19, name = "ascets_results")
write_outputs_to_file(ascets_output, location = "Plot_2/aSCNAs")
