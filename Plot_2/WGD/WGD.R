library(readxl)
library(facetsSuite)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# View table
table_s2a <- read_excel("Tables_S1-4/Table_S2.xlsx", sheet = 1, skip = 2)

msk_met <- read.table(file = "Plot_2/WGD/msk_met_2021_clinical_data.tsv",
                      sep = '\t', header = TRUE)

primary_samples <- msk_met[which(msk_met$Sample.Type == "Primary"), ]
metastatic_samples <- msk_met[which(msk_met$Sample.Type != "Primary"), ]

message(paste("Number of primary samples: ", nrow(primary_samples)))
message(paste("Number of metastatic samples: ", nrow(metastatic_samples)))

# Find samples with FACETS 
facets_fga_df <- table_s2a[which(table_s2a$alteration == "fga_facets"), ]
sum(facets_fga_df$primary_n) + sum(facets_fga_df$metastasis_n)

# Open mutation data
mutations <- read.table(file = "Plot_2/WGD/data_mutations.txt", 
                        header = TRUE, fill = TRUE)
# Remove empty columns
mutations <- mutations[ ,colSums(is.na(mutations)) < nrow(mutations)]

library(pctGCdata)
test_read_counts

test_facets_output <- run_facets(test_read_counts, cval = 500, genome = 'hg38')
calculate_fraction_cna(test_facets_output$segs, test_facets_output$ploidy, 'hg38', 'em')

arm_level_changes(test_facets_output$segs, test_facets_output$ploidy, 'hg38', 'em')

check_fit(test_facets_output)
