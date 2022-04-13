library(readxl)
library(facetsSuite)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# View S2A table
table_s2a <- read_excel("Tables_S1-4/Table_S2.xlsx", sheet = 1, skip = 2)

msk_met <- read.table(file = "Plot_2/WGD/msk_met_2021_clinical_data.tsv",
                      sep = '\t', header = TRUE)

primary_samples <- msk_met[which(msk_met$Sample.Type == "Primary"), ]
metastatic_samples <- msk_met[which(msk_met$Sample.Type != "Primary"), ]

message(paste("Number of primary samples: ", nrow(primary_samples)))
message(paste("Number of metastatic samples: ", nrow(metastatic_samples)))

# Find samples with FACETS 
facets_fga_df <- table_s2a[which(table_s2a$alteration == "fga_facets"), ]
message(paste("Found", sum(facets_fga_df$primary_n) + sum(facets_fga_df$metastasis_n),
              "samples with FACETS data"))

facets_fga_df$tumor_type

table_s2a



# Remove NA columns from dataframe
removeEmptyCols <- function(df) {
  return(df[ ,colSums(is.na(df)) < nrow(df)])
}



# Open mutation data
mutations <- removeEmptyCols(read.table(file = "Data/data_mutations.txt", 
                                        header = TRUE, fill = TRUE))

cna <- read.table(file = "Data/data_cna.txt", 
                        header = TRUE, fill = TRUE)

cna_hg19 <- read.table(file = "Data/data_cna_hg19.seg", 
                       header = TRUE)

sv <- read.table(file = "Data/data_sv.txt", sep = '\t', quote = "",
                       header = TRUE, fill = TRUE)

samples <- read.table(file = "Data/data_clinical_sample.txt", sep = '\t', quote = "",
                 header = TRUE, fill = TRUE)


cna_hg19[1:5, 1:5]
summary(cna[1:5, 1:5])

summary(cna_hg19$seg.mean)
summary(cna_hg19$seg.mean)

"P-0008840-T01-IM5" %in% cna_hg19$ID

unique(cna_hg19$chrom)

cna_hg19[which(cna_hg19$ID == "P-0008840-T01-IM5"), ]

library(pctGCdata)
test_read_counts

test_facets_output <- run_facets(test_read_counts, cval = 500, genome = 'hg38')
calculate_fraction_cna(test_facets_output$segs, test_facets_output$ploidy, 'hg38', 'em')

arm_level_changes(test_facets_output$segs, test_facets_output$ploidy, 'hg38', 'em')

check_fit(test_facets_output)
