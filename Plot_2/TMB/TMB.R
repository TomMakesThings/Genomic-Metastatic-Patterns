library(dplyr)
library(reshape2)

setwd("~/AllMphil/GenomicsAS3/")

## Tumor mutational burden (TMB) was calculated for each sample as the total number of nonsynonymous mutations, divided by the
## number of bases sequenced. 

##IMPACT468 - 1156473 base

#IMPACT341 - 901587 base

#IMPACT410 - 1021743 base

table_S1 <- read.csv("supplementary_table/table_S1.csv", header = TRUE, skip = 2) #25775 obs. of  56 variables
clinical_sample <- read.table("data/data_clinical_sample.txt", sep = "\t",header = TRUE)
mutations <- read.table("data/data_mutations.txt", sep = "\t",quote = "",header = TRUE)
#sv <- read.table("data/data_sv.txt", sep = "\t",quote = "",header = TRUE) #Structural Variant
#seg19 <- read.table("data/data_cna_hg19.seg", sep = "\t", header = TRUE)
panels <- read.table("data/data_gene_panel_matrix.txt", sep = "\t",header = TRUE)
cbioportal <- read.table("supplementary_table/cbioportal.tsv", sep = "\t", header = TRUE)


true_TMB <- data.frame(sampleID = cbioportal$Sample.ID, TMB = cbioportal$TMB..nonsynonymous.)
mut <- data.frame(sampleID = mutations$Tumor_Sample_Barcode, consequence = mutations$Consequence,
                  variant_classification = mutations$Variant_Classification,
                  variant_type = mutations$Variant_Type, Is_nonsynonymous=NA) #some of the sampleID do not have mutation data, in those cases, TMB=0



## set the total length (base) of sequence of each platform 
platforms <- unique(panels$mutations) #"IMPACT468" "IMPACT410" "IMPACT341"
sequence_length <- c(1156473, 1021743, 901587)
names(sequence_length) <- platforms


## set the variant_classification that count as "nonsynonymous mutations", in total 10 types

variant_classification <- unique(mut$variant_classification)

no_variant_class <- c("5'Flank", "Translation_Start_Site", "3'Flank",
                      "Silent", "5'UTR", "IGR","Intron") #does not count as "nonsynonymous mutations", 6 types
yes_variant_class <- variant_classification[!variant_classification %in% no_variant_class] # count as "nonsynonymous mutations", 10types


## count the number of "nonsynonymous mutations" of each sample
add.nonsynonymous.mut <- mutate(mut, Is_nonsynonymous = as.numeric(variant_classification %in% yes_variant_class)) #230419 obs. of  5 variables, 24755 unique sampleIDs

merged.mut <- data.frame(sampleID = true_TMB$sampleID, 
                        num_nonsynonymous = 0,
                        panel = NA,
                        total_sequence = NA, #25775 obs. of  5 variables
                        TMB = 0)

#add in panel info of each sample
map_mut_panel <- panels$mutations
names(map_mut_panel) <- panels$SAMPLE_ID

merged.mut$panel = unname(map_mut_panel[merged.mut$sampleID])

#add in sequence length info of each sample
merged.mut <- merged.mut %>%
    mutate(., total_sequence = unname(sequence_length[panel]))


#add in count of nonsynonymous mutations for each sample

count_name <- unique(add.nonsynonymous.mut$sampleID) #24755
count <- rep(0, length(count_name)) #24755

for (i in 1:length(count_name)){
    x <- subset(add.nonsynonymous.mut, sampleID == count_name[i])
    count[i] <- sum(x$Is_nonsynonymous)
} #1:29 starts 1:30 ends

names(count) <- count_name #24755

merged.mut <- mutate(merged.mut, num_nonsynonymous = count[merged.mut$sampleID]) %>%
    mutate_all(~replace(., is.na(.), 0))%>%
    mutate(TMB = num_nonsynonymous * 10^6 / total_sequence)

names(merged.mut$num_nonsynonymous) <- NULL
names(merged.mut$TMB) <- NULL


## output
our_TMB <- cbind.data.frame(clinical_sample[, 1:18], Our_TMB = merged.mut$TMB, 
                            Sup_TMB = table_S1$tmb,
                            clinical_sample[, 19:42])


write.csv(our_TMB, file = "tmb/calculated_TMB.csv",row.names = FALSE)



