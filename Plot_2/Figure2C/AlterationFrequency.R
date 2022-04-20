library(readxl)
library(ggplot2)
library(stringr)

setwd("C:/Users/redds/Documents/GitHub/Genomics-II-Group/")

# Open tables S1A and S2A 
table_s1a <- data.frame(read_excel("Tables_S1-4/Table_S1.xlsx", sheet = 1, skip = 2))
table_s2a <- data.frame(read_excel("Tables_S1-4/Table_S2.xlsx", sheet = 1, skip = 2))
# Open calculated TMB and FGA data
samples_data <- read.csv(file = "Plot_2/FGA/calculated_TMB_and_FGA.csv",
                         header = TRUE, fill = TRUE)

# Find rows with oncogenic alterations, i.e. amplifications, mutations, deletions, fusions
oncogenic_alterations <- table_s2a[which(table_s2a$alteration_type == "oncogenic alteration"),]

# Form a dataframe of recurrent alterations 
recurrent_alterations <- data.frame()

for (r in 1:nrow(oncogenic_alterations)) {
  # Get the row for an alteration
  alteration_row <- oncogenic_alterations[r,]
  
  # Record alteration if present in at least 5% of either primary or metastatic samples
  if (alteration_row$primary_pc > 0.05 | alteration_row$metastasis_pc > 0.05) {
    recurrent_alterations <- rbind(recurrent_alterations, alteration_row)
  }
}

# Remove columns of all NA
recurrent_alterations <- recurrent_alterations[, colMeans(is.na(recurrent_alterations)) != 1]

# Find significant alterations with q-value < 0.05
significant_alterations <- recurrent_alterations[which(recurrent_alterations$qval < 0.05),]

# List tumour types to plot
tumour_types <- unique(significant_alterations$tumor_type)

# Create dataframe to plot the alterations per tumour
alteration_plot_df <- data.frame()

for (tumour in tumour_types) {
  # Get the alterations for the tumour
  tumour_alts <- significant_alterations[which(significant_alterations$tumor_type == tumour), ]
  # Find the plot background colour
  tumour_colour <- table_s1a[which(table_s1a$curated_subtype == tumour), ]$color_subtype
  
  for (r in 1:nrow(tumour_alts)) {
    # Get the row for an alteration
    alteration_row <- tumour_alts[r,]
    
    # Get the alteration frequency for primary and metastatic samples
    primary_freq <- alteration_row$primary_pc
    metastasis_freq <- alteration_row$metastasis_pc
    
    if (primary_freq > metastasis_freq) {
      higher_in_metastasis <- FALSE
    } else {
      higher_in_metastasis = TRUE
    }
    
    # Split alteration into name and type
    # E.g. CCND1_FGF19_Amplification becomes "CCND1 FGF19" and "Amplification"
    alteration <- strsplit(alteration_row$alteration, split = '_', fixed = TRUE)
    alt_name <- paste(alteration[[1]][1:(length(alteration[[1]])-1)])
    alt_type <- alteration[[1]][length(alteration[[1]])]
    
    # Set the colour of the triangle to plot based on the alteration type
    if (tolower(alt_type) == "mut") {
      triangle_color <- "green"
    } else if (tolower(alt_type) == "amplification") {
      triangle_color <- "red"
    } else if (tolower(alt_type) == "deletion") {
      triangle_color <- "blue"
    } else if (tolower(alt_type) == "purple") {
    } else {
      triangle_color <- "purple"
    }
    
    # Append a row to the alteration dataframe
    alteration_plot_df <- rbind(alteration_plot_df,
                                data.frame(alteration_name = alt_name,
                                           alteration_type = alt_type,
                                           alternation_freq1 = min(primary_freq, metastasis_freq),
                                           alternation_freq2 = max(primary_freq, metastasis_freq),
                                           higher_in_metastasis = higher_in_metastasis,
                                           bg_color = tumour_colour,
                                           triangle_color = triangle_color))
  }
}

write.csv(alteration_plot_df, file = "Plot_2/Figure2C/figure2c_plot_info.csv", row.names = FALSE)

# tumour_types <- c("lung_adenocarcinoma", "prostate_adenocarcinoma",
#                   "IDC_HR+HER2-", "colorectal_cancer_mss", "ILC_HR+HER2-",
#                   "pancreatic_adenocarcinoma", "pancreatic_neuroendocrine",
#                   "uterine_endometrioid", "gastrointestinal_stromal_tumor",
#                   "head_and_neck_squamous", "uterine_hypermutant",
#                   "lung_squamous_cell_carcinoma", "renal_clear_cell",
#                   "melanoma_cutaneous", "thyroid_papillary",
#                   "esophageal_adenocarcinoma")

# for (tumour_type in all_tumour_types) {
#   # Get samples for the subtype
#   tumour_rows <- recurrent_alterations[which(recurrent_alterations$tumor_type == tumour_type), ]
# 
#   # Mann-Whitney U test to compare alteration frequencies between primary and metastatic samples
#   man_whitney_results <- wilcox.test(c(0.06153846, 0.0438))
#   tumour_pvalues <- c(tumour_pvalues, man_whitney_results$p.value)
# }
# 
# # Convert from p-values to q-values and return which subtypes are significant
# findSigniciant <- function(pvalues, subtype_list, threshold) {
#   qvalues <- p.adjust(tmb_high_pvalues, method = "fdr")
#   
#   return(qvalues < threshold)
# }
