library(stringr)
library(reshape)
library(brglm)

# Open calculated TMB and FGA data
samples_data <- read.csv(file = "Plot_2/FGA/calculated_TMB_and_FGA.csv",
                         header = TRUE, fill = TRUE)
# Open ASCETS arm-level CNA results per sample
subtype_mean_cna <- read.csv(file = "Plot_2/aSCNAs/subtype_mean_arm_level_cna.csv",
                             header = TRUE, fill = TRUE)
sample_cna <- read.csv(file = "Plot_2/aSCNAs/sample_arm_level_cna.csv",
                       header = TRUE, fill = TRUE)
arm_level_calls <- read.csv(file = "Plot_2/aSCNAs/sample_arm_level_calls.csv",
                            header = TRUE, fill = TRUE)

# Add FGA to CNA data
sample_cna$FGA <- samples_data$Our_FGA

# Use regex to find the column names for the chromosome arms, i.e. X10p, X10q, X11p...
chromosome_cols <- names(subtype_mean_cna)[str_extract(names(subtype_mean_cna),
                                                       "X\\d?\\d?p?q?") %in% names(subtype_mean_cna)]
# Get subtypes
all_subtypes <- unique(subtype_mean_cna$SUBTYPE)

# Create dataframe of amplification and deletion fractions per subtype
subtype_fraction_alteration <- data.frame()

for (subtype in all_subtypes) {
  # Find primary and metastasis samples for each subtype
  subtype_call_rows <- arm_level_calls[which(arm_level_calls$SUBTYPE == subtype), ]
  # Find primary and metastasis samples for each subtype
  primary_subtype <- subtype_call_rows[which(subtype_call_rows$SAMPLE_TYPE == "Primary"), ]
  metastasis_subtype <- subtype_call_rows[which(subtype_call_rows$SAMPLE_TYPE == "Metastasis"), ]
  
  # Calculate gain and loss frequency per chromosome arm
  alt_frequencies <- list(Gain = list(Primary = list(), Metastasis = list()),
                          Loss = list(Primary = list(), Metastasis = list()))
  
  for (arm_id in chromosome_cols) {
    # Calculate gain and loss frequencies for primary samples
    primary_alt_frequencies <- table(primary_subtype[arm_id])
    primary_amp_frequency <- primary_alt_frequencies["AMP"] / sum(primary_alt_frequencies) * 100
    primary_del_frequency <- primary_alt_frequencies["DEL"] / sum(primary_alt_frequencies) * 100
    
    # Calculate frequencies for metastatic samples
    metastasis_alt_frequencies <- table(metastasis_subtype[arm_id])
    metastasis_amp_frequency <- metastasis_alt_frequencies["AMP"] / sum(metastasis_alt_frequencies) * 100
    metastasis_del_frequency <- metastasis_alt_frequencies["DEL"] / sum(metastasis_alt_frequencies) * 100
    
    # Record the frequencies
    alt_frequencies[["Gain"]][["Primary"]][arm_id] <- primary_amp_frequency
    alt_frequencies[["Gain"]][["Metastasis"]][arm_id] <- metastasis_amp_frequency
    alt_frequencies[["Loss"]][["Primary"]][arm_id] <- primary_del_frequency
    alt_frequencies[["Loss"]][["Metastasis"]][arm_id] <- metastasis_del_frequency
  }
  
  # Iterate through alteration type, i.e. gain and loss
  for (alt_type in names(alt_frequencies)) {
    # Iterate through primary and metastasis sample types
    for (sample_type in names(alt_frequencies[[alt_type]])) {
      # Form row of the subtype frequencies for each sample type and alteration
      fraction_alt_row <- data.frame(alt_frequencies[[alt_type]][[sample_type]])
      # Replace any NA frequencies with 0
      fraction_alt_row[is.na(fraction_alt_row)] <- 0
      
      # Add subtype annotation columns
      fraction_alt_row <- cbind(data.frame(SUBTYPE = subtype_call_rows[1,]$SUBTYPE),
                                data.frame(SAMPLE_TYPE = sample_type, ALTERATION = alt_type),
                                fraction_alt_row)
      # Record in main dataframe
      subtype_fraction_alteration <- rbind(subtype_fraction_alteration, fraction_alt_row)
    }
  }
}

# Save frequency of arm-level copy number alterations between primary tumors and metastases
write.csv(subtype_fraction_alteration, "Plot_2/Figure2B/subtype_arm_level_fraction_alteration.csv",
          row.names = FALSE)

# Normalise a dataframe column between 0 - 1
normaliseColumn <-function(col_values) {
  x <- col_values[!is.na(col_values)]
  x <- (x - min(x)) / (max(x) - min(x))
  col_values[!is.na(col_values)] <- x
  
  return(col_values)
}

# Split by primary and metastatic samples
primary_fractions <- subtype_fraction_alteration[which(subtype_fraction_alteration$SAMPLE_TYPE == "Primary"), ]
metastasis_fractions <- subtype_fraction_alteration[which(subtype_fraction_alteration$SAMPLE_TYPE == "Metastasis"), ]

# Melt data so each gain/loss fraction is on a new row
primary_alt_frac_melt <- melt(primary_fractions, id = setdiff(names(primary_fractions), chromosome_cols))
names(primary_alt_frac_melt)[names(primary_alt_frac_melt) == "variable"] <- "ARM_ID"
names(primary_alt_frac_melt)[names(primary_alt_frac_melt) == "value"] <- "PRIMARY_FRACTION"

metastasis_fractions <- melt(metastasis_fractions, id = setdiff(names(metastasis_fractions), chromosome_cols))
names(metastasis_fractions)[names(metastasis_fractions) == "variable"] <- "ARM_ID"
names(metastasis_fractions)[names(metastasis_fractions) == "value"] <- "METASTASIS_FRACTION"

# Combine melted data
alt_frac_melt <- primary_alt_frac_melt[c("SUBTYPE", "ARM_ID", "PRIMARY_FRACTION")]
alt_frac_melt$PRIMARY_FRACTION <- primary_alt_frac_melt$PRIMARY_FRACTION / 100
alt_frac_melt$METASTASIS_FRACTION <- metastasis_fractions$METASTASIS_FRACTION / 100

lung_frac <- alt_frac_melt[which(alt_frac_melt$SUBTYPE == "Lung Adenocarcinoma"), ]

# Fit the logistic regression model
fitted_model <- brglm(METASTASIS_FRACTION ~ PRIMARY_FRACTION, 
                      family = binomial(), data = lung_frac)
regression_coefficients <- data.frame(coef(summary(fitted_model)))
regression_coefficients["PRIMARY_FRACTION",]$Pr...z.. < 0.05

predict(fitted_model, type = "response") < 0.05
