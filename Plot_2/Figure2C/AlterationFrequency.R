library(readxl)
library(ggplot2)
library(stringr)
library(ggpubr)

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
  tumour_display_name <- table_s1a[which(table_s1a$curated_subtype == tumour), ]$curated_subtype_display
  
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
                                data.frame(tumour_type = tumour_display_name,
                                           primary_n = tumour_alts$primary_n[1],
                                           alteration_name = alt_name,
                                           alteration_type = alt_type,
                                           alternation_freq1 = min(primary_freq, metastasis_freq),
                                           alternation_freq2 = max(primary_freq, metastasis_freq),
                                           higher_in_metastasis = higher_in_metastasis,
                                           bg_color = tumour_colour,
                                           triangle_color = triangle_color))
  }
}

# First reorder rows from low to high by lowest alteration frequency
alteration_plot_df <- alteration_plot_df[order(alteration_plot_df$alternation_freq1),]

# Then reorder rows by number of alterations per tumour
alteration_plot_df <- transform(alteration_plot_df, freq = ave(seq(nrow(alteration_plot_df)),
                                                               tumour_type, FUN = length))
alteration_plot_df <- alteration_plot_df[order(-alteration_plot_df$freq), ]
alteration_plot_df$freq <- NULL

# Capitalise first letter of each alteration type, e.g. "deletion" -> "Deletion"
alteration_plot_df$alteration_type <- str_to_title(alteration_plot_df$alteration_type)
# Replace "Mut" with "Mutation"
alteration_plot_df$alteration_type <- str_replace_all(alteration_plot_df$alteration_type,
                                                      "Mut", "Mutation")

# Save plot data to file
write.csv(alteration_plot_df, file = "Plot_2/Figure2C/figure2c_plot_info.csv", row.names = FALSE)

# Set list of formal tumour names
all_tumour_names <- unique(alteration_plot_df$tumour_type)

# Record each plot in a list
tumour_alt_plots <- list()

for (i in 1:length(all_tumour_names)) {
  # Get the name as an index
  tumour_name <- all_tumour_names[i]
  
  # Find the rows for the tumour type
  tumour_alts <- alteration_plot_df[which(alteration_plot_df$tumour_type == tumour_name),]
  # Convert columns to factors to keep order
  tumour_alts$triangle_color <- factor(tumour_alts$triangle_color,
                                          levels = unique(tumour_alts$triangle_color))
  tumour_alts$alteration_name  <- factor(tumour_alts$alteration_name,
                                         levels = unique(tumour_alts$alteration_name))
  tumour_alts$alteration_type  <- factor(tumour_alts$alteration_type,
                                         levels = unique(tumour_alts$alteration_type))
  tumour_alts$higher_in_metastasis  <- factor(tumour_alts$higher_in_metastasis,
                                         levels = unique(tumour_alts$higher_in_metastasis))
  
  if (tumour_alts$higher_in_metastasis[1]) {
    # Set right and left facing triangles
    triangle_shapes <- c("\u25C4", "\u25BA")
    triangle_direction <- c("Higher in\nmetastasis", "Lower in\nMetastasis")
  } else {
    # Set left and right facing triangles
    triangle_shapes <- c("\u25BA", "\u25C4")
    triangle_direction <- c("Higher in metastasis", "Lower in Metastasis")
  }
  
  # Create the plot of alteration frequency for the tumour type
  alt_plot <- ggplot(data = tumour_alts,
                     aes(x = alternation_freq1 * 100,
                         y = alteration_name,
                         xmin = alternation_freq1, xmax = alternation_freq2,
                         color = alteration_type, shape = higher_in_metastasis)) +
    geom_point(size = 6) +
    geom_text(aes(label = round(alternation_freq1, 2) * 100), # Add text to left and right of arrows
              color = "black", size = 3, hjust = 3, vjust = 0.5) +
    geom_text(aes(label = round(alternation_freq2, 2) * 100), 
              color = "black", size = 3, hjust = -3, vjust = 0.5) +
    scale_shape_manual(values = triangle_shapes, # Set triangle shapes
                       labels = triangle_direction) +
    scale_color_manual(values = as.vector(unique(tumour_alts$triangle_color))) +
    scale_x_continuous(position = 'top', limits = c(0, 100),
                       expand = c(0.1, 0.1), labels = function(x) paste0(x, "%")) +
    labs(title = paste(tumour_alts$tumour_type[1], " (", tumour_alts$primary_n[1], ")", sep = ""), 
         x = "Alternation frequency", y = "") +
    theme(plot.title = element_text(size = 12, color = tumour_alts$bg_color),
          panel.background = element_rect(fill = alpha(tumour_alts$bg_color, 0.5), # Set background colour
                                          color = tumour_alts$bg_color),
          panel.grid.major.x = element_blank(), # Hide x-axis background grid lines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dotted", size = 0.7), # Set style of y-axis background grid lines
          panel.grid.minor.y = element_line(linetype = "dotted", size = 0.7)) +
    guides(color = guide_legend(title = ""),
           shape = guide_legend(title = ""))
  alt_plot
  tumour_alt_plots[[tumour_name]] <- alt_plot
}

cairo_pdf("Plot_2/Figure2C/Figure2C.pdf", width = 14, height = 12)
ggarrange(plotlist = tumour_alt_plots, heights = c(4,2.5,1.5,1.5), 
          common.legend = TRUE, legend = "bottom")
graphics.off()
