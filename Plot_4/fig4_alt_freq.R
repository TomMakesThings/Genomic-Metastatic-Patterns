library(readxl)
library(DescTools)
library(ggplot2)
library(tidyverse)
samples_data <- read.csv(file = "Plot_2/FGA/calculated_TMB_and_FGA.csv",
                         header = TRUE, fill = TRUE)
table_s1a <- read_excel("Tables_S1-4/Table_S1.xlsx", sheet = 1, skip = 2)
table_s2a <- read_excel("Tables_S1-4/Table_S2.xlsx", sheet = 1, skip = 2)

#all cancer types
cancer_subtypes <- unique(table_s4b$tumor_type)
cancer_subtypes
#filter by qvalue
table_s2a_sig_fga<-table_s2a %>% filter(alteration=='fga') %>% filter(qval<0.05)
table_s2a_sig_fga


#tumor_type_sig_fga<-unique(table_s2a_sig_fga$tumor_type)
tumor_type_sig_fga<-c("lung_adenocarcinoma","colorectal_cancer_mss","prostate_adenocarcinoma",
                      "IDC_HR+HER2-","pancreatic_adenocarcinoma", "head_and_neck_squamous",
                      "melanoma_cutaneous","thyroid_papillary","esophageal_adenocarcinoma",
                      "adenoid_cystic_carcinoma","adenoid_cystic_carcinoma","lung_neuroendocrine" ,"mesothelioma_pleural",
                      "sarcoma_ups/mxf"
                      )
subtype_colour=c()
for (subtype in tumor_type_sig_fga) {
    subtype_colour<-c(subtype_colour,table_s1a[which(table_s1a$curated_subtype == subtype), ]$color_subtype)
}
# a dataframe of tumor_type(sig) and corresponding color 
significant_tumor_type=data.frame(tumor_type=tumor_type_sig_fga,colour=subtype_colour)
significant_tumor_type
subtype_colour

#input data for figure 4
table_s4b <- read_excel("Tables_S1-4/Table_S4.xlsx", sheet = 2, skip = 2)
target_organ<-unique(table_s4b$TO)
length(target_organ)
target_organ

#classify alternation
alternation<-unique(table_s4b$alteration)
library(stringr)
alternation_updated=alternation[which(str_detect(alternation,pattern="gain")==FALSE)]
alternation_updated=alternation_updated[which(str_detect(alternation_updated,pattern="loss")==FALSE)]
alternation_updated=alternation_updated[which(str_detect(alternation_updated,pattern="fga")==FALSE)]
alternation_updated=alternation_updated[which(str_detect(alternation_updated,pattern="tmb")==FALSE)]
mut<-alternation_updated[which(str_detect(alternation_updated,pattern="_mut")==TRUE)]
Amplification<-alternation_updated[which(str_detect(alternation_updated,pattern="_Amplification")==TRUE)]
Deletion<-alternation_updated[which(str_detect(alternation_updated,pattern="_Deletion")==TRUE)]
fusion<-alternation_updated[which(str_detect(alternation_updated,pattern="_fusion")==TRUE)]
alternation_updated

pathway=alternation_updated[!alternation_updated %in% c(mut,Amplification,Deletion,fusion)]
non_pathway=c(mut,Amplification,Deletion,fusion)

figure4_plot_info=data.frame()

alt=NA
TO=NA
alt_freq_low=NA
alt_freq_high=NA
bg_color=NA
triangle_color=NA
pval=NA

#assign color to mut,amp,del,fus,pathway
#green,red,blue,purple,navy respectively

map_color=function(s){
  if (s %in% mut){c='green'}
  if (s %in% Amplification){c='red'}
  if (s %in% Deletion){c='blue'}
  if (s %in% fusion){c='purple'}
  if (s %in% pathway){c='navy'}
  return(c)
}


for (i in 1:length(significant_tumor_type$tumor_type)){
  temp=table_s4b[(table_s4b$tumor_type==significant_tumor_type$tumor_type[i]),]
  temp=temp %>%filter(qval<0.05)
  #non_pathway come first(italic font in graph)
  for (p in 1:length(temp$tumor_type)){
    if (((temp$alteration[p] %in% non_pathway)==TRUE) && ((temp$TO[p] %in% target_organ)==TRUE)){
        alt=c(alt,temp$alteration[p])
        TO=c(TO,temp$TO[p])
        alt_freq_low=c(alt_freq_low,temp$noTO_pc[p])
        alt_freq_high=c(alt_freq_high,temp$TO_pc[p])
        bg_color=c(bg_color,significant_tumor_type$colour[i])
        triangle_color=c(triangle_color,map_color(temp$alteration[p]))
        pval=c(pval,temp$pval[p])
    }
    #pathway(regular font, under nonpathway)
  }
  for (p in 1:length(temp$tumor_type)){
    if (((temp$alteration[p] %in% pathway)==TRUE) && ((temp$TO[p] %in% target_organ)==TRUE)){
      alt=c(alt,temp$alteration[p])
      TO=c(TO,temp$TO[p])
      alt_freq_low=c(alt_freq_low,temp$noTO_pc[p])
      alt_freq_high=c(alt_freq_high,temp$TO_pc[p])
      bg_color=c(bg_color,significant_tumor_type$colour[i])
      triangle_color=c(triangle_color,map_color(temp$alteration[p]))
      pval=c(pval,temp$pval[p])
    }
  }
}


figure4_plot_info=data.frame(alternation=alt,TO=TO,alternation_freq1=alt_freq_low,alternation_freq2=alt_freq_high,
                             bg_color=bg_color,triangle_color=triangle_color,pval=pval)
figure4_plot_info
write.csv(figure4_plot_info, file = "Plot_4/figure4_plot_info.csv",row.names = FALSE)

#file sorted using excel, this process is done in excel, so we just open the sorted file directly in the next step
figure4_plot_info_sorted <- read.csv(file = "Plot_4/figure4_plot_info_sorted.csv",
                         header = TRUE, fill = TRUE)


library(randomcoloR)
#randomly generate color for each of the target organ
target_organ_color=randomcoloR::distinctColorPalette(length(target_organ))

TO_color=c()
for (i in 1:nrow(figure4_plot_info_sorted)){
  for (j in 1:length(target_organ)){
    if (figure4_plot_info_sorted$TO[i]==target_organ[j]){
      TO_color=c(TO_color,target_organ_color[j])
    }
  }
}
figure4_plot_info_sorted$TO_color=TO_color
figure4_plot_info_sorted$alt_TO=paste(figure4_plot_info_sorted$alternation,'     ',figure4_plot_info_sorted$TO)
figure4_plot_info_sorted$alt_TO

#update sorted file(with generated color for each TO)
write.csv(figure4_plot_info_sorted, file = "Plot_4/figure4_plot_info_sorted.csv",row.names = FALSE)
## plot for lung adeno
plot_subfigure<-function(colourtype,titlename,r){
  figure4_plot_info_sorted %>% filter(bg_color==colourtype) %>%  
    mutate(plotdata=alt_TO) %>%
    ggplot(aes(x=alternation_freq1*100,y=plotdata,xmin=alternation_freq1,xmax=alternation_freq2,color=factor(triangle_color)))+
    geom_point(aes(x=alternation_freq1*100,y=plotdata),size=4,shape=20)+
    geom_text(aes(label=round(alternation_freq1,2)*100),hjust=-1,vjust=0.5)+
    geom_point(aes(x=alternation_freq2*100,y=plotdata),size=5,shape=20)+
    geom_text(aes(label=round(alternation_freq2,2)*100),hjust=-3,vjust=0.5)+
    theme(panel.background = element_rect(fill = alpha(colourtype,0.5), color =colourtype))+
    scale_x_continuous(position='top',limits = c(0, 100))+
    xlab(label='alternation frequency')+
    ylab(label='alternation')+
    scale_colour_identity()+
    labs(title = titlename)+
    theme(plot.title = element_text(size = 15, color = colourtype)) +
    #dev.new(width = 8, height = 8, unit = "in")
    theme(aspect.ratio = r )
}
plot_subfigure('#b368d9','Lung Adenocarcinoma',2.5)
plot_subfigure('#007eb5','Colorectal MSS',1.5)
plot_subfigure('#be1e2d','Prostate Adenocarcinoma',1.8)

plot_subfigure('#355a91','pancreatic_adenocarcinoma',0.3)
plot_subfigure('#e6308e','Breast Ductal HR+HER2−',0.5)
plot_subfigure('#049347','Head and Neck Squamousa',0.08)

plot_subfigure('#4ebe03','Melanoma Cutaneous',0.08)
plot_subfigure('#cccc33','thyroid_papillary',0.15)
plot_subfigure('#7aa8d9','esophageal_adenocarcinoma',0.08)
plot_subfigure('#026631','adenoid_cystic_carcinoma',0.05)

plot_subfigure('#310725','lung_neuroendocrine',0.05)
plot_subfigure('#682767','mesothelioma_pleural',0.05)
plot_subfigure('#00dbca','sarcoma_ups/mxf',0.05)
