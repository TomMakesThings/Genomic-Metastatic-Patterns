#######   figure2_A.R written by Zhan Shi(zs384)

library("readxl")
library(tidyverse)
tableS2A <- read_excel("Tables_S1-4/Table_S2.xlsx",skip = 1)
tableS2A<-data.frame(tableS2A)

tableS1A<-read_excel("Tables_S1-4/Table_S1.xlsx",skip = 1)
tableS1A<-data.frame(tableS1A)

#our calculated fga and tmb
calculated_fga_tmb<-read.csv(file = "Plot_2/FGA/calculated_TMB_and_FGA.csv",
                             header = TRUE, fill = TRUE)

subtypes=unique(calculated_fga_tmb$SUBTYPE_ABBREVIATION)
subtypes_fullname=unique(calculated_fga_tmb$SUBTYPE)
# get median of fga 
get_median_fga<-function(subtype,sample_type){
  if (sample_type==0){
  temp=calculated_fga_tmb %>% filter(SUBTYPE_ABBREVIATION==subtype) %>% filter(SAMPLE_TYPE=='Primary')
  }
  if (sample_type==1){
  temp=calculated_fga_tmb %>% filter(SUBTYPE_ABBREVIATION==subtype) %>% filter(SAMPLE_TYPE=='Metastasis')
  }
  median=median(temp$Our_FGA)
}
fga_p_median<-c()
fga_m_median<-c()
for (i in 1:length(subtypes)){
  fga_p_median<-c(fga_p_median,get_median_fga(subtypes[i],0))
  fga_m_median<-c(fga_m_median,get_median_fga(subtypes[i],1))
}
fga_data<-data.frame(subtypes_fullname=subtypes_fullname,subtypes=subtypes,p_median=fga_p_median,m_median=fga_m_median)

# get median of tmb
get_median_tmb<-function(subtype,sample_type){
  if (sample_type==0){
    temp=calculated_fga_tmb %>% filter(SUBTYPE_ABBREVIATION==subtype) %>% filter(SAMPLE_TYPE=='Primary')
  }
  if (sample_type==1){
    temp=calculated_fga_tmb %>% filter(SUBTYPE_ABBREVIATION==subtype) %>% filter(SAMPLE_TYPE=='Metastasis')
  }
  median=median(temp$Our_TMB)
}
tmb_p_median<-c()
tmb_m_median<-c()
for (i in 1:length(subtypes)){
  tmb_p_median<-c(tmb_p_median,get_median_tmb(subtypes[i],0))
  tmb_m_median<-c(tmb_m_median,get_median_(subtypes[i],1))
}

tmb_data<-data.frame(subtypes_fullname=subtypes_fullname,subtypes=subtypes,p_median=tmb_p_median,m_median=tmb_m_median)





#### assign color
data2<-tableS2A %>% filter(alteration=='fga')
data2$colour=NA
assign_color=function(data){
  for (i in 1:length(data[,1])){
    for (j in 1:length(tableS1A[,1])){
      
    
      if (data[i,1]==tableS1A[j,1]){
        data$colour[i]=tableS1A$color_subtype[j]
      }
    }
  }
  return(data)
}

assign_col=function(data){
  for (i in 1:length(data[,1])){
    for (j in 1:length(tableS1A[,1])){
      
      
      if (data[i,1]==tableS1A[j,9]){
        data$colour[i]=tableS1A$color_subtype[j]
      }
    }
  }
  return(data)
}

data2<-assign_color(data2)
fga_data<-assign_col(fga_data)

#threshold used for figure2A
data2_sig=data2 %>% filter(qval<0.05)
data2_sig
data2

#plot figure2A_fga origin
data2$tumor_type=factor(data2$tumor_type,levels=data2$tumor_type)
data2$tumor_type
ggplot(data2,aes(primary_median,metastasis_median,label=tumor_type,color=factor(colour)))+
  xlim(-0.02,0.52)+
  ylim(-0.02,0.52)+
  theme_minimal()+
  theme(legend.position = "none")+
  geom_point(size=4)+
  geom_point(data=data2 %>% filter(qval<0.05),alpha=0.5,size=7,aes(color=data2_sig$colour))+
  #scale_color_manual(values=data2$colour)+
  geom_text(data=data2 %>% filter(qval<0.05),hjust=0.5,vjust=1.5)+
  geom_abline(slope=1,color='grey',size=1)+
  labs(y='Metastatsis median ',x='Primary median ',title = 'Fraction of genome altered',size=15)+
  theme(plot.title = element_text(size = 20),axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  scale_colour_identity()

###
#plot figure2A_fga our result
fga_data

ggplot(fga_data,aes(p_median,m_median,label=subtypes,color=colour))+
  xlim(-0.02,0.52)+
  ylim(-0.02,0.52)+
  theme_minimal()+
  theme(legend.position = "none")+
  geom_point(size=2)+
  geom_point(data=fga_data,alpha=0.5,size=3,aes(color=colour))+
  #scale_color_manual(values=data2$colour)+
  geom_text(data=fga_data,hjust=0.5,vjust=1.5)+
  geom_abline(slope=1,color='grey',size=1)+
  labs(y='Metastatsis median ',x='Primary median ',title = 'Fraction of genome altered',size=15)+
  theme(plot.title = element_text(size = 20),axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  scale_colour_identity()


###############################
#plot figure 2A_wgd
data_wgd<-tableS2A %>%  filter(alteration=='WGD')
data_wgd=assign_color(data_wgd)
data_wgd %>% filter(qval<0.05)

ggplot(data_wgd,aes(primary_pc,metastasis_pc,label=tumor_type,color=factor(colour)))+
  xlim(-0.02,1.02)+
  ylim(-0.02,1.02)+
  theme_minimal()+
  theme(legend.position = "none")+
  geom_point(size=4)+
  geom_point(data=data_wgd %>% filter(qval<0.05),alpha=0.5,size=7,aes(color=colour))+
  scale_color_manual(values=data_wgd$colour)+
  geom_text(data=data_wgd %>% filter(qval<0.05),hjust=0.5,vjust=1.5)+
  geom_abline(slope=1,color='grey',size=1)+
  labs(y='Metastatsis median ',x='Primary median ',title = 'Whole genome duplication',size=15)+
  theme(plot.title = element_text(size = 20),axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  scale_colour_identity()


###############################
#plot figure 2A_tmb origin
data_tmb<-tableS2A %>%  filter(alteration=='tmb')
data_tmb=assign_color(data_tmb)
data_tmb%>% filter(qval<0.05)

ggplot(data_tmb,aes(primary_median,metastasis_median,label=tumor_type,color=colour))+
  xlim(-0.02,10)+
  ylim(-0.02,10)+
  theme_minimal()+
  theme(legend.position = "none")+
  geom_point(size=4)+
  geom_point(data=data_tmb %>% filter(qval<0.05),alpha=0.5,size=7,aes(color=colour))+
  scale_color_manual(values=data_tmb$colour)+
  geom_text(data=data_tmb %>% filter(qval<0.05),hjust=0.5,vjust=1.5)+
  geom_abline(slope=1,color='grey',size=1)+
  labs(y='Metastatsis median mut/Mb',x='Primary median mut/Mb',title = 'Tumor mutational burden',size=15)+
  theme(plot.title = element_text(size = 20),axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  scale_colour_identity()
#plot figure 2A_tmb our result
tmb_data=assign_col(tmb_data)
ggplot(tmb_data,aes(p_median,m_median,label=subtypes,color=colour))+
  xlim(-0.02,10)+
  ylim(-0.02,10)+
  theme_minimal()+
  theme(legend.position = "none")+
  geom_point(size=2)+
  geom_point(data=tmb_data,alpha=0.5,size=3,aes(color=colour))+
  #scale_color_manual(values=data2$colour)+
  geom_text(data=tmb_data,hjust=0.5,vjust=1.5)+
  geom_abline(slope=1,color='grey',size=1)+
  labs(y='Metastatsis median ',x='Primary median ',title = 'Tumor mutational burden',size=15)+
  theme(plot.title = element_text(size = 20),axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  scale_colour_identity()
