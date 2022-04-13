library("readxl")
tableS2A <- read_excel("Table_S2.xlsx",skip = 1)
tableS2A<-data.frame(tableS2A)

tableS1A<-read_excel("Table_S1.xlsx",skip = 1)
tableS1A<-data.frame(tableS1A)

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
data2<-assign_color(data2)
#threshold used for figure2A
data2 %>% filter(qval<0.05)
data2$colour
data2$tumor_type

#plot figure2A_fga
ggplot(data2,aes(primary_median,metastasis_median,label=tumor_type,color=colour))+
  xlim(-0.02,0.52)+
  ylim(-0.02,0.52)+
  theme_minimal()+
  theme(legend.position = "none")+
  geom_point()+
  geom_point(data=data2 %>% filter(qval<0.05),alpha=0.5,size=5,aes(color=colour))+
  scale_color_manual(values=data2$colour)+
  geom_text(data=data2 %>% filter(qval<0.05),hjust=0.5,vjust=1.5)+
  geom_abline(slope=1,color='grey',size=1)+
  labs(y='Metastatsis median ',x='Primary median ',title = 'Fraction of genome altered')
?geom_text
  

###############################
#plot figure 2A_wgd
data_wgd<-tableS2A %>%  filter(alteration=='WGD')
data_wgd=assign_color(data_wgd)
data_wgd %>% filter(qval<0.05)

ggplot(data_wgd,aes(primary_pc,metastasis_pc,label=tumor_type,color=colour))+
  xlim(-0.02,1.02)+
  ylim(-0.02,1.02)+
  theme_minimal()+
  theme(legend.position = "none")+
  geom_point()+
  geom_point(data=data_wgd %>% filter(qval<0.05),alpha=0.5,size=5,aes(color=colour))+
  scale_color_manual(values=data_wgd$colour)+
  geom_text(data=data_wgd %>% filter(qval<0.05),hjust=0.5,vjust=1.5)+
  geom_abline(slope=1,color='grey',size=1)+
  labs(y='Metastatsis median %',x='Primary median %',title = 'Whole genome duplication')


data2
