library(ggplot2)
library(egg)
library(cowplot)

setwd("~/AllMphil/GenomicsAS3/")

out_amp <- read.table("plotdata/out_amp.txt",sep = "\t")
out_del <- read.table("plotdata/out_del.txt",sep = "\t")
out_mut <- read.table("plotdata/out_mut.txt",sep = "\t")

#get useful columns
useful_columns <- c("tumor_type", "genes", "alteration", "n_pts","rho","lwr.ci","upr.ci")
df <- rbind.data.frame(out_amp[useful_columns],
                       out_del[useful_columns],
                       out_mut[useful_columns])
df$alt <- sapply(df$alteration, function(x){
    list <- strsplit(x, "_(?!.*_)", perl=TRUE)
    return(unlist(list)[2])
})

df$color.code <- sapply(df$alt, function(x){
    if (x=="Amplification"){
        y <- "red"
    }else if (x=="Mutation"){
        y <- "forestgreen"
    }else{
        y <- "blue"
    }
    return(y)
})

#plot function
plot3c <- function(df, tumor, themecolor){
    df <- subset(df, df$tumor_type == tumor)
    df <- df[order(df$rho, decreasing = FALSE), ]
    df$alteration <- factor(df$alteration, levels = df$alteration)
    #df$genes <- factor(df$genes, levels = df$genes)

    plot <- ggplot(data = df, aes(x=alteration, y=rho, colour=color.code,
                                  ymin = lwr.ci, ymax = upr.ci
                                  ))+
        geom_pointrange(size=0.8)+
        coord_flip()+
        scale_colour_identity()+
        geom_hline(yintercept = 0, color = "white", size = 1)+
        scale_y_continuous(position = "right", limits = c(-0.4, 0.4), 
                           breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4))+
        theme(plot.title = element_text(size = 12, hjust = 0, vjust = 0), # Set title position
              plot.caption.position = "plot", # Set caption position and size
              plot.caption = element_text(hjust = 1, size = 8),
              axis.text.y = element_text(hjust = 1), # Colour labels
              legend.position = "bottom",
              panel.grid.major.x = element_blank(), # Hide x-axis background grid lines
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_line(linetype = "dotted", size = 0.8), # Set style of y-axis background grid lines
              panel.grid.minor.y = element_line(linetype = "dotted", size = 0.8),
              panel.background = element_rect(fill = alpha(themecolor,0.2), color =themecolor))+
        labs(title = paste0(df$tumor_type[1], "(",df$n_pts[1],")"),x="",y="")+
        scale_x_discrete(labels= df$genes)

    return(plot)
}

PA <- plot3c(df, tumor = "Prostate Adenocarcinoma",themecolor = "#be1e2d")
LA <- plot3c(df, tumor = "Lung Adenocarcinoma",themecolor = "#b368d9")
BD <- plot3c(df, tumor = "Breast Ductal HR+HER2-",themecolor = "#e6308e")
CM <- plot3c(df, tumor = "Colorectal MSS",themecolor = "#007eb5")
TP <- plot3c(df, tumor = "Thyroid Papillary",themecolor = "#cccc33")
BDP <- plot3c(df, tumor = "Colorectal MSS",themecolor = "#e6308e")
BL <- plot3c(df, tumor = "Breast Lobular HR+",themecolor = "#e6308e")
BU <- plot3c(df, tumor = "Bladder Urothelial",themecolor = "#F8766D")




final <- plot_grid(PA,LA,BD,CM,TP,BDP,BL,BU, nrow = 8, rel_heights = 6*c(1,4/7,3/7,2/7,3/7,2/7,2/7,2/7))

ggsave("plot/fig3c.png", plot = final, width = 7, height = 14)
ggsave("plot/fig3c.pdf", plot = final, width = 7, height = 14)


