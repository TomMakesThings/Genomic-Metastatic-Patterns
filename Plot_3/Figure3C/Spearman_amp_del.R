library(readxl)
library(DescTools)
library(ggplot2)

setwd("~/AllMphil/GenomicsAS3/")

table_s3b <- read_excel("supplementary_table/1-s2.0-S0092867422000034-mmc3.xlsx", sheet = 2, skip = 2) #my target

##############################use cbioportal cna  to try amp & del#####################
cna3c <- read.table("cBioPortal/cna_3c.txt", sep = "\t", header = TRUE) #25775    19
mut3c <- read.table("cBioPortal/mut_3c.txt", sep = "\t", header = TRUE) #25775     19




#map tumor subtype info into cna10 & mut10 using table "data_clinical_sample.txt"

clinical_sample <- read.table("data/data_clinical_sample.txt", sep = "\t",header = TRUE)
sampleID_subtype <- clinical_sample[c("SAMPLE_ID", "PATIENT_ID","SAMPLE_TYPE","SUBTYPE", "SUBTYPE_ABBREVIATION",
                                      "MET_SITE_COUNT")] ##get useful column

add_subtype <- function(df, sampleID_subtype){
    if(all(df[,2]==sampleID_subtype$SAMPLE_ID)){
        merged_df <- cbind.data.frame(df[,1:2], SUBTYPE_ABBREVIATION = sampleID_subtype$SUBTYPE_ABBREVIATION,
                                      SUBTYPE = sampleID_subtype$SUBTYPE,
                                      SAMPLE_TYPE = sampleID_subtype$SAMPLE_TYPE,
                                      MET_SITE_COUNT = sampleID_subtype$MET_SITE_COUNT,
                                      PATIENT_ID = sampleID_subtype$PATIENT_ID,
                                      df[,-c(1,2)])
        return(merged_df)
    }else{
        stop("sampleID and subtype order not mapped! sort and check!")
    }
    
}

cna3c <- add_subtype(cna3c, sampleID_subtype)
mut3c <- add_subtype(mut3c, sampleID_subtype)




##make new amp dataframe: in all genes, only keep count cna > 0 data, cna <0 (del) converted to 0

amp_generate <- function(df){
    
    w <- ncol(df)
    
    for (i in 8:w){ #gene starts from 8th column
        df[which(df[, i]<0), i] <- 0
    }
    
    return(df)
}

del_generate <- function(df){
    w <- ncol(df)
    
    for (i in 8:w){ #gene starts from 8th column
        df[which(df[, i]>0), i] <- 0
        df[,i] <- abs(df[,i])
    }
    
    return(df)
}

cna3c_amp <- amp_generate(cna3c) #all del cna converted to 0
cna3c_del <- del_generate(cna3c)

##merge CCND1/FGF 19
cna3c_amp <- mutate(cna3c_amp, CCND1_FGF19 = CCND1 + FGF19)
cna3c_del <- mutate(cna3c_del, CCND1_FGF19 = CCND1 + FGF19)

subtype_spearman <- function(df, alt="Amplification"){
    
    df <- subset(df, df$MET_SITE_COUNT > 0) #Inferred from "n_pts", so we're only looking at patients with MB
    
    #get the list of cancer subtypes
    all_subtypes <- unique(df$SUBTYPE)
    all_genes <- colnames(df[ ,-(1:7)])
    
    #spearman returns NA when two lists are the same / one of them is pure 0.
    output <- data.frame(tumor_type=rep(all_subtypes, each=length(all_genes)), 
                         genes=all_genes,
                         alteration=paste0(all_genes,"_",alt),
                         n_pts=NA,
                         rho=NA,
                         lwr.ci=NA,
                         upr.ci=NA,
                         pvalue=NA,
                         adj.qvalue=NA,
                         sig.afterSpearman=NA)
    
    
    #calculate spearman coef for each subtype
    for (i in 1:nrow(output)){
        
        
        tumor <- output$tumor_type[i]
        gene <- output$genes[i]
        
        tumor_df <- subset(df, df$SUBTYPE == tumor) 
        
        #timer
        if (i %% 100 == 0){
            print(paste0("current gene:",gene))
        }
        #run spearman
        run <- SpearmanRho(x=tumor_df[,gene], y=tumor_df$MET_SITE_COUNT, conf.level = 0.95)
        output[i, c("rho", "lwr.ci", "upr.ci")] <- unname(run)
        
        #extract p-value
        pvalue <- cor.test(tumor_df[,gene],
                           tumor_df$MET_SITE_COUNT,
                           method = "spearman",
                           exact = FALSE)$p.value
        output[i, "pvalue"] <- pvalue
        
        #add n_pts
        output[i, "n_pts"] <- nrow(tumor_df)
        
    }
    
    # Adjust p-values for false discovery rate
    output$adj.qvalue <- p.adjust(output$pvalue, method = "fdr")
    
    # Mark significant genes for each subtype
    output$sig.afterSpearman <- ((output$pvalue < 0.05) & (output$adj.qvalue < 0.05))
    
    return(output)
    
}

##Run linear model for each subtype: MET_SITE_COUNT ~ SAMPLE_TYPE + CNA_AMP

run_lm <- function(df_spearman, cna_df){ #now cna_df is cna10_amp; the df you download from Cbioportal
    
    all_subtypes <- unique(cna_df$SUBTYPE) #50 subtypes
    all_genes <- unique(df_spearman$genes)
    
    df_spearman$lm.pvalue <- NA
    df_spearman$lm.sig <- NA
    
    
    for (i in all_subtypes){
        subtype_cna <- subset(cna_df, cna_df$SUBTYPE == i)
        #only MB>0 #new
        subtype_cna <- subset(subtype_cna, subtype_cna$MET_SITE_COUNT>0)
        
        regression_coefficients <- rep(NA, length(all_genes))
        coefficient_p_value <- rep(NA, length(all_genes))
        
        for (j in 1:length(all_genes)){ 
            
            
            regression_results <- lm(as.formula(paste("MET_SITE_COUNT ~ SAMPLE_TYPE +", all_genes[j])),
                                     data = subtype_cna) 
            
            # Get the p-value for the CNA of gene coefficient and check if significant
            
            regression_coefficients <- data.frame(summary(regression_results)$coefficients)
            coefficient_p_value[j] <- regression_coefficients[all_genes[j],]$Pr...t.. 
            
            df_spearman[which(df_spearman$tumor_type==i), ]$lm.pvalue <- coefficient_p_value
            
            
        }
    }
    
    df_spearman$lm.sig <- df_spearman$lm.pvalue < 0.05
    
    #adjust lm.pvalue #new
    #df_spearman$lm.adj.pvalue <- p.adjust(df_spearman$lm.pvalue, method = "fdr")
    
    df_spearman$final.sig <- df_spearman$lm.sig & df_spearman$sig.afterSpearman
    
    return(df_spearman)
}


amp_spearman <- subtype_spearman(cna3c_amp, alt = "Amplification")
amp_lm <- run_lm(amp_spearman, cna3c_amp)

del_spearman <- subtype_spearman(cna3c_del, alt = "Deletion")
del_lm <- run_lm(del_spearman, cna3c_del)



amp_spearman_sig <- subset(amp_lm, amp_lm$final.sig)
del_spearman_sig <- subset(del_lm, del_lm$final.sig)

#remove merged gene : CCND1 and FGF19
out_amp <- subset(amp_spearman_sig,!amp_spearman_sig$genes %in% c("CCND1","FGF19"))
out_del <- subset(del_spearman_sig,!del_spearman_sig$genes %in% c("CCND1","FGF19"))


write.table(out_amp, file="plotdata/out_amp.txt",sep = "\t")

write.table(out_del, file = "plotdata/out_del.txt",sep = "\t")






