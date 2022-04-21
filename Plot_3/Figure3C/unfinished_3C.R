library(readxl)
library(DescTools)
library(ggplot2)

setwd("~/AllMphil/GenomicsAS3/")


##############################use 10 genes to try amp & del#####################

cna10 <- read.table("cBioPortal/cna_10.txt", sep = "\t", header = TRUE) #25775    11
mut10 <- read.table("cBioPortal/mut.txt", sep = "\t", header = TRUE) #25775     4



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

cna10 <- add_subtype(cna10, sampleID_subtype)
mut10 <- add_subtype(mut10, sampleID_subtype)




##make new amp dataframe: in all genes, only keep count cna > 0 data, cna <0 (del) converted to 0

amp_generate <- function(df){
    
    w <- ncol(df)
    
    for (i in 8:w){ #gene starts from 8th column
        df[which(df[, i]<0), i] <- 0
    }
    
    return(df)
}

cna10_amp <- amp_generate(cna10) #all del cna converted to 0

subtype_spearman <- function(df){
    
    df <- subset(df, df$MET_SITE_COUNT > 0) #Inferred from "n_pts", so we're only looking at patients with MB
    
    #get the list of cancer subtypes
    all_subtypes <- unique(df$SUBTYPE)
    all_genes <- colnames(df[ ,-(1:7)])
    
    #如果两个list中出现全0， 或者两个list 完全一样，spearman会NA
    output <- data.frame(tumor_type=rep(all_subtypes, each=length(all_genes)), 
                         genes=all_genes,
                         alteration=paste0(all_genes,"_","Amplification"),
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

df_spearman <- subtype_spearman(cna10_amp)

##Run linear model for each subtype: MET_SITE_COUNT ~ SAMPLE_TYPE + CNA_AMP

run_lm <- function(df_spearman, cna_df){ #now cna_df is cna10_amp; the df you download from Cbioportal
    
    all_subtypes <- unique(cna_df$SUBTYPE) #50 subtypes
    all_genes <- unique(df_spearman$genes)
    
    df_spearman$lm.pvalue <- NA
    df_spearman$lm.sig <- NA
    
    
    for (i in all_subtypes){
        subtype_cna <- subset(cna_df, cna_df$SUBTYPE == i)
        
        regression_coefficients <- rep(NA, length(all_genes))
        coefficient_p_value <- rep(NA, length(all_genes))
        
        for (j in 1:length(all_genes)){ #注意这种做法，前面不要remove NA了
            
            
            regression_results <- lm(as.formula(paste("MET_SITE_COUNT ~ SAMPLE_TYPE +", all_genes[j])),
                                     data = subtype_cna) 
            
            # Get the p-value for the CNA of gene coefficient and check if significant
            
            regression_coefficients <- data.frame(summary(regression_results)$coefficients)
            coefficient_p_value[j] <- regression_coefficients[all_genes[j],]$Pr...t.. 
            
            df_spearman[which(df_spearman$tumor_type==i), ]$lm.pvalue <- coefficient_p_value
            
            
        }
    }
    
    df_spearman$lm.sig <- df_spearman$lm.pvalue < 0.05
    
    return(df_spearman)
}


df_lm <- run_lm(df_spearman,cna10_amp)


#Problem: 
#lung adenocarcinoma, FOXA1 (3349) should be: 
#lm.adj.coef	lm.adj.SE	lm.adj.pval
#0.11	0.22	6.3E-01
#But I got a pvalue of 1.115378e-04







