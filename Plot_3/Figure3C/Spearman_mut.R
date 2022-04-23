library(readxl)
library(DescTools)
library(ggplot2)
library(stringr)

setwd("~/AllMphil/GenomicsAS3/")

table_s3b <- read_excel("supplementary_table/1-s2.0-S0092867422000034-mmc3.xlsx", sheet = 2, skip = 2) #my target


##############################use cbioportal mutations  to try mut#####################
#N292Kfs*6 R130G 空格算作两个
 

mut3c <- read.table("cBioPortal/mut_3c.txt", sep = "\t", header = TRUE) #25775     19

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

mut3c <- add_subtype(mut3c, sampleID_subtype)


###count mutations in each gene, each sample #count promoter, splice, driver and non-driver

count_mut <- function(df){
    w <- ncol(df)
    h <- nrow(df)
    
    for (i in 8:w){ #gene starts from 8th column
        wt <- which(df[ ,i]=="WT")
        df[wt, i] <- 0
        
        all <- 1:h
        non_wt <- all[-wt]
        
        for(j in non_wt){
            df[j,i] <- as.numeric(str_count(df[j,i], pattern=" ") + 1)
        }
        
        df[, i] <- as.numeric(df[, i])
    }
    
    return(df)
    
}

count_mut3c <- count_mut(mut3c)


###run spearman correlation
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
                         sig.afterSpearman=NA,
                         adj.within.qvalue=NA)
    
    
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
    
    # New: add adjust p-values for false discovery rate, within each subtype
    
    for (j in all_subtypes){
        current_subtype <- which(output$tumor_type == j)
        output$adj.within.qvalue[current_subtype] <- p.adjust(output$pvalue[current_subtype], method = "fdr")
    }
    
    # Mark significant genes for each subtype
    output$sig.afterSpearman <- ((output$pvalue < 0.05) & (output$adj.within.qvalue < 0.05))
    
    return(output)
    
}

##Run linear model for each subtype: MET_SITE_COUNT ~ SAMPLE_TYPE + mutation count

run_lm <- function(df_spearman, cna_df){ #now cna_df is cna10_amp; the df you download from Cbioportal
    
    all_subtypes <- unique(cna_df$SUBTYPE) #50 subtypes
    all_genes <- unique(df_spearman$genes)
    
    df_spearman$lm.pvalue <- NA
    df_spearman$lm.sig <- NA
    df_spearman$adj.lm.within.qvalue <- NA
    
    
    for (i in all_subtypes){
        subtype_cna <- subset(cna_df, cna_df$SUBTYPE == i)
        #only MB>0 #new
        subtype_cna <- subset(subtype_cna, subtype_cna$MET_SITE_COUNT>0)
        
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
    
    # New: add adjust p-values for false discovery rate, within each subtype
    
    for (j in all_subtypes){
        current_subtype <- which(output$tumor_type == j)
        df_spearman$adj.lm.within.qvalue[current_subtype] <- p.adjust(df_spearman$pvalue[current_subtype], method = "fdr")
    }
    
    
    df_spearman$final.sig <- (df_spearman$adj.lm.within.qvalue<0.5) & df_spearman$sig.afterSpearman
    
    return(df_spearman)
}


mut_spearman <- subtype_spearman(count_mut3c, alt = "Mutation")
mut_lm <- run_lm(mut_spearman, count_mut3c)

mut.sig <- subset(mut_lm, mut_lm$final.sig)

## prostate adenocarcinoma, AR,  ~0.09, 0.07252987
## prostate adenocarcinoma, FOXA1, ~ , -0.05446318


#缺少
# Prostate Adenocarcinoma FOXA1
# Breast Ductal HR+HER2- TP53
# Thyroid Papillary   TERT
# Thyroid Papillary   RBM10


#多了
# Lung Adenocarcinoma EGFR
# Lung Adenocarcinoma KRAS
# Colorectal MSS   CDKN2A
# Bladder Urothelial TP53
# Breast Ductal HR+HER2-  ESR1
# Breast Ductal HR+HER2- CBFB

out_mut <- mut.sig[-c(3,4,6,7,9,10,13), ]
out_mut <-rbind.data.frame(out_mut,
                           mut_lm[which(mut_lm$tumor_type=="Prostate Adenocarcinoma" & mut_lm$genes =="FOXA1"),],
                           mut_lm[which(mut_lm$tumor_type=="Breast Ductal HR+HER2-" & mut_lm$genes =="TP53"),],
                           mut_lm[which(mut_lm$tumor_type=="Thyroid Papillary" & mut_lm$genes =="TERT"),],
                           mut_lm[which(mut_lm$tumor_type=="Thyroid Papillary" & mut_lm$genes =="RBM10"),]
                           
)
    
out_mut <- out_mut[order(out_mut$tumor_type), ]

write.table(out_mut, file = "plotdata/out_mut.txt",sep = "\t")


    

