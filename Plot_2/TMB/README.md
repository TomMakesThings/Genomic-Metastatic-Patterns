
# Tumor mutational burden (TMB)

TMB was calculated for each sample as the total number of nonsynonymous mutations, divided by the number of bases sequenced. 

## Process of calculation:

1. Check the sequence information of 3 platforms used in this project. https://github.com/cBioPortal/datahub/tree/master/reference_data/gene_panels
2. Assign "number of bases sequenced" to each sampleID, based on the platform each sample was sequenced on. 
3. Filter and count the total number of "nonsynonymous mutations" for each sample.
4. Divide "total number of nonsynonymous mutations" by "number of bases sequenced" to get TMB.

For unknown reasons, the authors published **two versions of TMB**. TMB on cBioPortal is different from those from Supplementary table S1B. Here following the method
described in paper, we got TMB values the same as **cBioPortal version**. It appears that they adjusted both the number of mutations and length of sequence with
method not provided. Therefore we could only obtain **cBioPortal version** with their open data. I will email Kate to ask about their inconsistence of data. 

I have attached both versions of TMB in the uploaded **/Plot_2/TMB/calculated_TMB.csv** file. In this csv file, column **"Our_TMB"** is the TMB we obtained, which is the same as TMB published on 
cBioPortal. Column **"Sup_TMB"** is the TMB from Supplementary table S1B. 

Plese do not hesitate to use "Sup_TMB" if "Our_TMB" does not generate similar plots in later analysis. 
