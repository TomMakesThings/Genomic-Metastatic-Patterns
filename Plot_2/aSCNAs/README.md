# ASCETS
Using segmented copy number data, chromosome arm-level copy number alterations were computed with the tool ASCETS. This was run with default parameters and hg19 chromosome arm genomic coordinates.

ASCETS can be run from the script RunASCETS.R or via command line:
```
Rscript ASCETS/run_ascets.R -i Data/data_cna_hg19.seg -c genomic_arm_coordinates_hg19.txt -o sample_output
```

The annotated results of ASCETS are saved to files `sample_arm_level_cna.csv` and `subtype_arm_level_cna.csv`, where the former contains all tumour samples and the latter is averaged across all samples belonging to a cancer subtype. WGD estimates are also included in these files. This was estimated based on the percentage of amplification across a sample's chromosomes. If for a sample no amplification was found and over half the chromsomes did not have a score, then WGD is set as NA.
