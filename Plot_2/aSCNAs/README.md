# ASCETS
Using segmented copy number data, chromosome arm-level copy number alterations were computed with the tool ASCETS. This was run with default parameters and hg19 chromosome arm genomic coordinates.

ASCETS can be run from the script RunASCETS.R or via command line:
```
Rscript ASCETS/run_ascets.R -i Data/data_cna_hg19.seg -c genomic_arm_coordinates_hg19.txt -o sample_output
```
