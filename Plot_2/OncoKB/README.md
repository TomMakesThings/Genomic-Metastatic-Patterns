# OncoKB Actionability Scores
First OncoTree codes were mapped to mutations in the MAF mutation file [Data/data_mutations.txt](https://github.com/TomMakesThings/Genomics-II-Group/blob/main/Data/data_mutations.txt) to create [mutations_with_oncotree_codes.txt](https://github.com/TomMakesThings/Genomics-II-Group/blob/main/Plot_2/OncoKB/mutations_with_oncotree_codes.txt). Then mutations in this file were annotated with actionability scores using OncoKB Annotator. The annotated results are saved in [oncokb_annotated_mutations.zip](https://github.com/TomMakesThings/Genomics-II-Group/blob/main/Plot_2/OncoKB/oncokb_annotated_mutations.zip). For information about the added columns, [click here](https://github.com/TomMakesThings/Genomics-II-Group/tree/main/OncoKB_Annotator#columns-added-in-the-annotation-files).

#### To run OncoKB code on anaconda:

Create conda environment:

`conda create -n oncokb python=3.8`

`conda activate oncokb`

Install dependencies:

`cd C:\Users\redds\Documents\GitHub\Genomics-II-Group`

`pip install -r requirements/common.txt -r requirements/pip3.txt`

Annotate mutations in the somatic MAF file "mutations_with_oncotree_codes.txt"

`python OncoKB_Annotator\MafAnnotator.py -i Plot_2\OncoKB\mutations_with_oncotree_codes.txt -o Plot_2\OncoKB\oncokb_annotated_mutations.txt -b 44ceacf7-c7d2-4b85-ae72-e7ef71c0da55`
