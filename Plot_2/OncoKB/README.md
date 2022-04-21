#### To run OncoKB code on anaconda:

Create conda environment:

`conda create -n oncokb python=3.8`

`conda activate oncokb`

Install dependencies:

`cd C:\Users\redds\Documents\GitHub\Genomics-II-Group`

`pip install -r requirements/common.txt -r requirements/pip3.txt`

Annotate mutations in the somatic MAF file "mutations_with_oncotree_codes.txt"

`python OncoKB_Annotator\MafAnnotator.py -i Plot_2\OncoKB\mutations_with_oncotree_codes.txt -o Plot_2\OncoKB\oncokb_annotated_mutations.txt -b 44ceacf7-c7d2-4b85-ae72-e7ef71c0da55`
