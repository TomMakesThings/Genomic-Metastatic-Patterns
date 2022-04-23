import pandas as pd
import numpy as np
import os
import glob
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerExtractor import sigpro as sig

# Set working directory
directory = "C:\\Users\\redds\\Documents\\GitHub\\Genomics-II-Group"
os.chdir(directory)
# Set directory to store all subtype MAFs
all_subtypes_directory = directory + "\\COSMIC\\MAFs\\Subtypes\\"

# Set the reference genome from the locally downloaded file
genInstall.install("GRCh37", offline_files_path = directory + "\\COSMIC\\")

# Find all subtype output SBS 96 mutational matrices
for file_path in glob.glob(all_subtypes_directory + '**/output/SBS/*.SBS96.all'):
    # Get the file name and directory
    split_path = file_path.split("\\")
    matrix_file = split_path[-1].split(".")[0]
    matrix_directory = 'C:\\' + '\\'.join(split_path[1:len(split_path) - 1])

# Hypermutated
file_path = all_subtypes_directory + 'Colorectal_Hypermutated\\output\\SBS\\CRC_HM_MAF.SBS96.all'

sig.sigProfilerExtractor("matrix", "results", file_path,
                         reference_genome = "GRCh37",
                         minimum_signatures = 1, maximum_signatures = 10,
                         nmf_replicates = 100, cpu = -1, gpu = False)

# wget.download("ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerExtractor/Example_data/21BRCA.zip",
#               directory + "\\COSMIC\\21BRCA.zip")
