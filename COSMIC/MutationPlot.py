import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# Set working directory
directory = "C:\\Users\\redds\\Documents\\GitHub\\Genomics-II-Group"

# Open the example data
example_maf_df = pd.read_table(directory + "\\COSMIC\\MAFs\\Example_MAF\\example.maf")
example_maf_df.head(5)

# Open the mutation data
annotated_mutation_df = pd.read_table(directory + "\\Plot_2\\OncoKB\\oncokb_annotated_mutations.txt")
mutation_df = pd.read_table(directory + "\\Data\\data_mutations.txt")
# Add labels to mutation data
mutation_df = mutation_df.join(annotated_mutation_df.loc[:,["Sample_ID", "Subtype", "Sample_Type"]])

mutation_df.head(5)

# List of all unique tumour subtypes
all_subtypes = [x for x in np.unique(list(mutation_df["Subtype"])) if x != 'nan']

# Calculate the average number of mutations in primary/metastatic samples for each subtype
subtype_primary_mutations = {}
subtype_metastasis_mutations = {}

for subtype in all_subtypes:
    # Get the rows in sample data for the subtype
    subtype_samples = mutation_df.loc[mutation_df["Subtype"] == subtype]
    
    primary_samples = subtype_samples.loc[subtype_samples["Sample_Type"] == "Primary"]
    metastasis_samples = subtype_samples.loc[subtype_samples["Sample_Type"] == "Metastasis"]
    
    # Count the average number of mutations per sample in the subtype
    subtype_primary_mutations[subtype] = np.mean(list(Counter(list(primary_samples["Sample_ID"])).values()))
    subtype_metastasis_mutations[subtype] = np.mean(list(Counter(list(metastasis_samples["Sample_ID"])).values()))

mutation_count_df = pd.DataFrame({"Subtype": all_subtypes,
                                  "Primary": list(subtype_primary_mutations.values()),
                                  "Metastasis": list(subtype_metastasis_mutations.values())})
mutation_count_df.head(5)

mutation_count_bar = mutation_count_df.plot(x = "Subtype",
                                            kind = "bar",
                                            stacked = False,
                                            title = "Average Number of Mutations per Sample")
plt.xlabel("Tumour Type")  
plt.ylabel("Number of Mutations")
plt.savefig(directory + '\\COSMIC\\Mutation_Count_Barchart.pdf', bbox_inches = 'tight')
plt.show()
