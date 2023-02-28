# Script to extract PLSDB meta info

import pandas as pd
import sys
import os

mash_result_file = sys.argv[1]
plsdb_meta_file = sys.argv[2]
output_file = sys.argv[3]
fasta_file = sys.argv[4]

sample = os.path.basename(mash_result_file).split("_")[0]
plas = os.path.basename(mash_result_file).split("_")[1]
with open(fasta_file, 'r') as f:
    first_line = f.readline()
    query_len = first_line.strip().split(" ")[1].replace("length=", "")

# pandas.read_csv() will throw an error with empty file as input
if os.stat(mash_result_file).st_size == 0:
    mash_result_df = pd.DataFrame(columns=["Identity", "Shared_hashes", "Median_multiplicity", "P-value",
                                           "Plasmid_accession", "Description"])
else:
    mash_result_df = pd.read_csv(mash_result_file, sep="\t", header=None)
    mash_result_df.columns = ["Identity", "Shared_hashes", "Median_multiplicity", "P-value",
                              "Plasmid_accession", "Description"]

plsdb_meta_df = pd.read_csv(plsdb_meta_file, sep="\t", header=0)

# Merge and output
merge_df = mash_result_df.merge(plsdb_meta_df, how="left",
                                left_on=["Plasmid_accession", "Description"],
                                right_on=["ACC_NUCCORE", "Description_NUCCORE"])
merge_df = merge_df.drop(labels=["UID_NUCCORE", "ACC_NUCCORE", "Description_NUCCORE"], axis=1)
merge_df.insert(0, "Query_length", query_len)
merge_df.insert(0, "Plasmid_ID", plas)
merge_df.insert(0, "Sample_ID", sample)

merge_df.sort_values(by=["Identity"], axis=0, ascending=False, inplace=True)
merge_df.to_csv(output_file, sep="\t", header=True, index=False)


