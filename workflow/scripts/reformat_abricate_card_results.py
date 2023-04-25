# Script to reformat Abricate (CARD) results
# 1. Replace "_" with " ".
# 2. Some AMR gene names are the model name in CARD index file.
#    Replace these gene names with ARO names to match up with RGI results.

import sys
import pandas as pd


def check_gene(gene, accession, card_index_df):
    if gene not in list(card_index_df["ARO Name"]):
        # check if gene is Model Name
        if gene in list(card_index_df["Model Name"]):
            n = 0
            for a in list(card_index_df.loc[card_index_df["Model Name"] == gene]["DNA Accession"]):
                if a == accession:
                    return list(card_index_df.loc[card_index_df["Model Name"] == gene]["ARO Name"])[n]
                n += 1
            return ""
        else:
            return ""
    else:
        return gene


# Input
abricate_results_file = sys.argv[1]
card_index_file = sys.argv[2]

# Output
output_file = sys.argv[3]

# Read CARD index file
card_index_df = pd.read_csv(card_index_file, sep="\t", header=0)

# Reformat Abricate (CARD) results and output
with open(output_file, 'w') as fout:
    with open(abricate_results_file, 'r') as fin:
        # header line
        buf = fin.readline()
        fout.write(buf)
        buf = fin.readline()
        while buf:
            items = buf.split("\t")
            gene = items[5]
            accession = items[12].split(":")[0]
            gene_new = check_gene(gene, accession, card_index_df)
            if gene_new:
                items[5] = gene_new
                fout.write("\t".join(items))
            buf = fin.readline()

