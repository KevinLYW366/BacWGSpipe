# Script to merge rgi and abricate card results

import sys
import pandas as pd
import os
import numpy as np

rgi_result_file = sys.argv[1]
abricate_result_file = sys.argv[2]
output_file = sys.argv[3]
card_aro_index_file = sys.argv[4]

# Extract info from abricate card results
abricate_result_file_tmp = abricate_result_file + ".tmp"
with open(abricate_result_file_tmp, 'w') as fout:
    fout.write("sample\tgene\tcontig\tstart\tend\tidentity\taccession\tresistance\n")
    with open(abricate_result_file) as fin:
        buf = fin.readline()
        buf = fin.readline()
        while buf:
            items = buf.strip().split("\t")
            sample = items[0]
            gene = items[5]
            contig = items[1]
            start = items[2]
            end = items[3]
            identity = items[10]
            accession = items[12].split(":")[0]
            resistance = items[14]
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, gene, contig, start, end,
                                                             identity, accession, resistance))
            buf = fin.readline()

# Extract info from rgi card results
card_aro_index_df = pd.read_csv(card_aro_index_file, header=0, sep="\t")


def convert_aro2accession(aro):
    aro_reformat = "ARO:{}".format(aro)
    return list(card_aro_index_df.loc[card_aro_index_df['ARO Accession'] == aro_reformat]["DNA Accession"])[0]


rgi_result_file_tmp = rgi_result_file + ".tmp"
with open(rgi_result_file_tmp, 'w') as fout:
    fout.write("sample\tgene\tcontig\tstart\tend\tidentity\taccession\tresistance\n")
    with open(rgi_result_file) as fin:
        buf = fin.readline()
        buf = fin.readline()
        while buf:
            items = buf.strip().split("\t")
            sample = items[0]
            gene = items[9]
            contig = items[2]
            start = items[3]
            end = items[4]
            identity = items[10]
            aro = items[11]
            accession = convert_aro2accession(aro)
            resistance = items[15]
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, gene, contig, start, end,
                                                                 identity, accession, resistance))
            buf = fin.readline()

# Merge rgi and abricate card results
rgi_result_df = pd.read_csv(rgi_result_file_tmp, sep="\t", header=0)
#rgi_result_df.contig = rgi_result_df.contig.astype(int)
merge_result_df = rgi_result_df.copy()
merge_result_df["tool"] = "RGI"

with open(abricate_result_file_tmp, 'r') as f:
    buf = f.readline()
    buf = f.readline()
    while buf:
        line = buf.strip()
        items = line.split("\t")
        # two types of contig format: "1" "2"... or "chr1" "chr2"...
        if merge_result_df["contig"].dtype == np.int64:
            items[2] = int(items[2])
        items[3] = int(items[3])
        items[4] = int(items[4])
        sample = items[0]
        gene = items[1]
        accession = items[6]
        if (accession in list(rgi_result_df.loc[rgi_result_df["sample"] == sample]["accession"])) & (gene in list(rgi_result_df.loc[rgi_result_df["sample"] == sample]["gene"])):
            merge_result_df.loc[(merge_result_df["sample"] == sample) & (merge_result_df["accession"] == accession) & (merge_result_df["gene"] == gene), "tool"] = "RGI; Abricate"
        else:
            #items[2] = int(items[2])
            items.append("Abricate")
            merge_result_df.loc[len(merge_result_df)] = items
        buf = f.readline()

merge_result_df = merge_result_df.sort_values(by=["sample", "contig", "start"], ascending=True)
merge_result_df.to_csv(output_file, sep="\t", header=True, index=False)

# Remove tmp files
os.remove(abricate_result_file_tmp)
os.remove(rgi_result_file_tmp)

