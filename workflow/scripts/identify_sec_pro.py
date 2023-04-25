# Script to identify secretory proteins based on SignalP and TMHMM results

import sys
import pandas as pd

signalp_results = sys.argv[1]
tmhmm_results = sys.argv[2]
output = sys.argv[3]

# load data
signalp_df = pd.read_csv(signalp_results, comment="#", header=None, sep="\t")
tmhmm_df = pd.read_csv(tmhmm_results, header=None, sep="\t")

signalp_df_sub = signalp_df.iloc[:, [0, 8]]
signalp_df_sub.columns = ["gene_id", "signal_peptides_signalp"]
tmhmm_df_sub = tmhmm_df.iloc[:, [0, 5]]
tmhmm_df_sub.columns = ["gene_id", "topology_tmhmm"]

# modify gene ids in signalp results
gene_id_new = signalp_df_sub.loc[:, "gene_id"].apply(lambda x: x.split(" ")[0])
signalp_df_sub_new = signalp_df_sub.copy()
signalp_df_sub_new.loc[:, "gene_id"] = gene_id_new

# merge two results
merged_df = signalp_df_sub_new.merge(tmhmm_df_sub, on="gene_id", sort=True)

# Secretory proteins:
# 1. with signal peptides predicted in SignalP results
# 2. Topology=o (without transmembrane helices and outside the membrane) in TMHMM results
merged_df["sectory_protein"] = merged_df.apply(lambda row: "sec_pro" if (row.signal_peptides_signalp != "" and row.topology_tmhmm == "Topology=o") else "non-sec_pro", axis=1)

merged_df.to_csv(output, header=True, index=False, sep="\t")
