#!/home/lyw/software/anaconda3/bin/python

##############################################################
# Script to annotate Taxonomy terms in eggnog-mapper results #
#   based on NCBI Taxonomy database and python ETE3 package  #
# by Yewei on 20221026                                       #
##############################################################

# Load library
import sys
import pandas as pd
import os
from ete3 import NCBITaxa
from pathlib import Path

# Input
# NCBI Taxonomy database directory
tax_db_dir = sys.argv[1]
# NCBI Taxonomy database dfile
db_file = os.path.join(tax_db_dir, "taxa.sqlite")
# eggnog-mapper results file
eggnog = sys.argv[2]
# output directory
outdir = sys.argv[3]
# sample name
sample = sys.argv[4]

ncbi = NCBITaxa(dbfile=db_file)
# Parse eggnog-mapper results
eggnog_df = pd.read_csv(eggnog, sep="\t", header=None, comment='#')
gene_id_list = list(eggnog_df[eggnog_df.columns[0]])
eggnog_seed_list = list(eggnog_df[eggnog_df.columns[1]])
# extract NCBI Taxonomy ID from eggNOG seed ortholog (e.g. "470" from "470.IX87_14445")
tax_id_list = [seed.split('.')[0] for seed in eggnog_seed_list]

# Convert all tax_id to tax_name which stores in a dictionary
tax_id2name_dict = ncbi.get_taxid_translator(tax_id_list)

# Convert NCBI Taxonomy id to name and output
fout = os.path.join(outdir, "{}_taxonomy.xls".format(sample))
with open(fout, 'w') as f:
    f.write("{}\t{}\t{}\t{}\n".format("gene_id", "eggNOG_seed_ortholog",
                                      "NCBI_Taxonomy_id", "NCBI_Taxonomy_name"))
    for i in range(len(gene_id_list)):
        gene = gene_id_list[i]
        seed = eggnog_seed_list[i]
        tax_id = tax_id_list[i]
        if int(tax_id) in tax_id2name_dict.keys():
            tax_name = tax_id2name_dict[int(tax_id)]
            f.write("{}\t{}\t{}\t{}\n".format(gene, seed, tax_id, tax_name))

