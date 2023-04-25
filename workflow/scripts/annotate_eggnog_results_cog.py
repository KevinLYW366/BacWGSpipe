#!/home/lyw/software/anaconda3/bin/python

#########################################################
# Script to annotate COG terms in eggnog-mapper results #
# by Yewei on 20221025                                  #
#########################################################

# Load library
import sys
import pandas as pd
import os

# Input
# COG database directory
cog_db_dir = sys.argv[1]
# COG database file with COG descriptions (cog-20.def.tab)
cog_def_file = os.path.join(cog_db_dir, "cog-20.def.tab")
# COG database file with descriptions of COG functional categories (fun-20.tab)
cog_fun_file = os.path.join(cog_db_dir, "fun-20.tab")
# eggnog-mapper results file
eggnog = sys.argv[2]
# output directory
outdir = sys.argv[3]
# sample name
sample = sys.argv[4]

# Read COG database file with COG descriptions (cog-20.def.tab)
cog_def_df = pd.read_csv(cog_def_file, sep="\t", header=None, encoding='cp1252', index_col=0).iloc[:, 0:2]
cog_def_df.columns = ["cog_category", "cog_description"]

# Read COG database file with descriptions of COG functional categories (fun-20.tab)
cog_fun_df = pd.read_csv(cog_fun_file, sep="\t", header=None, index_col=0).iloc[:, 1]

# Parse eggnog-mapper results into a dictionary: gene_id - COG_id
eggnog_df = pd.read_csv(eggnog, sep="\t", header=None, comment='#')
gene_id_list = list(eggnog_df[eggnog_df.columns[0]])
eggNOG_OGs_list = list(eggnog_df[eggnog_df.columns[4]])
# extract root COG IDs from eggNOG_OGs for each gene
cog_list = []
for og in eggNOG_OGs_list:
    root_og = og.split(",")[0]
    # extract something like COG2059@1|root
    if root_og.startswith("COG"):
        cog = root_og.split("@")[0]
    else:
        cog = ""
    cog_list.append(cog)
gene_cog_dict = {gene_id_list[i]: cog_list[i] for i in range(len(gene_id_list))}

# Fetch information for each cog id and output
fout1 = os.path.join(outdir, "{}_cog.xls".format(sample))
with open(fout1, 'w') as f:
    # header line
    f.write("{}\t{}\t{}\t{}\t{}\n".format("gene_id", "COG_id", "COG_description",
                                      "COG_functional_category", "COG_functional_category_description"))
    # fetch information from COG database
    for gene in list(gene_cog_dict.keys()):
        cog = gene_cog_dict[gene]
        if cog:
            # in case some COG IDs are obsolete
            if cog in cog_def_df.index:
                cog_desc = cog_def_df.loc[cog, "cog_description"]
                cog_cat = cog_def_df.loc[cog, "cog_category"]
                cog_cat_desc = ";".join([cog_fun_df[c] for c in cog_cat])
                f.write("{}\t{}\t{}\t{}\t{}\n".format(gene, cog, cog_desc,
                                                      cog_cat, cog_cat_desc))


