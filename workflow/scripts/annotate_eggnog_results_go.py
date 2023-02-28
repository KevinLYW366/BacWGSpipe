#!/home/lyw/software/anaconda3/bin/python

########################################################
# Script to annotate GO terms in eggnog-mapper results #
#  based on python package goatools                    #
# by Yewei on 20221024                                 #
########################################################

# conda install -c bioconda goatools
# conda install docopt
# conda activate goatools

# Load library
import sys
from goatools import obo_parser
import pandas as pd
from goatools.gosubdag.gosubdag import GoSubDag
import os

# Input
# core GO ontology (OBO Format) file
go_obo = sys.argv[1]
# eggnog-mapper results file
eggnog = sys.argv[2]
# output directory
outdir = sys.argv[3]
# sample name
sample = sys.argv[4]


# function to query a GO term in goatools
def query_go(go_id, query_dict):
    if go_id in list(query_dist.keys()):
        term_query = query_dist[go_id]
    else:
        term_query = goadag.query_term(go_id)
        query_dist[go_id] = term_query
    return term_query, query_dict


# function to extract info from goatools query term
def extract_term(term, bp_list, cc_list, mf_list):
    go_level = str(term.level)
    go_class = str(term.namespace)
    go_description = str(term.name)
    # group all GO terms into three level-0 GO terms
    if go_class == "biological_process":
        bp_list.append("{} {}".format(go, go_description))
    elif go_class == "cellular_component":
        cc_list.append("{} {}".format(go, go_description))
    elif go_class == "molecular_function":
        mf_list.append("{} {}".format(go, go_description))
    return go_level, go_class, go_description


# Parse eggnog-mapper results into a dictionary: gene_id - list of GO terms
eggnog_df = pd.read_csv(eggnog, sep="\t", header=None, comment='#')
gene_id_list = list(eggnog_df[eggnog_df.columns[0]])
go_list = list(eggnog_df[eggnog_df.columns[9]])
gene_go_dict = {gene_id_list[i]:go_list[i] for i in range(len(gene_id_list))}

# Read OBO file
goadag = obo_parser.GODag(go_obo, load_obsolete=True, optional_attrs={'consider', 'replaced_by'})

# Annotate GO terms and output
fout1 = os.path.join(outdir, "{}_go.xls".format(sample))
#fout2 = "{}/GO_annotation_level_1.xls".format(outdir)
fout2 = os.path.join(outdir, "{}_go_with_level.xls".format(sample))
with open(fout1, 'w') as f1:
    with open(fout2, 'w') as f2:
        # header line
        f1.write("{}\t{}\t{}\t{}\n".format("gene_id", "biological_process",
                                           "cellular_component", "molecular_function"))
        #f2.write("{}\t{}\t{}\t{}\n".format("gene_id", "GO_id", "GO_class", "GO_description"))
        f2.write("{}\t{}\t{}\t{}\t{}\n".format("gene_id", "GO_id", "GO_class", "GO_level", "GO_description"))
        # annotate GO terms
        query_dist = {}
        #ancestor_dict = {}
        for gene in list(gene_go_dict.keys()):
            go_list = gene_go_dict[gene].split(",")
            bp = []
            cc = []
            mf = []
            for go in go_list:
                if go != "-":
                    term, query_dist = query_go(go, query_dist)
                    if term:
                        # Some GO IDs might be obsolete because of the GO database update
                        # Replace these GO IDs with new GO IDs in term.replaced_by and term.consider
                        if term.is_obsolete:
                            if term.consider or term.replaced_by:
                                go_list_new = list(term.consider)
                                if term.replaced_by and term.replaced_by not in go_list_new:
                                    go_list_new.append(term.replaced_by)
                                for go_new in go_list_new:
                                    term, query_dist = query_go(go_new, query_dist)
                                    if term:
                                        go_level, go_class, go_description = extract_term(term, bp, cc, mf)
                                        # output all GO terms with level
                                        f2.write("{}\t{}\t{}\t{}\t{}\n".format(gene, go_new, go_class,
                                                                               go_level, go_description))
                        else:
                            go_level, go_class, go_description = extract_term(term, bp, cc, mf)
                            # output all GO terms with level
                            f2.write("{}\t{}\t{}\t{}\t{}\n".format(gene, go, go_class, go_level, go_description))
            f1.write("{}\t{}\t{}\t{}\n".format(gene, ";".join(bp), ";".join(cc), ";".join(mf)))
