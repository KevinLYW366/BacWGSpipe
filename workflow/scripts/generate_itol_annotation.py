# Script to generate itol annotation files

# Library
import sys
import os
import random

# Input
amr_input_file = sys.argv[1]
vf_input_file = sys.argv[2]
mlst_input_file = sys.argv[3]

# Output
amr_output_file = sys.argv[4]
vf_output_file = sys.argv[5]
mlst_output_file = sys.argv[6]

# AMR annotation
with open(amr_input_file, 'r') as fin:
    # header
    buf = fin.readline()
    amr_list = buf.strip().split("\t")
    amr_n = len(amr_list)
    # matrix
    amr_matrix_dict = {}
    buf = fin.readline()
    while buf:
        items = buf.split("\t")
        sample = items[0]
        binary_list = ["0" if x in ["\n", ""] else "1" for x in items[1:]]
        amr_matrix_dict[sample] = binary_list
        buf = fin.readline()

with open(amr_output_file, 'w') as fout:
    fout.write("DATASET_BINARY\n")
    fout.write("SEPARATOR TAB\n")
    fout.write("DATASET_LABEL\t{}\n".format("AMR"))
    fout.write("COLOR\t{}\n".format("#20854eff"))
    fout.write("\n")
    fout.write("SHOW_LABELS\t1\n")
    fout.write("FIELD_SHAPES\t{}\n".format("\t".join(["2"] * amr_n)))
    fout.write("FIELD_COLORS\t{}\n".format("\t".join(["#f95700ff"] * amr_n)))
    fout.write("FIELD_LABELS\t{}\n".format("\t".join(amr_list)))
    fout.write("\n")
    fout.write("LEGEND_TITLE\tAMR\n")
    fout.write("LEGEND_SHAPES\t2\n")
    fout.write("LEGEND_LABELS\tExist\n")
    fout.write("LEGEND_COLORS\t#f95700ff\n")
    fout.write("\n")
    fout.write("LABEL_ROTATION\t45\n")
    fout.write("LABEL_SHIFT\t10\n")
    fout.write("\n")
    fout.write("DATA\n")
    for sample in amr_matrix_dict:
        fout.write("{}\t{}\n".format(sample, "\t".join(amr_matrix_dict[sample])))

# VF annotation
with open(vf_input_file, 'r') as fin:
    # header
    buf = fin.readline()
    vf_list = buf.strip().split("\t")
    vf_n = len(vf_list)
    # matrix
    vf_matrix_dict = {}
    buf = fin.readline()
    while buf:
        items = buf.split("\t")
        sample = items[0]
        binary_list = ["0" if x in ["\n", ""] else "1" for x in items[1:]]
        vf_matrix_dict[sample] = binary_list
        buf = fin.readline()

with open(vf_output_file, 'w') as fout:
    fout.write("DATASET_BINARY\n")
    fout.write("SEPARATOR TAB\n")
    fout.write("DATASET_LABEL\t{}\n".format("VF"))
    fout.write("COLOR\t{}\n".format("#20854eff"))
    fout.write("\n")
    fout.write("SHOW_LABELS\t1\n")
    fout.write("FIELD_SHAPES\t{}\n".format("\t".join(["2"] * vf_n)))
    fout.write("FIELD_COLORS\t{}\n".format("\t".join(["#00a4ccff"] * vf_n)))
    fout.write("FIELD_LABELS\t{}\n".format("\t".join(vf_list)))
    fout.write("\n")
    fout.write("LEGEND_TITLE\tVF\n")
    fout.write("LEGEND_SHAPES\t2\n")
    fout.write("LEGEND_LABELS\tExist\n")
    fout.write("LEGEND_COLORS\t#f95700ff\n")
    fout.write("\n")
    fout.write("LABEL_ROTATION\t45\n")
    fout.write("LABEL_SHIFT\t10\n")
    fout.write("\n")
    fout.write("DATA\n")
    for sample in vf_matrix_dict:
        fout.write("{}\t{}\n".format(sample, "\t".join(vf_matrix_dict[sample])))

# MLST annotation
with open(mlst_input_file, 'r') as fin:
    # header
    buf = fin.readline()
    # matrix
    mlst_dict = {}
    buf = fin.readline()
    while buf:
        items = buf.split("\t")
        sample = items[0]
        scheme = items[1]
        st = "{} - {}".format(scheme, "ST" + items[2])
        mlst_dict[sample] = st
        buf = fin.readline()

# generate colors palette
mlst_dedup = list(set(mlst_dict.values()))
n_mlst_dedup = len(mlst_dedup)
color_palette_list = []
for j in range(n_mlst_dedup):
    color_palette_list.append("#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]))
mlst_color_dict = {mlst_dedup[i]: color_palette_list[i] for i in range(n_mlst_dedup)}

with open(mlst_output_file, 'w') as fout:
    fout.write("DATASET_COLORSTRIP\n")
    fout.write("SEPARATOR TAB\n")
    fout.write("DATASET_LABEL\t{}\n".format("MLST"))
    fout.write("COLOR\t{}\n".format("#ff0000"))
    fout.write("\n")
    fout.write("LEGEND_TITLE\tMLST\n")
    fout.write("LEGEND_SHAPES\t1\n")
    fout.write("LEGEND_LABELS\t{}\n".format("\t".join(mlst_color_dict.keys())))
    fout.write("LEGEND_COLORS\t{}\n".format("\t".join(mlst_color_dict.values())))
    fout.write("\n")
    fout.write("BORDER_WIDTH\t2\n")
    fout.write("BORDER_COLOR\t#ffffff\n")
    fout.write("COMPLETE_BORDER\t1\n")
    fout.write("\n")
    fout.write("DATA\n")
    for sample in mlst_dict.keys():
        fout.write("{}\t{}\n".format(sample, mlst_color_dict[mlst_dict[sample]]))
