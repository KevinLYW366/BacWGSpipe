# Script to extract statistics from fastp json output

import sys
import json
import os

# Input
input_dir = sys.argv[1]
sample_list_file = sys.argv[2]

# Output
output_file = sys.argv[3]

# Read sample list
with open(sample_list_file, 'r') as f:
    sample_list = [s.strip() for s in f.readlines()]

# Read fastp json output
with open(output_file, 'w') as fout:
    fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("sample", "clean_reads_length_bp", "raw_data_bp",
                                                         "filtered_reads_percentage", "clean_data_bp",
                                                         "clean_data_gc", "clean_data_q20", "clean_data_q30"))
    for s in sample_list:
        json_file = os.path.join(input_dir, "{}/fastp/{}_fastp.json".format(s, s))
        with open(json_file, 'r') as j:
            data = json.load(j)
            mean_length = "({}:{})".format(data["summary"]["after_filtering"]["read1_mean_length"],
                                           data["summary"]["after_filtering"]["read2_mean_length"])
            raw_data_read = int(data["summary"]["before_filtering"]["total_reads"])
            raw_data_bp = int(data["summary"]["before_filtering"]["total_bases"])
            clean_data_read = int(data["summary"]["after_filtering"]["total_reads"])
            clean_data_bp = int(data["summary"]["after_filtering"]["total_bases"])
            clean_data_gc = round(float(data["summary"]["after_filtering"]["gc_content"])*100, 2)
            clean_data_q20 = round(float(data["summary"]["after_filtering"]["q20_rate"])*100, 2)
            clean_data_q30 = round(float(data["summary"]["after_filtering"]["q30_rate"])*100, 2)
            filter_read = round((raw_data_read - clean_data_read) / raw_data_read * 100, 2)
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(s, mean_length, raw_data_bp,
                                                                 filter_read, clean_data_bp,
                                                                 clean_data_gc, clean_data_q20, clean_data_q30))




