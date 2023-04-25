# Script to extract info from NanoStat results and merge

import sys
import re

sample_list_file = sys.argv[1]
indir = sys.argv[2]
outdir = sys.argv[3]
# mode: raw or clean
mode = sys.argv[4]

with open(sample_list_file, 'r') as f:
    sample_list = [s.strip() for s in f.readlines()]

outfile = "{}/all_{}_NanoStats.xls".format(outdir, mode)
with open(outfile, 'w') as fout:
    fout.write("Sample ID\tNumber of Reads\tTotal bases(bp)\tMean Read Length(bp)\tN50 Read Length(bp)\tMean Read Quality\n")
    for s in sample_list:
        infile = "{}/{}_{}_NanoStats.txt".format(indir, s, mode)
        with open(infile, 'r') as fin:
            content = fin.read()
            # Number of Reads
            m = re.search(r'Number of reads:(.*)\n', content)
            n_reads = m.group(1).strip()
            # Total bases(G)
            m = re.search(r'Total bases:(.*)\n', content)
            total_bases = m.group(1).strip()
            # Mean Read Length(bp)
            m = re.search(r'Mean read length:(.*)\n', content)
            mean_len = m.group(1).strip()
            # N50 Read Length(bp)
            m = re.search(r'Read length N50:(.*)\n', content)
            n50_len = m.group(1).strip()
            # Mean Read Quality
            m = re.search(r'Mean read quality:(.*)\n', content)
            mean_qual = m.group(1).strip()
            fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(s, n_reads, total_bases, mean_len, n50_len, mean_qual))

