# Script to summarize genome assembly statistics

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import os
import re

sample_list_file = sys.argv[1]
indir = sys.argv[2]
outdir = sys.argv[3]
# hybrid,long,short
mode = sys.argv[4]

with open(sample_list_file, 'r') as f:
    sample_list = [s.strip() for s in f.readlines()]

outfile = "{}/all_assembly_stats.xls".format(outdir)
if mode == "hybrid":
    with open(outfile, 'w') as fout:
        fout.write("sample_id\ttype\tcontig_id\tsize_bp\tgc_percent\tcircular\n")
        for sample in sample_list:
            infile = "{}/{}/assembly/{}.fasta".format(indir, sample, sample)
            with open(infile) as fin:
                records = SeqIO.parse(fin, "fasta")
                for record in records:
                    if record.id.startswith("chr") or record.id.startswith("plas"):
                        items = record.description.split(" ")
                        if record.id.startswith("chr"):
                            type = "Chromosome"
                        else:
                            type = "Plasmid"
                        contig_id = items[0]
                        size_bp = items[1].replace("length=", "")
                        gc_percent = str(round(GC(record.seq), 2))
                        circular = "non-circular"
                        if len(items) >= 4:
                            if items[3] == "circular=true":
                                circular = "circular"
                        fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, type, contig_id, size_bp, gc_percent, circular))
elif mode == "short" or mode == "long":
    with open(outfile, 'w') as fout:
        fout.write("sample\tcontig_num\ttotal_length\tN50\tN90\tmax_length\tgc\n")
        for sample in sample_list:
            quast_report_file = os.path.join(indir, "quast/{}/report.tsv".format(sample))
            with open(quast_report_file, 'r') as fin:
                content = fin.read()
                # Number of contigs
                m = re.search(r'# contigs\s+([0-9]+)\n', content)
                contig_num = m.group(1)
                # Total length
                m = re.search(r'Total length\s+([0-9]+)\n', content)
                total_length = m.group(1)
                # N50
                m = re.search(r'N50\s+([0-9]+)\n', content)
                n50 = m.group(1)
                # N90
                m = re.search(r'N90\s+([0-9]+)\n', content)
                n90 = m.group(1)
                # MAX length
                m = re.search(r'Largest contig\s+([0-9]+)\n', content)
                max_length = m.group(1)
                # GC percentage
                m = re.search(r'GC \(%\)\s+([0-9.]+)\n', content)
                gc = m.group(1)
                fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, contig_num, total_length, n50, n90, max_length, gc))
