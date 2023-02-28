#!/home/lyw/software/anaconda3/bin/python

# Script to extract the following stats from Bracken report in each sample
#   1. Percentage of reads from target species
#   2. Species with the highest percentage of reads (Species & Reads percentage)

import sys

bracken_report = sys.argv[1]
output = sys.argv[2]
sample_name = sys.argv[3]
genus = sys.argv[4]
species = sys.argv[5]

genus_species = "{} {}".format(genus, species)

#f_out = open(output, 'w')
# without overwriting
f_out = open(output, 'a')
cov_max = 0
cov_target = 0

f_in = open(bracken_report, 'r')
species_max = ""

for line in f_in:
    fields = line.strip().split("\t")
    cov = float(fields[0].strip())
    lev = fields[3].strip()
    spe = fields[5].strip()
    if lev == "S":
        if cov > cov_max:
            cov_max = cov
            species_max = spe
        if spe == genus_species:
            cov_target = cov
f_in.close()
# Mix of different bacteria and there is no target:
if genus_species == "Genus species":
    cov_target = "-"

f_out.write("%s\t%s\t%s\t%s\n" % (sample_name, str(cov_target), species_max, str(cov_max)))
f_out.close()
