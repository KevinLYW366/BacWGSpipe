# Script to rename separated assembly fasta files based on viralVerify results

import sys
import os

work_dir = sys.argv[1]
sample = sys.argv[2]
#mob_recon_result = sys.argv[3]
viralverify_result = sys.argv[3]


def replace_fasta_contig_name(fasta, old_string, new_string):
    # Read in the file
    with open(fasta, 'r') as file:
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace(old_string, new_string)
    # Write the file out again
    with open(fasta, 'w') as file:
        file.write(filedata)


chr_n = 0
plas_n = 0
contig_n = 0
#for line in [line.strip() for line in open(mob_recon_result)][1:]:
for line in [line.strip() for line in open(viralverify_result)][1:]:
    #contig = line.split("\t")[4].split(" ")[0]
    contig = line.split(',')[0]
    #type = line.split("\t")[1]
    type = line.split(',')[1]
    if type == "Chromosome":
        chr_n += 1
        contig_fasta = os.path.join(work_dir, "{}.part_{}.fasta".format(sample, contig))
        new_contig_fasta = os.path.join(work_dir, "{}{}.fasta".format("chr", str(chr_n)))
        # modify each contig fasta filename
        cmd = "mv {} {}".format(contig_fasta, new_contig_fasta)
        os.system(cmd)
        # replace contig name in genome assembly fasta
        assembly_fasta = os.path.join(work_dir, "{}.fasta".format(sample))
        replace_fasta_contig_name(assembly_fasta, ">{} ".format(contig), ">{}{} ".format("chr", str(chr_n)))
        replace_fasta_contig_name(new_contig_fasta, ">{} ".format(contig), ">{}{} ".format("chr", str(chr_n)))
    # "Plasmid" = not "Chromosome"
    ## Note:
    ## 1. Viralverify might categorize short plasmid fragments into "Virus", "Unknown" or something else
    ##   other than "Chromosome" since not enough genes exist in these short fragments.
    ## 2. MOB-recon is a good replacement since it depends on searching database. However,
    ##   plasmid on the genome of rare bacteria strains might be classified as "Chromosome" if no match
    ##   could be found in database.
    else:
        plas_n += 1
        contig_fasta = os.path.join(work_dir, "{}.part_{}.fasta".format(sample, contig))
        new_contig_fasta = os.path.join(work_dir, "{}{}.fasta".format("plas", str(plas_n)))
        # modify each contig fasta filename
        cmd = "mv {} {}".format(contig_fasta, new_contig_fasta)
        os.system(cmd)
        # replace contig name in genome assembly fasta
        assembly_fasta = os.path.join(work_dir, "{}.fasta".format(sample))
        replace_fasta_contig_name(assembly_fasta, ">{} ".format(contig), ">{}{} ".format("plas", str(plas_n)))
        replace_fasta_contig_name(new_contig_fasta, ">{} ".format(contig), ">{}{} ".format("plas", str(plas_n)))
    # elif type == "plasmid":
    #     plas_n += 1
    #     contig_fasta = os.path.join(work_dir, "{}.part_{}.fasta".format(sample, contig))
    #     new_contig_fasta = os.path.join(work_dir, "{}{}.fasta".format("plas", str(plas_n)))
    #     # modify each contig fasta filename
    #     cmd = "mv {} {}".format(contig_fasta, new_contig_fasta)
    #     os.system(cmd)
    #     # replace contig name in genome assembly fasta
    #     assembly_fasta = os.path.join(work_dir, "{}.fasta".format(sample))
    #     replace_fasta_contig_name(assembly_fasta, ">{} ".format(contig), ">{}{} ".format("plas", str(plas_n)))
    #     replace_fasta_contig_name(new_contig_fasta, ">{} ".format(contig), ">{}{} ".format("plas", str(plas_n)))
    # else:
    #     contig_n += 1
    #     contig_fasta = os.path.join(work_dir, "{}.part_{}.fasta".format(sample, contig))
    #     new_contig_fasta = os.path.join(work_dir, "{}{}.fasta".format("contig", str(contig_n)))
    #     # modify each contig fasta filename
    #     cmd = "mv {} {}".format(contig_fasta, new_contig_fasta)
    #     os.system(cmd)
    #     # replace contig name in genome assembly fasta
    #     assembly_fasta = os.path.join(work_dir, "{}.fasta".format(sample))
    #     replace_fasta_contig_name(assembly_fasta, ">{} ".format(contig), ">{}{} ".format("contig", str(contig_n)))
    #     replace_fasta_contig_name(new_contig_fasta, ">{} ".format(contig), ">{}{} ".format("contig", str(contig_n)))

