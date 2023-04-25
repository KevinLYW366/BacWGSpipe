##### Common Module #####

#######################
# Load python library #
#######################
from snakemake.utils import min_version
import pandas as pd
import os

###########################
# Check snakemake version #
###########################
min_version("6.13.1")

##############################
# Identify working directory #
##############################
WORKDIR = os.getcwd()

####################
# Read config file #
####################
# Load configs
configfile: "config/config.yaml"

################
# Load options #
################
# Load a sample list from sample_list.txt
with open(config["sample_list"], 'r') as f:
    SAMPLES = f.read().strip().split('\n')

# If kraken filter is turned on, an updated sample list would be generated and used in some analysis steps
if config["kraken2_filter"]:
    SAMPLE_LIST_UPDATE = "results/01.data_clean/short_read/all.kraken2.pass.list"
else:
    SAMPLE_LIST_UPDATE = config["sample_list"]

# Input fastq file name options
## Illumina short reads
if config["illumina_fastq_read_id_format"]:
    FASTQREAD2 = "R"
else:
    FASTQREAD2 = ""
if config["illumina_fastq_suffix_format"]:
    FASTQSUFFIX2 = "fq.gz"
else:
    FASTQSUFFIX2 = "fastq.gz"
# Long reads (either PacBio or Nanopore)
if config["long_fastq_suffix_format"]:
    FASTQSUFFIX3 = "fq.gz"
else:
    FASTQSUFFIX3 = "fastq.gz"

# Input data directory format options
# Illumina short reads
if config["illumina_data_dir_format"]:
    DATADIRINDEX2 = "/{sample}"
else:
    DATADIRINDEX2 = ""
# Long reads (either PacBio or Nanopore)
if config["long_data_dir_format"]:
    DATADIRINDEX3 = "/{sample}"
else:
    DATADIRINDEX3 = ""

# Assembly Strategy
if config["seqdata_source"] == 0:
    ASSEMBLY_STRATEGY = "hybrid"
elif config["seqdata_source"] == 2:
    ASSEMBLY_STRATEGY = "short"
elif config["seqdata_source"] == 1:
    ASSEMBLY_STRATEGY = "long"

####################
# Helper functions #
####################
def extract_absolute_dirname(wildcards, file_path):
    """
    Extract ABSOLUTE directory path from a file path.
    Usage example:
    Directory path extracted from input or output section could be used in params section.
    """
    absolute_dir_name = os.path.dirname(os.path.abspath(file_path))
    return absolute_dir_name


def get_qualified_results(path, wildcards):
    """
    Input a file path with sample wildcards such as "results/04.AMR/rgi/{sample}/{sample}.report.txt".
    Return required file path of qualified samples based on kraken qc result.
    This function will be used to generate an output file list, which is helpful to
      prevent unqualified samples from running specific analysis steps of the workflow.
    """
    thres = config["kraken2_threshold"]
    if config["kraken2_filter"]:
        qc = pd.read_csv(checkpoints.kraken2_qc.get().output[0], sep="\t", header=None, comment='#')
        qc.columns = ["sample", "target_percent", "species", "species_percent"]
        return expand(path, sample=qc[qc["target_percent"] > thres]["sample"])
    else:
        return expand(path, sample=SAMPLES)

