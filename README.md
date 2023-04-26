# BacWGSpipe

A Snakemake workflow for a complete analysis of bacterial whole-genome sequencing data.

## Data Sources

BacWGSpipe supports three types of WGS data sources:
    
1. Illumina short reads 
2. Long reads (either PacBio or Nanopore) 
3. Illumina short reads + Long reads (either PacBio or Nanopore) from the same isolate
   - This hybrid mode enables us to get complete assembly sequences of both chromosome and plasmid.
	
## Workflow Structure

<p align="center">
  <img src="https://github.com/KevinLYW366/BacWGSpipe/blob/main/BacWGSpipe_workflow.png" width="55%">
</p>

## Dependencies

1. Kraken2 database: 
   - Download PlusPF database from https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20230314.tar.gz and decompress
   - Modify "kraken2_db" in config/config.yaml
1. viralVerify HMM database:
   - Download from https://figshare.com/s/f897d463b31a35ad7bf0 and decompress
   - Modify "viralverify_hmm_path" in config/config.yaml
1. PLSDB Mash sketches file and meta archive file:
   - Download from https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plasmids_meta.tar.bz2 and https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2 
   - Decompress and put all files in the same directory
   - Modify "plsdb_mash_sketches" in config/config.yaml (path to plsdb.msh)
   - Modify "plsdb_meta_archive" in config/config.yaml (path to plsdb.tsv)
1. Platon database:
   - Download from https://zenodo.org/record/4066768/files/db.tar.gz and decompress
   - Modify "platon_db" in config/config.yaml
1. ICEfinder local version:
   - Visit https://bioinfo-mml.sjtu.edu.cn/ICEfinder/download.html to download
   - Modify "icefinder_dir" in config/config.yaml
1. EggNOG-mapper database:
   - Read "Installation - Setup" section in https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.10#user-content-Overview
   - Modify "eggnog_mapper_database" in config/config.yaml
1. COG database:
   - Download by `wget -r -e robots=off https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/`
   - Modify "cog_database_directory" in config/config.yaml
1. NCBI Taxonomy database:
   - Download from "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
   - Modify "taxonomy_database_directory" in config/config.yaml
1. CARD database:
   - Read "RGI Usage Documentation - RGI Databases" section in https://github.com/arpcard/rgi
   - Modify "card_database_directory" in config/config.yaml
1. signalp local version:
   - Visit https://services.healthtech.dtu.dk/services/SignalP-6.0/ to download
   - Create a local conda environment using workflow/env/tmhmm.yaml by `conda create -n signalp --file workflow/env/signalp.yaml`
   - Copy downloaded model files to the local conda environment by 
   ```bash
   conda activate signalp
   SIGNALP_DIR=$(python -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
   cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
   ```
   - Modify values of $envDIR in workflow/scripts/signalp_path.sh
1. TMHMM local version
   - Visit https://services.healthtech.dtu.dk/services/TMHMM-2.0/ to download
   - Create a local conda environment using workflow/env/tmhmm.yaml by `conda create -n tmhmm --file workflow/env/tmhmm.yaml`
   - Modify values of $envDIR and $binDIR in workflow/scripts/tmhmm_path.sh
1. Other dependencies:
   - Unzip resources/*.zip: `unzip resources/*.zip`

## Workflow Usage

1. Clone the repository:

    ```bash
    git clone git@github.com:KevinLYW366/BacWGSpipe.git
    ```

2. Modify following items in Config file (config/config.yaml) based on project information every time before running the workflow.
   - 3.1 everything in "Data input" section;
     - Note: "seqdata_source" should match the data source you would input to the workflow. If the mode of Illumina short reads ONLY was selected, nothing in "Long reads" needs to be modified, and vice versa.
   - 3.2 "threads" in "Analysis - Global" section (Threads used by tools for each sample);
   - 3.3 feel free to modify some analysis thresholds if you know what they mean, such as "plsdb_max_pvalue", "plasmidfinder_mincov" and ... in "Analysis" section.

3. Set up the snakemake environment:

    ```bash
    # Add conda channels
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
      
    # Install mamba first (mamba provides a faster and more roboust way for conda packages installation)
    conda install -n base -c conda-forge mamba
      
    # Install snakemake using mamba
    mamba create -c conda-forge -c bioconda -n snakemake snakemake
    ```
4. Run the workflow:

    ```bash
    # Move to the directory of BacWGSpipe
    cd /path/to/BacWGSpipe
   
    # A dry-run is recommended at first to check if everything is okay (feel free to modify the number of CPU cores)
    snakemake --cores 64 --use-conda -r -p -n
   
    # If no error message shows up, let's do a formal run
    snakemake --cores 64 --use-conda -r -p
    ```

# Test dataset
V. Murigneux et al., “MicroPIPE: Validating an end-to-end workflow for high-quality complete bacterial genome construction,” BMC genomics, vol. 22, no. 1, p. 474, 2021.

- 12 ST131 Escherichia coli strains with both Nanopore long-read sequencing data (SRA accession: SRP293329) and Illumina short-read sequencing data (SRA accession: ERP001354).

## Author

Yewei Lu (lyw@cwmda.com)

Xiangchen Li (lxc@cwmda.com)

