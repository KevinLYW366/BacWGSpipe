# BacWGSpipe

A Snakemake workflow for a complete analysis of bacterial whole-genome sequencing data.

## Data Sources

BacWGSpipe supports three types of WGS data sources:
    
1. Illumina short reads 
2. Long reads (either PacBio or Nanopore) 
3. Illumina short reads + Long reads (either PacBio or Nanopore) from the same isolate
   - This hybrid mode enables us to get complete assembly sequences of both chromosome and plasmid.
	
## Workflow Structure

![Workflow](BacWGSpipe_workfflow.png)

## Workflow Usage

1. Clone the repository:

    ```bash
    git clone git@github.com:KevinLYW366/BacWGSpipe.git
    ```

2. Unzip some resource files:

    ```bash
    cd BacWGSpipe/resources
    unzip digIS-digISv1.2.zip
    unzip SPAdes-3.15.4-Linux.zip
    ```

3. Modify following items in Config file (config/config.yaml) based on project information every time before running the workflow.
   - 3.1 everything in "Data input" section;
     - Note: "seqdata_source" should match the data source you would input to the workflow. If the mode of Illumina short reads ONLY was selected, nothing in "Long reads (either PacBio or Nanopore)" needs to be modified, and vice versa.
   - 3.2 "threads" in "Analysis - Global" section (Threads used by tools for each sample);
   - 3.3 feel free to modify some analysis thresholds if you know what they mean, such as "plsdb_max_pvalue", "plasmidfinder_mincov" and ... in "Analysis" section.

4. Set up the snakemake environment:

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
5. Running the workflow:

    ```bash
    # Move to the directory of BacWGSpipe
    cd /path/to/BacWGSpipe
   
    # A dry-run is recommended at first to check if everything is okay (feel free to modify the number of CPU cores)
    snakemake --cores 64 --use-conda -r -p -n
   
    # If no error message shows up, let's do a formal run
    snakemake --cores 64 --use-conda -r -p
    ```

## Author

Yewei Lu (lyw@cwmda.com)

Xiangchen Li (lxc@cwmda.com)

