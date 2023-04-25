# pass parameters
args <- commandArgs(trailingOnly=TRUE)
rmd_script <- args[1]
output <- args[2]
sample_list <- args[3]
pass_sample_list <- args[4]
long_read_qc_raw <- args[5]
long_read_qc_clean <- args[6]
long_read_qc_dir <- args[7]
assembly_stats <- args[8]
genome_component_dir <- args[9]
plasmidfinder_results <- args[10]
phigaro_results <- args[11]
icefinder_results <- args[12]
digis_results <- args[13]
islandpath_results <- args[14]
cctyper_results <- args[15]
gene_function_dir <- args[16]
mlst_results <- args[17]

# absolute path for these files otherwise Rmarkdown will raise an error
sample_list <- paste0(normalizePath(dirname(sample_list)), '/',  basename(sample_list))
pass_sample_list <- paste0(normalizePath(dirname(pass_sample_list)), '/',  basename(pass_sample_list))
long_read_qc_raw <- paste0(normalizePath(dirname(long_read_qc_raw)), '/',  basename(long_read_qc_raw))
long_read_qc_clean <- paste0(normalizePath(dirname(long_read_qc_clean)), '/',  basename(long_read_qc_clean))
long_read_qc_dir <- paste0(normalizePath(dirname(long_read_qc_dir)), '/',  basename(long_read_qc_dir))
assembly_stats <- paste0(normalizePath(dirname(assembly_stats)), '/',  basename(assembly_stats))
genome_component_dir <- paste0(normalizePath(dirname(genome_component_dir)), '/',  basename(genome_component_dir))
plasmidfinder_results <- paste0(normalizePath(dirname(plasmidfinder_results)), '/',  basename(plasmidfinder_results))
phigaro_results <- paste0(normalizePath(dirname(phigaro_results)), '/',  basename(phigaro_results))
icefinder_results <- paste0(normalizePath(dirname(icefinder_results)), '/',  basename(icefinder_results))
digis_results <- paste0(normalizePath(dirname(digis_results)), '/',  basename(digis_results))
islandpath_results <- paste0(normalizePath(dirname(islandpath_results)), '/',  basename(islandpath_results))
cctyper_results <- paste0(normalizePath(dirname(cctyper_results)), '/',  basename(cctyper_results))
gene_function_dir <- paste0(normalizePath(dirname(gene_function_dir)), '/',  basename(gene_function_dir))
mlst_results <- paste0(normalizePath(dirname(mlst_results)), '/',  basename(mlst_results))

# run Rmarkdown script
rmarkdown::render(
  rmd_script,
  output_file = basename(output),
  output_dir = dirname(output),
  params=list(
    sample_list = sample_list,
    pass_sample_list = pass_sample_list,
    long_read_qc_raw = long_read_qc_raw,
    long_read_qc_clean = long_read_qc_clean,
    long_read_qc_dir = long_read_qc_dir,
    assembly_stats = assembly_stats,
    genome_component_dir = genome_component_dir,
    plasmidfinder_results = plasmidfinder_results,
    phigaro_results = phigaro_results,
    icefinder_results = icefinder_results,
    digis_results = digis_results,
    islandpath_results = islandpath_results,
    cctyper_results = cctyper_results,
    gene_function_dir = gene_function_dir,
    mlst_results = mlst_results
  )
)
