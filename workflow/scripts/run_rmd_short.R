# pass parameters
args <- commandArgs(trailingOnly=TRUE)
rmd_script <- args[1]
output <- args[2]
sample_list <- args[3]
pass_sample_list <- args[4]
fastp_qc_stats <- args[5]
kraken2_results <- args[6]
assembly_stats <- args[7]
genome_component_dir <- args[8]
plasmidfinder_results <- args[9]
platon_results <- args[10]
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
fastp_qc_stats <- paste0(normalizePath(dirname(fastp_qc_stats)), '/',  basename(fastp_qc_stats))
kraken2_results <- paste0(normalizePath(dirname(kraken2_results)), '/',  basename(kraken2_results))
assembly_stats <- paste0(normalizePath(dirname(assembly_stats)), '/',  basename(assembly_stats))
genome_component_dir <- paste0(normalizePath(dirname(genome_component_dir)), '/',  basename(genome_component_dir))
plasmidfinder_results <- paste0(normalizePath(dirname(plasmidfinder_results)), '/',  basename(plasmidfinder_results))
platon_results <- paste0(normalizePath(dirname(platon_results)), '/',  basename(platon_results))
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
    fastp_qc_stats = fastp_qc_stats,
    kraken2_results = kraken2_results,
    assembly_stats = assembly_stats,
    genome_component_dir = genome_component_dir,
    plasmidfinder_results = plasmidfinder_results,
    platon_results = platon_results,
    phigaro_results = phigaro_results,
    icefinder_results = icefinder_results,
    digis_results = digis_results,
    islandpath_results = islandpath_results,
    cctyper_results = cctyper_results,
    gene_function_dir = gene_function_dir,
    mlst_results = mlst_results
  )
)
