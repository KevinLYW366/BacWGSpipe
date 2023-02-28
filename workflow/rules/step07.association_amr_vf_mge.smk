###############################################
# Step 07. Association between AMR/VF and MGE #
###############################################

# Match MGE results into AMR/VF matrix
## Place MGE position after AMR/VF position
rule match_mge:
    input:
        ## AMR/VF results matrix
        amr_matrix = "results/04.gene_function/all_card_results_matrix.xls",
        vf_matrix = "results/04.gene_function/all_vfdb_results_matrix.xls",
        ## merged MGE results for each MGE category
        phigaro = "results/03.genome_component/all_phigaro_results.xls",
        icefinder = "results/03.genome_component/all_icefinder_results.xls",
        #plasmidfinder = "results/08.plasmid/plasmidfinder/all_plasmidfinder_result.txt",
        #platon = "results/08.plasmid/platon/all_platon_result.txt",
        digis = "results/03.genome_component/all_digis_results.xls",
        islandpath = "results/03.genome_component/all_islandpath_results.xls"
    output:
        ## AMR/VF results matrix with matched MGE results
        amr_mge_matrix = "results/07.association_amr_vf_mge/all_amr_mge_matrix.xls",
        vf_mge_matrix = "results/07.association_amr_vf_mge/all_vf_mge_matrix.xls"
    params:
        vfdb_result_dir = "results/04.gene_function/",
        prokka_result_dir = "results/03.genome_component",
        sample_list = config["sample_list"],
        script = config["script_match_mge"],
        distance = config["distance_amr_vf_mge"]
    log:
        "logs/07.association_amr_vf_mge/match_mge.log"
    shell:
        """
        python {params.script} {input.amr_matrix} {input.vf_matrix} {output.amr_mge_matrix} \
        {output.vf_mge_matrix} {input.phigaro} {input.icefinder} {input.digis} {input.islandpath} \
         {params.distance} > {log} 2>&1
        """