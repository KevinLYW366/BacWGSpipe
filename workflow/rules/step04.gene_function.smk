###################################
# Step 04. Gene function analysis #
###################################

##### 1. Gene function annotation #####
# 1.1 EggNOG-mapper: Fast genome-wide functional annotation through orthology assignment
rule eggnog_mapper:
    input:
        # prokka annotated protein FASTA file of the translated CDS sequences
        "results/03.genome_component/{sample}/prokka/{sample}.faa"
    output:
        "results/04.gene_function/{sample}/eggnog-mapper/{sample}.emapper.annotations"
    params:
        # eggNOG-mapper databases
        db = config["eggnog_mapper_database"],
        outdir = "results/04.gene_function/{sample}/eggnog-mapper"
    log:
        "logs/04.gene_function/eggnog-mapper/{sample}.log"
    threads:
        config["threads"]
    conda:
        "../envs/eggnog-mapper.yaml"
    shell:
        """
        emapper.py --cpu {threads} --override -i {input} --itype proteins \
        --data_dir {params.db} -o {wildcards.sample} --excel --go_evidence all \
        --output_dir {params.outdir} > {log} 2>&1
        """

# 1.2 GO: Gene Ontology database annotation
## Perform GO annotation on eggNOG-mapper results
rule go_annotate:
    input:
        "results/04.gene_function/{sample}/eggnog-mapper/{sample}.emapper.annotations"
    output:
        go = "results/04.gene_function/{sample}/GO/{sample}_go.xls",
        go_level = "results/04.gene_function/{sample}/GO/{sample}_go_with_level.xls"
    params:
        # GO database core GO ontology (OBO Format) file
        obo = config["go_obo"],
        outdir = "results/04.gene_function/{sample}/GO",
        # Python script to annotate GO terms
        script = config["go_annotate_script"]
    log:
        "logs/04.gene_function/GO/go_annotate/{sample}.log"
    conda:
        "../envs/goatools.yaml"
    shell:
        """
        python {params.script} {params.obo} {input} {params.outdir} {wildcards.sample} > {log} 2>&1
        """

## Plot GO annotation results
rule go_plot:
    input:
        "results/04.gene_function/{sample}/GO/{sample}_go_with_level.xls"
    output:
        pdf = "results/04.gene_function/{sample}/GO/{sample}_go.pdf",
        png = "results/04.gene_function/{sample}/GO/{sample}_go.png"
    params:
        # R script to plot GO annotation results
        script = config["go_plot_script"],
        outdir = "results/04.gene_function/{sample}/GO"
    log:
        "logs/04.gene_function/GO/go_plot/{sample}.log"
    conda:
        "../envs/R_plot.yaml"
    shell:
        """
        Rscript {params.script} {input} {wildcards.sample} {params.outdir} > {log} 2>&1
        """

# 1.3 KEGG: KEGG PATHWAY database
## Perform KEGG PATHWAY database annotation and plotting on eggNOG-mapper results
rule kegg_annotate_plot:
    input:
        "results/04.gene_function/{sample}/eggnog-mapper/{sample}.emapper.annotations"
    output:
        kegg_annotation = "results/04.gene_function/{sample}/KEGG/{sample}_kegg.xls",
        kegg_plot_pdf = "results/04.gene_function/{sample}/KEGG/{sample}_kegg.pdf",
        kegg_plot_png = "results/04.gene_function/{sample}/KEGG/{sample}_kegg.png",
    params:
        # R script to annotate and plot KEGG PATHWAY info
        script = config["kegg_annotate_plot_script"],
        # Path to KEGG PATHWAY database file
        db = config["kegg_db"],
        outdir = "results/04.gene_function/{sample}/KEGG"
    log:
        "logs/04.gene_function/KEGG/{sample}.log"
    conda:
        "../envs/R_plot.yaml"
    shell:
        """
        Rscript {params.script} {input} {params.db} {wildcards.sample} {params.outdir} > {log} 2>&1
        """

# 1.4 COG: NCBI Clusters of Orthologous Genes database
## Perform COG annotation on eggNOG-mapper results
rule cog_annotate:
    input:
        "results/04.gene_function/{sample}/eggnog-mapper/{sample}.emapper.annotations"
    output:
        "results/04.gene_function/{sample}/COG/{sample}_cog.xls"
    params:
        # Python script to annotate COG terms
        script = config["cog_annotate_script"],
        # Path to COG database files directory
        db = config["cog_database_directory"],
        outdir = "results/04.gene_function/{sample}/COG"
    log:
        "logs/04.gene_function/COG/cog_annotate/{sample}.log"
    shell:
        """
        python {params.script} {params.db} {input} {params.outdir} {wildcards.sample} > {log} 2>&1
        """

## Plot COG annotation results
rule cog_plot:
    input:
        "results/04.gene_function/{sample}/COG/{sample}_cog.xls"
    output:
        pdf = "results/04.gene_function/{sample}/COG/{sample}_cog.pdf",
        png = "results/04.gene_function/{sample}/COG/{sample}_cog.png"
    params:
        # R script to plot COG annotation results
        script = config["cog_plot_script"],
        # Path to COG database files directory
        db = config["cog_database_directory"],
        outdir = "results/04.gene_function/{sample}/COG"
    log:
        "logs/04.gene_function/COG/cog_plot/{sample}.log"
    conda:
        "../envs/R_plot.yaml"
    shell:
        """
        Rscript {params.script} {input} {params.db} {wildcards.sample} {params.outdir} > {log} 2>&1
        """

# 1.5 NCBI Taxonomy
## Upgrade the local NCBI Taxonomy database in ETE toolkit
rule taxonomy_db_upgrade:
    input:
        expand("results/04.gene_function/{sample}/eggnog-mapper/{sample}.emapper.annotations", sample=SAMPLES)
    output:
        "results/workflow_signal/04.gene_function/taxonomy_db_upgrade.done"
    params:
        # Python script to annotate NCBI Taxonomy terms
        script = config["taxonomy_annotate_script"],
        # Path to NCBI Taxonomy database files directory
        db = config["taxonomy_database_directory"],
        outdir = "results/workflow_signal/04.gene_function"
    log:
        "logs/04.gene_function/Taxonomy/taxonomy_db_upgrade.log"
    conda:
        "../envs/ete3.yaml"
    shell:
        """
        python {params.script} {params.db} empty {params.outdir} empty y > {log} 2>&1
        """


## Perform NCBI Taxonomy annotation on eggNOG-mapper results
rule taxonomy_annotate:
    input:
        emapper = "results/04.gene_function/{sample}/eggnog-mapper/{sample}.emapper.annotations",
        flag = "results/workflow_signal/04.gene_function/taxonomy_db_upgrade.done"
    output:
        "results/04.gene_function/{sample}/Taxonomy/{sample}_taxonomy.xls"
    params:
        # Python script to annotate NCBI Taxonomy terms
        script = config["taxonomy_annotate_script"],
        # Path to NCBI Taxonomy database files directory
        db = config["taxonomy_database_directory"],
        outdir = "results/04.gene_function/{sample}/Taxonomy"
    log:
        "logs/04.gene_function/Taxonomy/taxonomy_annotate/{sample}.log"
    conda:
        "../envs/ete3.yaml"
    shell:
        """
        python {params.script} {params.db} {input.emapper} {params.outdir} {wildcards.sample} n > {log} 2>&1
        """

## Plot NCBI Taxonomy annotation results
rule taxonomy_plot:
    input:
        "results/04.gene_function/{sample}/Taxonomy/{sample}_taxonomy.xls"
    output:
        pdf = "results/04.gene_function/{sample}/Taxonomy/{sample}_taxonomy.pdf",
        png = "results/04.gene_function/{sample}/Taxonomy/{sample}_taxonomy.png"
    params:
        # R script to plot NCBI Taxonomy annotation results
        script = config["taxonomy_plot_script"],
        # Path to NCBI Taxonomy database files directory
        db = config["taxonomy_database_directory"],
        outdir = "results/04.gene_function/{sample}/Taxonomy"
    log:
        "logs/04.gene_function/Taxonomy/taxonomy_plot/{sample}.log"
    conda:
        "../envs/R_plot.yaml"
    shell:
        """
        Rscript {params.script} {input} {wildcards.sample} {params.outdir} > {log} 2>&1
        """

##### 2. Virulence or pathogenicity analysis #####
# 2.1 Antimicrobial resistance (AMR) genes
## 2.1.1 RGI: predict resistome with reference data from CARD database
### RGI is a tool to predict resistome(s) from protein or nucleotide data based on homology and SNP models
### RGI uses reference data from the Comprehensive Antibiotic Resistance Database (CARD)
### Run RGI for each sample
rule rgi_run:
    input:
        # protein FASTA file of the translated CDS sequences
        faa = "results/03.genome_component/{sample}/prokka/{sample}.faa"
    output:
        txt = "results/04.gene_function/{sample}/CARD/rgi/{sample}.report.txt"
    params:
        input_dir = lambda wildcards, input: extract_absolute_dirname(wildcards, input.faa),
        output_dir = lambda wildcards, output: extract_absolute_dirname(wildcards, output.txt),
        card_dir = config['card_database_directory'],
        workdir = WORKDIR
    threads:
        config["threads"]
    log:
        "logs/04.gene_function/CARD/rgi_run/{sample}.log"
    conda:
        "../envs/rgi.yaml"
    shell:
        """
        # have to run RGI local mode in CARD database directory
        cd {params.card_dir}
        
        # run RGI
        rgi main --input_sequence {params.input_dir}/{wildcards.sample}.faa \
        --output_file {params.output_dir}/{wildcards.sample}.report \
        --input_type protein --local --clean > {params.workdir}/{log} 2>&1
        
        # go back to snakemake workflow directory
        cd {params.workdir}
        """

### Merge RGI results for all samples
rule rgi_merge:
    input:
        expand("results/04.gene_function/{sample}/CARD/rgi/{sample}.report.txt", sample=SAMPLES)
    output:
        "results/04.gene_function/all_rgi_results.xls"
    params:
        input_dir = "results/04.gene_function",
        output_dir = "results/04.gene_function",
        sample_list = config["sample_list"],
        script = config["merge_results_script"]
    log:
        "logs/04.gene_function/CARD/rgi_merge.log"
    shell:
        """
        python {params.script} -m rgi -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
        """

## 2.1.2 Abricate (CARD): Mass screening of contigs for antimicrobial resistance or virulence genes
### Note:
### Abricate only detects acquired resistance genes, NOT point mutations.
### Abricate uses a DNA sequence database, not protein.
### Since we use RGI to seach AMR genes in protein sequence, we might lose some AMR genes because of poor assembly
### quality or bad gene prediction results. Abricate provides methods to catch up on the missing AMR genes.
### Run Abricate (CARD) for each sample.
rule abricate_card_run:
    input:
        # Nucleotide FASTA file of the assembly contig sequences
        fna = "results/03.genome_component/{sample}/prokka/{sample}.fna"
    output:
        "results/04.gene_function/{sample}/CARD/abricate/{sample}_abricate_card_results.xls"
    params:
        ### Python script to reformat Abricate (CARD) results
        script = config["script_reformat_abricate_card"],
        card_aro_index = os.path.join(config["card_database_directory"],"card_database/aro_index.tsv")
    threads:
        config["threads"]
    log:
        "logs/04.gene_function/CARD/abricate_card_run/{sample}.log"
    conda:
        "../envs/abricate.yaml"
    shell:
        """
        abricate --db card --threads {threads} {input.fna} > {output}.tmp 2> {log}
        python {params.script} {output}.tmp {params.card_aro_index} {output} >> {log} 2>&1
        """

### Merge Abricate (CARD) results for all samples
rule abricate_card_merge:
    input:
        expand("results/04.gene_function/{sample}/CARD/abricate/{sample}_abricate_card_results.xls", sample=SAMPLES)
    output:
        "results/04.gene_function/all_abricate_card_results.xls"
    params:
        input_dir = "results/04.gene_function",
        output_dir = "results/04.gene_function",
        sample_list = config["sample_list"],
        script = config["merge_results_script"]
    log:
        "logs/04.gene_function/CARD/abricate_card_merge.log"
    shell:
        """
        python {params.script} -m abricate_card -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
        """

## 2.1.3 Merge RGI and Abricate (CARD) results
rule rgi_abricate_card_merge:
    input:
        rgi = "results/04.gene_function/all_rgi_results.xls",
        abricate_card = "results/04.gene_function/all_abricate_card_results.xls"
    output:
        "results/04.gene_function/all_card_results.xls"
    params:
        script = config["script_merge_rgi_abricate_card"],
        card_aro_index = os.path.join(config["card_database_directory"], "card_database/aro_index.tsv")
    log:
        "logs/04.gene_function/CARD/rgi_abricate_card_merge.log"
    shell:
        """
        python {params.script} {input.rgi} {input.abricate_card} {output} {params.card_aro_index} > {log} 2>&1
        """

## 2.1.4 Extract AMR results of all samples into a matrix with a specific format
### Output format:
###          Best_Hit_ARO_1    Best_Hit_ARO_2  ...
### Sample1  Contig:Start:End  ...
### Sample2  ...
### ...
rule extract_card_matrix:
    input:
        ## merged AMR prediction report generated by RGI and Abricate (CARD)
        "results/04.gene_function/all_card_results.xls",
    output:
        matrix = "results/04.gene_function/all_card_results_matrix.xls",
        aro_info = "results/04.gene_function/card_aro_info.xls",
        binary_matrix = "results/04.gene_function/all_card_results_matrix_binary.xlsx",
    params:
        card_aro_index = os.path.join(config["card_database_directory"], "card_database/aro_index.tsv"),
        script = config["script_extract_card"]
    log:
        "logs/04.gene_function/CARD/extract_card_matrix.log"
    conda:
        "../envs/xlsxwriter.yaml"
    shell:
        """
        python {params.script} {input} {params.card_aro_index} {output.matrix} {output.aro_info} {output.binary_matrix} > {log} 2>&1
        """

# 2.2 Virulence factors (VF)
## Blast + VFDB: Search and annotate virulence genes
### The virulence factor database (VFDB) is an integrated and comprehensive online resource
###   for curating information about virulence factors of bacterial pathogens
rule vfdb_annotate:
    input:
        # genes predicted by Prokka (nucleotide seq)
        "results/03.genome_component/{sample}/prokka/{sample}.ffn"
    output:
        # with predicted gene sequences
        "results/04.gene_function/{sample}/VFDB/{sample}_vfdb_blastn_onGenes.txt"
    params:
        # script to run blast related work (from https://github.com/fmalmeida/bacannot)
        script = config["blast_script"],
        # VFDB full nucleotide database
        db = config["vfdb_database_directory"],
        # blast parameters
        ## VFDB minimum coverage
        minid = config["vfdb_blast_minid"],
        ## VFDB minimum hreshold for identity
        mincov = config["vfdb_blast_mincov"]
    log:
        "logs/04.gene_function/VFDB/vfdb_annotate/{sample}.log"
    threads:
        config["threads"]
    conda:
        "../envs/blast.yaml"
    shell:
        """
        # With predicted gene sequences
        python {params.script} blastn --query {input} --db {params.db} --minid {params.minid} \
        --mincov {params.mincov} --threads {threads} --out {output} --2way > {log} 2>&1
        """

## merge VFDB annotation results for all samples
rule vfdb_merge:
    input:
        expand("results/04.gene_function/{sample}/VFDB/{sample}_vfdb_blastn_onGenes.txt", sample=SAMPLES)
    output:
        "results/04.gene_function/all_vfdb_results.xls"
    params:
        input_dir = "results/04.gene_function",
        output_dir = "results/04.gene_function",
        sample_list = config["sample_list"],
        script = config["merge_results_script"]
    log:
        "logs/04.gene_function/VFDB/vfdb_merge.log"
    shell:
        """
        python {params.script} -m vfdb -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
        """

## Extract VF results of all samples into a matrix with a specific format
### Output format:
###          VF1               VF2  ...
### Sample1  Contig:Start:End  ...
### Sample2  ...
### ...
rule extract_vfdb_matrix:
    input:
        ## merged AMR prediction report generated by RGI
        "results/04.gene_function/all_vfdb_results.xls",
    output:
        matrix = "results/04.gene_function/all_vfdb_results_matrix.xls",
        vf_info = "results/04.gene_function/vfdb_vf_info.xls",
        binary_matrix = "results/04.gene_function/all_vfdb_results_matrix_binary.xlsx",
    params:
        script = config["script_extract_vfdb"]
    log:
        "logs/04.gene_function/VFDB/extract_vfdb_matrix.log"
    conda:
        "../envs/xlsxwriter.yaml"
    shell:
        """
        python {params.script} {input} {output.matrix} {output.vf_info} {output.binary_matrix} > {log} 2>&1
        """

# 2.3 Secretory proteins
## SignalP: predict secretory signal peptides in protein sequences
rule signalp:
    input:
        # protein FASTA file of the translated CDS sequences
        "results/03.genome_component/{sample}/prokka/{sample}.faa"
    output:
        # with predicted gene sequences
        "results/04.gene_function/{sample}/secretory_protein/signalp/prediction_results.txt"
    params:
        # script to load a local conda environment settings
        env = config["signalp_env_script"],
        outdir = "results/04.gene_function/{sample}/secretory_protein/signalp"
    log:
        "logs/04.gene_function/secretory_protein/signalp/{sample}.log"
    threads:
        config["threads"]
    shell:
        """
        # load local conda environment
        source {params.env}
        
        # run SignalP
        signalp6 --fastafile {input} --organism other --output_dir {params.outdir} \
        --format none --mode fast -wp {threads} > {log} 2>&1
        """

## TMHMM: predict transmembrane helices
rule tmhmm:
    input:
        # protein FASTA file of the translated CDS sequences
        "results/03.genome_component/{sample}/prokka/{sample}.faa"
    output:
        # with predicted gene sequences
        "results/04.gene_function/{sample}/secretory_protein/tmhmm/prediction_results.txt"
    params:
        # script to load a local conda environment settings
        env = config["tmhmm_env_script"],
        outdir = "results/04.gene_function/{sample}/secretory_protein/tmhmm"
    log:
        "logs/04.gene_function/secretory_protein/tmhmm/{sample}.log"
    shell:
        """
        # load local conda environment
        source {params.env}
    
        # run SignalP
        tmhmm -workdir {params.outdir} -short {input} > {output} 2> {log}
        
        # remove tmp file directory
        rm -rf {params.outdir}/TMHMM* 2>> {log}
        """

## Secretory proteins:
##   1. with signal peptides predicted in SignalP results
##   2. Topology=o (without transmembrane helices and outside the membrane) in TMHMM results
rule sec_pro:
    input:
        # SignalP results
        signalp = "results/04.gene_function/{sample}/secretory_protein/signalp/prediction_results.txt",
        # TMHMM results
        tmhmm = "results/04.gene_function/{sample}/secretory_protein/tmhmm/prediction_results.txt"
    output:
        # Secretory proteins identification results
        "results/04.gene_function/{sample}/secretory_protein/{sample}_secretory_protein_results.xls"
    params:
        outdir = "results/04.gene_function/{sample}/secretory_protein",
        # Python script to identify secretory proteins based on SignalP and TMHMM results
        script = config["script_identify_sec_pro"]
    log:
        "logs/04.gene_function/secretory_protein/sec_pro/{sample}.log"
    shell:
        """
        python {params.script} {input.signalp} {input.tmhmm} {output} > {log} 2>&1
        """
