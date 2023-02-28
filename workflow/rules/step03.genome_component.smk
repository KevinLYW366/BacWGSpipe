#############################################################
# Step 03. Genome annotation and Genomic component analysis #
#############################################################

##### 1. Genes & ncRNA #####
# Prokka: Annotate assembly genome using
## Genes & ncRNA INFO could be found in prokka annotation results
rule prokka:
    input:
        "results/02.assembly/{sample}/assembly/{sample}.fasta"
    output:
        gff = "results/03.genome_component/{sample}/prokka/{sample}.gff",
        faa = "results/03.genome_component/{sample}/prokka/{sample}.faa",
        fna = "results/03.genome_component/{sample}/prokka/{sample}.fna",
        ffn = "results/03.genome_component/{sample}/prokka/{sample}.ffn",
        gbk = "results/03.genome_component/{sample}/prokka/{sample}.gbk",
        tsv = "results/03.genome_component/{sample}/prokka/{sample}.tsv"
    params:
        outdir = "results/03.genome_component/{sample}/prokka",
        prefix = "{sample}",
        genus = config["genus"],
        species = config["species"],
        strain = "{sample}",
        # Force overwriting existing output folder
        force = "--force"
    log:
        "logs/03.genome_component/prokka/{sample}.log"
    threads:
        config["threads"]
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        # Format date information in gbk file in English
        export LANG=en_US.UTF-8
        
        # Run prokka
        prokka --cpus {threads} --outdir {params.outdir} --locustag {params.strain} --strain {params.strain} \
        --genus {params.genus} --species {params.species} --prefix {params.prefix} {params.force} {input} > {log} 2>&1
        """

## Plot gene length distribution based on prokka prediction results
rule plot_gene_length:
    input:
        "results/03.genome_component/{sample}/prokka/{sample}.tsv"
    output:
        "results/03.genome_component/{sample}/prokka/{sample}_gene_length.png"
    params:
        outdir = "results/03.genome_component/{sample}/prokka",
        # R script to plot gene length distribution
        script = config["script_plot_gene_length"]
    log:
        "logs/03.genome_component/prokka/plot_gene_length/{sample}.log"
    conda:
        # Conda env R ggplot2
        "../envs/R_plot.yaml"
    shell:
        """
        Rscript {params.script} {input} {params.outdir} {wildcards.sample} > {log} 2>&1
        """

##### 2. Mobile genetic elements (MGEs) #####
# 2.1 Plasmid
if ASSEMBLY_STRATEGY == "hybrid":
    ## 2.1.1 PLSDB annotation
    ### Run PLSDB annotation using Mash screen (Which plasmids are contained in the input sample?)
    rule plsdb_screen:
        input:
            "results/02.assembly/{sample}/assembly/separate_assembly.done"
        output:
            "results/03.genome_component/{sample}/plsdb/plsdb_screen.done"
        params:
            outdir = "results/03.genome_component/{sample}/plsdb",
            # PLSDB Mash sketches file
            db = config["plsdb_mash_sketches"],
            # Maximal p-value to report in Mash screen
            max_pvalue = config["plsdb_max_pvalue"],
            # Minimal identity in Mash screen in Mash screen
            min_ident = config["plsdb_min_ident"],
        log:
            "logs/03.genome_component/plsdb/plsdb_screen/{sample}.log"
        threads:
            config["threads"]
        conda:
            "../envs/plsdb.yaml"
        shell:
            """
            # Screen PLSDB on each plasmid fasta of each sample
            shopt -s nullglob
            for f in results/02.assembly/{wildcards.sample}/assembly/plas*.fasta
            do
              plas_num=`basename $f .fasta`
              
              # Run mash screen
              mash screen {params.db} $f -v {params.max_pvalue} -i {params.min_ident} \
              -p {threads} > {params.outdir}/{wildcards.sample}_${{plas_num}}_plsdb_result.xls.tmp 2> {log}
              
            done
            
            # Create a flag file when everything done
            touch {output} 2>> {log}
            """

    ### Extract PLSDB meta info
    rule plsdb_meta:
        input:
            "results/03.genome_component/{sample}/plsdb/plsdb_screen.done"
        output:
            "results/03.genome_component/{sample}/plsdb/plsdb_meta.done"
        params:
            outdir = "results/03.genome_component/{sample}/plsdb",
            # Python script to extract matched plasmid info from PLSDB meta archive
            script = config["plsdb_script_meta"],
            # Python script to extract matched plasmid info from PLSDB meta archive
            meta = config["plsdb_meta_archive"]
        log:
            "logs/03.genome_component/plsdb/plsdb_meta/{sample}.log"
        shell:
            """
            shopt -s nullglob
            for f in results/02.assembly/{wildcards.sample}/assembly/plas*.fasta
            do
              plas_num=`basename $f .fasta`
              
              python {params.script} {params.outdir}/{wildcards.sample}_${{plas_num}}_plsdb_result.xls.tmp \
              {params.meta} {params.outdir}/{wildcards.sample}_${{plas_num}}_plsdb_result.xls.tmp2 $f > {log} 2>&1
    
              head -n21 {params.outdir}/{wildcards.sample}_${{plas_num}}_plsdb_result.xls.tmp2 > \
              {params.outdir}/{wildcards.sample}_${{plas_num}}_plsdb_result.xls 2>> {log} 
            done
            
            # Clean up temporary files
            rm -rf {params.outdir}/{wildcards.sample}_*_plsdb_result.xls.tmp* 2>> {log}
            
            # Create a flag file when everything done
            touch {output} 2>> {log}
            """

    ### Merge PLSDB annotation results for all samples
    rule plsdb_merge:
        input:
            expand("results/03.genome_component/{sample}/plsdb/plsdb_meta.done", sample=SAMPLES)
        output:
            "results/03.genome_component/all_plsdb_results.xls"
        params:
            indir = "results/03.genome_component",
            outdir = "results/03.genome_component",
            sample_list = config["sample_list"],
            script = config["merge_results_script"]
        log:
            "logs/03.genome_component/plsdb/plsdb_merge.log"
        shell:
            """
            python {params.script} -m plsdb -i {params.indir} -o {params.outdir} -s {params.sample_list} > {log} 2>&1
            """

    ## 2.1.2 Mobility prediction (MOB-suite)
    ### Run MOB-suite
    rule mob_suite_run:
        input:
            "results/02.assembly/{sample}/assembly/separate_assembly.done"
        output:
            flag = "results/03.genome_component/{sample}/mob-suite/mob_suite_run.done",
            result = "results/03.genome_component/{sample}/mob-suite/{sample}_mobtyper_result.xls"
        params:
            outdir = "results/03.genome_component/{sample}/mob-suite",
            # MOB-suite database
            db = config["mob_suite_database"]
        log:
            "logs/03.genome_component/mob-suite/mob_suite_run/{sample}.log"
        threads:
            config["threads"]
        conda:
            "../envs/mob-suite.yaml"
        shell:
            """
            # Run MOB-suite on each plasmid fasta of each sample
            shopt -s nullglob
            for f in results/02.assembly/{wildcards.sample}/assembly/plas*.fasta
            do
              plas_num=`basename $f .fasta`
              mob_typer --infile $f --out_file {params.outdir}/{wildcards.sample}_${{plas_num}}_mobtyper_result.xls \
              -d {params.db} --num_threads {threads} >> {log} 2>&1
            done
            
            # Merge all MOB-suite results for each sample
            n=0
            touch {output.result} 2>> {log}
            shopt -s nullglob
            for f in {params.outdir}/{wildcards.sample}_plas*_mobtyper_result.xls
            do
              if [ $n -eq 0 ]
              then
                cat $f >> {output.result} 2>> {log}
                rm -rf $f
                n=1
              else
                sed '1d' $f >> {output.result} 2>> {log}
                rm -rf $f
              fi
            done
            
            # If no pladmid is identified 
            if [ $n -eq 0 ]
              then
                echo "sample_id\tnum_contigs\tsize\tgc\tmd5\trep_type(s)\trep_type_accession(s)\trelaxase_type(s)\trelaxase_type_accession(s)\tmpf_type	mpf_type_accession(s)\torit_type(s)\torit_accession(s)\tpredicted_mobility\tmash_nearest_neighbor\tmash_neighbor_distance\tmash_neighbor_identification\tprimary_cluster_id\tsecondary_cluster_id\tpredicted_host_range_overall_rank\tpredicted_host_range_overall_name	observed_host_range_ncbi_rank\tobserved_host_range_ncbi_name\treported_host_range_lit_rank\treported_host_range_lit_name\tassociated_pmid(s)" >> {output.result} 2>> {log}
            fi
            
            # Create flag file
            touch {output.flag} 2>> {log}
            """

    ### Merge MOB-suite results for all samples
    rule mob_suite_merge:
        input:
            expand("results/03.genome_component/{sample}/mob-suite/{sample}_mobtyper_result.xls", sample=SAMPLES)
        output:
            "results/03.genome_component/all_mob-suite_results.xls"
        params:
            indir = "results/03.genome_component",
            outdir = "results/03.genome_component",
            sample_list = config["sample_list"],
            script = config["merge_results_script"]
        log:
            "logs/03.genome_component/mob-suite/mob_suite_merge.log"
        shell:
            """
            python {params.script} -m mob_suite -i {params.indir} -o {params.outdir} -s {params.sample_list} > {log} 2>&1
            """
elif ASSEMBLY_STRATEGY == "short":
    ## 2.1.1 Plasmidfinder: In silico detection of plasmids
    ### Run Plasmidfinder
    rule plasmidfinder_run:
        input:
            # genome assembly annotated by Prokka
            "results/03.genome_component/{sample}/prokka/{sample}.fna"
        output:
            # Plasmidfinder results
            "results/03.genome_component/{sample}/plasmidfinder/results_tab.tsv"
        params:
            outdir = "results/03.genome_component/{sample}/plasmidfinder",
            # Plasmidfinder database
            db = config["plasmidfinder_db"],
            # Plasmidfinder minimum coverage
            mincov = config["plasmidfinder_mincov"],
            # Plasmidfinder minimum hreshold for identity
            minid = config["plasmidfinder_minid"]
        log:
            "logs/03.genome_component/plasmidfinder/plasmidfinder_run/{sample}.log"
        conda:
            "../envs/plasmidfinder.yaml"
        shell:
            """
            plasmidfinder.py -i {input} -p {params.db} -o {params.outdir} \
            -l {params.mincov} -t {params.minid} -x > {log} 2>&1

            # clean up tmp files
            rm -rf {params.outdir}/tmp 2>> {log}
            """

    ### Merge Plasmidfinder results for all samples
    rule plasmidfinder_merge:
        input:
            # Plasmidfinder results
            expand("results/03.genome_component/{sample}/plasmidfinder/results_tab.tsv", sample=SAMPLES)
        output:
            "results/03.genome_component/all_plasmidfinder_results.xls"
        params:
            input_dir = "results/03.genome_component",
            output_dir = "results/03.genome_component",
            sample_list = config["sample_list"],
            script = config["merge_results_script"]
        log:
            "logs/03.genome_component/plasmidfinder/plasmidfinder_merge.log"
        shell:
            """
            python {params.script} -m plasmidfinder -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
            """

    ## 2.1.2 Platon: Identification and characterization of bacterial plasmid contigs
    ### Run Platon
    rule platon_run:
        input:
            # genome assembly annotated by Prokka
            "results/03.genome_component/{sample}/prokka/{sample}.fna"
        output:
            # Platon results
            "results/03.genome_component/{sample}/platon/{sample}.tsv"
        params:
            outdir = "results/03.genome_component/{sample}/platon",
            # Platon database
            db = config["platon_db"],
        log:
            "logs/03.genome_component/platon/platon_run/{sample}.log"
        threads:
            config["threads"]
        conda:
            "../envs/platon.yaml"
        shell:
            """
            platon --db {params.db} -o {params.outdir} -p {wildcards.sample} --threads {threads} {input} > {log} 2>&1
            """

    ### Merge Platon results for all samples
    rule platon_merge:
        input:
            # Platon results
            expand("results/03.genome_component/{sample}/platon/{sample}.tsv", sample=SAMPLES)
        output:
            "results/03.genome_component/all_platon_results.xls"
        params:
            input_dir = "results/03.genome_component",
            output_dir = "results/03.genome_component",
            sample_list = config["sample_list"],
            script = config["merge_results_script"]
        log:
            "logs/03.genome_component/platon/platon_merge.log"
        shell:
            """
            python {params.script} -m platon -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
            """
# remove Platon module in long reads ONLY mode 
elif ASSEMBLY_STRATEGY == "long" or ASSEMBLY_STRATEGY == "skip":
    ## 2.1.1 Plasmidfinder: In silico detection of plasmids
    ### Run Plasmidfinder
    rule plasmidfinder_run:
        input:
            # genome assembly annotated by Prokka
            "results/03.genome_component/{sample}/prokka/{sample}.fna"
        output:
            # Plasmidfinder results
            "results/03.genome_component/{sample}/plasmidfinder/results_tab.tsv"
        params:
            outdir = "results/03.genome_component/{sample}/plasmidfinder",
            # Plasmidfinder database
            db = config["plasmidfinder_db"],
            # Plasmidfinder minimum coverage
            mincov = config["plasmidfinder_mincov"],
            # Plasmidfinder minimum hreshold for identity
            minid = config["plasmidfinder_minid"]
        log:
            "logs/03.genome_component/plasmidfinder/plasmidfinder_run/{sample}.log"
        conda:
            "../envs/plasmidfinder.yaml"
        shell:
            """
            plasmidfinder.py -i {input} -p {params.db} -o {params.outdir} \
            -l {params.mincov} -t {params.minid} -x > {log} 2>&1

            # clean up tmp files
            rm -rf {params.outdir}/tmp 2>> {log}
            """

    ### Merge Plasmidfinder results for all samples
    rule plasmidfinder_merge:
        input:
            # Plasmidfinder results
            expand("results/03.genome_component/{sample}/plasmidfinder/results_tab.tsv", sample=SAMPLES)
        output:
            "results/03.genome_component/all_plasmidfinder_results.xls"
        params:
            input_dir = "results/03.genome_component",
            output_dir = "results/03.genome_component",
            sample_list = config["sample_list"],
            script = config["merge_results_script"]
        log:
            "logs/03.genome_component/plasmidfinder/plasmidfinder_merge.log"
        shell:
            """
            python {params.script} -m plasmidfinder -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
            """

# 2.2 Prophage
## Run phigaro
rule phigaro_run:
    input:
        "results/03.genome_component/{sample}/prokka/{sample}.fna"
    output:
        "results/03.genome_component/{sample}/phigaro/{sample}.phigaro.tsv"
    params:
        outdir = "results/03.genome_component/{sample}/phigaro"
    log:
        "logs/03.genome_component/phigaro/phigaro_run/{sample}.log"
    threads:
        config["threads"]
    conda:
        "../envs/phigaro.yaml"
    shell:
        "phigaro -f {input} -o {params.outdir} -p --not-open -e tsv html gff bed -d -t {threads} > {log} 2>&1"

## Merge Phigaro results for all samples
rule phigaro_merge:
    input:
        expand("results/03.genome_component/{sample}/phigaro/{sample}.phigaro.tsv", sample=SAMPLES)
    output:
        "results/03.genome_component/all_phigaro_results.xls"
    params:
        indir = "results/03.genome_component",
        outdir = "results/03.genome_component",
        sample_list = config["sample_list"],
        script = config["merge_results_script"]
    log:
        "logs/03.genome_component/phigaro/phigaro_merge.log"
    shell:
        """
        python {params.script} -m phigaro -i {params.indir} -o {params.outdir} -s {params.sample_list} > {log} 2>&1
        """

# 2.3 Integrative and conjugative elements (ICEs)
## ICEfinder local version - Prediction of ICE and IME in Genome Sequences of Bacteria
### Setup ICEfinder working environment
rule icefinder_setup:
    input:
        expand("results/03.genome_component/{sample}/prokka/{sample}.gbk", sample=SAMPLES)
    output:
        # signal file to indicate the end of ICEfinder environment setup
        "results/workflow_signal/03.genome_component/icefinder_setup.done"
    log:
        "logs/03.genome_component/icefinder/icefinder_setup.log"
    conda:
        "../envs/icefinder.yaml"
    shell:
        """
        # preparation for ICEfinder
        ## install a perl module
        cpanm Linux::Distribution > {log} 2>&1
        ## generate the signal file
        touch {output} 2>> {log}
        """

### Run ICEfinder
rule icefinder_run:
    input:
        # signal file to indicate the step of ICEfinder working environment setup has been done
        signal = "results/workflow_signal/03.genome_component/icefinder_setup.done",
        # standard Genbank file derived from the master .gff generated by Prokka
        gbk = "results/03.genome_component/{sample}/prokka/{sample}.gbk"
    output:
        # prediction of ICE and IME in genome sequences of bacteria by ICEfinder for each sample
        "results/03.genome_component/{sample}/icefinder/{sample}_icefinder_result.txt"
    params:
        workdir = WORKDIR,
        outdir = "results/03.genome_component/{sample}/icefinder",
        # ICEfinder resources directory
        icefinder_dir = config["icefinder_dir"],
        # script to split genbank files (from https://github.com/fmalmeida/bacannot)
        splitgenbank_script = config["splitgenbank_script"]
    log:
        "logs/03.genome_component/icefinder/icefinder_run/{sample}.log"
    conda:
        "../envs/icefinder.yaml"
    shell:
        """
        # preparation for ICEfinder
        ## copy ICEfinder resources to the working directory
        cp -r {params.icefinder_dir}/data {params.icefinder_dir}/scripts {params.icefinder_dir}/tools {params.outdir}/ 2> {log}
        cp {params.icefinder_dir}/ICEfinder_local.pl {params.outdir}/ 2>> {log}
        ln -fs $(which seqret) {params.outdir}/tools/seqret_7
        ln -fs $(which transeq) {params.outdir}/tools/transeq_7
        ## make some folders
        mkdir -p {params.outdir}/tmp 2>> {log}
        mkdir -p {params.outdir}/result 2>> {log}
        mkdir -p {params.outdir}/gbk 2>> {log}

        # split gbk file
        cp {input.gbk} {params.outdir}/gbk/annotation.gbk 2>> {log}
        cd {params.outdir}/gbk
        python {params.workdir}/{params.splitgenbank_script} annotation.gbk && rm annotation.gbk >> {params.workdir}/{log} 2>&1

        cd {params.workdir}/{params.outdir}

        # run ICEfinder
        ls ./gbk > gbk_list.txt 2>> {params.workdir}/{log}
        perl ICEfinder_local.pl gbk_list.txt >> {params.workdir}/{log} 2>&1

        # generate result summary
        echo -e "ICEfinder_Job_id\tStrain\tGenome_len\tICEfinder_output_name\tDescription\tCoordinate\tLength\toriT\tGC\tGenome_GC\t|Delta_GC|\tARG\tVF" > {params.workdir}/{output} 2>> {params.workdir}/{log} 
        shopt -s nullglob
        for f in ./result/*/*_summary.txt
        do
          cat ${{f}} >> {params.workdir}/{output} 2>> {params.workdir}/{log}
        done

        # clean up intermediate files
        rm -rf tmp data scripts tools 2>> {params.workdir}/{log}
        rm ICEfinder_local.pl 2>> {params.workdir}/{log}
        rm gbk_list.txt 2>> {params.workdir}/{log}
        rm ICEfinder.log 2>> {params.workdir}/{log}

        cd {params.workdir} 
        """

### Merge ICEfinder results for all samples
rule icefinder_merge:
    input:
        # prediction of ICE and IME in genome sequences of bacteria by ICEfinder for each sample
        expand("results/03.genome_component/{sample}/icefinder/{sample}_icefinder_result.txt", sample=SAMPLES)
    output:
        "results/03.genome_component/all_icefinder_results.xls"
    params:
        input_dir = "results/03.genome_component",
        output_dir = "results/03.genome_component",
        sample_list = config["sample_list"],
        script = config["merge_results_script"]
    log:
        "logs/03.genome_component/icefinder/icefinder_merge.log"
    shell:
        """
        python {params.script} -m icefinder -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
        """

# 2.4 Insertion sequences (IS)
## digIS: Focused detection of insertion sequences
### Note: resources/digIS-digISv1.2/src/common/grange.py has been replaced with an
###   updated version from https://github.com/fmalmeida/bacannot to fix a digIS bug
### Some digIS parameters are defined in resources/digIS-digISv1.2/definitions.py, such as NUM_THREADS = 4
### Run digIS
rule digis_run:
    input:
        # genome assembly annotated by Prokka
        fna = "results/03.genome_component/{sample}/prokka/{sample}.fna",
        # standard Genbank file derived from the master .gff generated by Prokka
        gbk = "results/03.genome_component/{sample}/prokka/{sample}.gbk"
    output:
        # digIS CSV output
        csv = "results/03.genome_component/{sample}/digis/results/{sample}.csv",
        # digIS GFF output
        gff = "results/03.genome_component/{sample}/digis/results/{sample}.gff",
        # digIS summary statistics output
        sum = "results/03.genome_component/{sample}/digis/results/{sample}.sum",
    params:
        outdir = "results/03.genome_component/{sample}/digis",
        # digISscript
        digis_script = config["digis_script"]
    log:
        "logs/03.genome_component/digis/digis_run/{sample}.log"
    threads:
        4
    conda:
        "../envs/digis.yaml"
    shell:
        """
        python {params.digis_script} -i {input.fna} -g {input.gbk} -o {params.outdir} > {log} 2>&1
    
        # remove intermediate output files
        rm -rf {params.outdir}/hmmer 2>> {log}
        rm -rf {params.outdir}/logs 2>> {log}
        rm -rf {params.outdir}/pep 2>> {log}
        """

### Merge digIS results for all samples
rule digis_merge:
    input:
        # digIS results
        expand("results/03.genome_component/{sample}/digis/results/{sample}.csv", sample=SAMPLES)
    output:
        "results/03.genome_component/all_digis_results.xls"
    params:
        input_dir = "results/03.genome_component",
        output_dir = "results/03.genome_component",
        sample_list = config["sample_list"],
        script = config["merge_results_script"]
    log:
        "logs/03.genome_component/digis/digis_merge.log"
    shell:
        """
        python {params.script} -m digis -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
        """

# 2.5 Genomic islands (GIs)
## IslandPath-DIMOB: Prediction and visualization of genomic islands
### Note: IslandPath-DIMOB ONLY works with single contig inputs
### Thus use a custom script (from https://github.com/fmalmeida/bacannot) to split GBK file generated by Prokka
### Run IslandPath-DIMOB
rule islandpath_run:
    input:
        # standard Genbank file derived from the master .gff generated by Prokka
        "results/03.genome_component/{sample}/prokka/{sample}.gbk"
    output:
        # bed file of genomic islands predicted by IslandPath-DIMOB
        "results/03.genome_component/{sample}/islandpath/{sample}_GIs.bed",
    params:
        workdir = WORKDIR,
        outdir = "results/03.genome_component/{sample}/islandpath",
        # script to split genbank files (from https://github.com/fmalmeida/bacannot)
        splitgenbank_script = config["splitgenbank_script"]
    log:
        "logs/03.genome_component/islandpath/islandpath_run/{sample}.log"
    conda:
        "../envs/islandpath.yaml"
    shell:
        """  
        # split genbank files
        mkdir -p {params.outdir}/gbk_input 
        cp {input} {params.outdir}/gbk_input/annotation.gbk
        cd {params.outdir}/gbk_input
        python {params.workdir}/{params.splitgenbank_script} annotation.gbk && rm annotation.gbk 2> {params.workdir}/{log}
        cd {params.workdir}

        # run islandpath in each
        touch {output}
        # have to run islandpath in this directory since it will generate a log file
        cd {params.outdir}
        for file in $(ls gbk_input/*.gbk)
        do
            filename=$(basename $file)
            touch ${{filename%%.gbk}}_GIs.txt
            grep -q "CDS" $file && islandpath $file ${{filename%%.gbk}}_GIs.txt >> {params.workdir}/{log} 2>&1
            name="${{filename%%.gbk}}"
            awk -v contig=$name 'BEGIN {{ FS = "\\t"; OFS="\\t" }} {{ print contig,$2,$3 }}' ${{filename%%.gbk}}_GIs.txt >> {wildcards.sample}_GIs.bed 2>> {params.workdir}/{log}
        done
        sort -nk 1,2 -o {wildcards.sample}_GIs.bed {wildcards.sample}_GIs.bed 2>> {params.workdir}/{log}
        cd {params.workdir}
        """

### Merge IslandPath-DIMOB results for all samples
rule islandpath_merge:
    input:
        # IslandPath-DIMOB results
        expand("results/03.genome_component/{sample}/islandpath/{sample}_GIs.bed", sample=SAMPLES)
    output:
        "results/03.genome_component/all_islandpath_results.xls"
    params:
        input_dir = "results/03.genome_component",
        output_dir = "results/03.genome_component",
        sample_list = config["sample_list"],
        script = config["merge_results_script"]
    log:
        "logs/03.genome_component/islandpath/islandpath_merge.log"
    shell:
        """
        python {params.script} -m islandpath -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
        """

# 2.6 CRISPRs
## CRISPRCasTyper: Detect CRISPR-Cas genes and arrays,
##   and predict the subtype based on both Cas genes and CRISPR repeat sequence
rule cctyper:
    input:
        "results/02.assembly/{sample}/assembly/{sample}.fasta"
    output:
        "results/03.genome_component/{sample}/cctyper/cctyper.done"
    params:
        outdir = "results/03.genome_component/{sample}/cctyper/results"
    log:
        "logs/03.genome_component/cctyper/cctyper_run/{sample}.log"
    conda:
        "../envs/cctyper.yaml"
    threads:
        16
    shell:
        """
        cctyper {input} {params.outdir} --circular -t {threads} > {log} 2>&1
        touch {output} 2>> {log}
        """

### Merge CRISPRCasTyper results for all samples
rule cctyper_merge:
    input:
        expand("results/03.genome_component/{sample}/cctyper/cctyper.done", sample=SAMPLES)
    output:
        "results/03.genome_component/all_cctyper_results.xls"
    params:
        input_dir = "results/03.genome_component",
        output_dir = "results/03.genome_component",
        sample_list = config["sample_list"],
        script = config["merge_results_script"]
    log:
        "logs/03.genome_component/cctyper/cctyper_merge.log"
    shell:
        """
        python {params.script} -m cctyper -i {params.input_dir} -o {params.output_dir} -s {params.sample_list} > {log} 2>&1
        """