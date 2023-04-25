#####################
# Step 02. Assembly #
#####################

##### 1. Genome assembly #####
# Illumina short reads + Nanopore or PacBio Long reads
if ASSEMBLY_STRATEGY == "hybrid":
    # 1.1 Unicycler hybrid assembly
    ## Hybrid assembly (using both Illumina read and long reads) is where Unicycler really shines.
    ##   Like with the Illumina-only pipeline described above, Unicycler will produce an Illumina assembly graph.
    ##   It then uses long reads to build bridges, which often allows it to resolve all repeats in the genome, resulting in a complete genome assembly.
    ## https://github.com/rrwick/Unicycler#method-hybrid-assembly
    rule unicycler:
        input:
            short_fq1 = "results/01.data_clean/short_read/cleandata/{sample}_1.clean.fastq.gz",
            short_fq2 = "results/01.data_clean/short_read/cleandata/{sample}_2.clean.fastq.gz",
            long_fq = "results/01.data_clean/long_read/cleandata/{sample}.clean.fastq.gz"
        output:
            assembly = "results/02.assembly/{sample}/unicycler/assembly.fasta"
        params:
            outdir = "results/02.assembly/{sample}/unicycler",
            # Path to the SPAdes 3.15.4 executable since in unicycler 0.5.0 conda env
            #   SPAdes version is 3.15.3 which is not compatible with Python 3.10.0
            spades = config["spades_path"]
        log:
            "logs/02.assembly/unicycler/{sample}.log"
        threads:
            config["threads"]
        conda:
            "../envs/assembly.yaml"
        shell:
            """
            unicycler -1 {input.short_fq1} -2 {input.short_fq2} -l {input.long_fq} \
            -o {params.outdir} --spades_path {params.spades} -t {threads} \
            --min_fasta_length 1000 > {log} 2>&1
            """
elif ASSEMBLY_STRATEGY == "long":
    # 1.1 Unicycler long-read-only assembly
    rule unicycler:
        input:
            long_fq = "results/01.data_clean/long_read/cleandata/{sample}.clean.fastq.gz"
        output:
            assembly = "results/02.assembly/{sample}/unicycler/assembly.fasta"
        params:
            outdir = "results/02.assembly/{sample}/unicycler",
        log:
            "logs/02.assembly/unicycler/{sample}.log"
        threads:
            config["threads"]
        conda:
            "../envs/assembly.yaml"
        shell:
            """
            unicycler -l {input.long_fq} -o {params.outdir} -t {threads} > {log} 2>&1
            """
elif ASSEMBLY_STRATEGY == "short":
    # 1.1 Unicycler illumina-read-only assembly
    rule unicycler:
        input:
            short_fq1 = "results/01.data_clean/short_read/cleandata/{sample}_1.clean.fastq.gz",
            short_fq2 = "results/01.data_clean/short_read/cleandata/{sample}_2.clean.fastq.gz",
        output:
            assembly = "results/02.assembly/{sample}/unicycler/assembly.fasta"
        params:
            outdir = "results/02.assembly/{sample}/unicycler",
            # Path to the SPAdes 3.15.4 executable since in unicycler 0.5.0 conda env
            #   SPAdes version is 3.15.3 which is not compatible with Python 3.10.0
            spades = config["spades_path"]
        log:
            "logs/02.assembly/unicycler/{sample}.log"
        threads:
            config["threads"]
        conda:
            "../envs/assembly.yaml"
        shell:
            """
            unicycler -1 {input.short_fq1} -2 {input.short_fq2} \
            -o {params.outdir} --spades_path {params.spades} -t {threads} > {log} 2>&1
            """

##### 2. Separate chromosome and plasmid #####
if ASSEMBLY_STRATEGY == "hybrid":
    # 2.1 Identify plasmid contig in genome assembly using viralVerify
    ## https://github.com/ablab/viralVerify
    rule viralverify:
        input:
            "results/02.assembly/{sample}/unicycler/assembly.fasta"
        output:
            "results/02.assembly/{sample}/viralverify/assembly_result_table.csv"
        params:
            outdir = "results/02.assembly/{sample}/viralverify/",
            # Path to viralVerify HMM database
            hmm = config["viralverify_hmm_path"]
        log:
            "logs/02.assembly/viralverify/{sample}.log"
        threads:
            config["threads"]
        conda:
            "../envs/assembly.yaml"
        shell:
            """
            viralverify -f {input} -o {params.outdir} --hmm {params.hmm} -t {threads} > {log} 2>&1
            """

    # # 2.1 Identify plasmid contig in genome assembly using MOB-recon
    # rule mob_recon:
    #     input:
    #         "results/02.assembly/{sample}/unicycler/assembly.fasta"
    #     output:
    #         "results/02.assembly/{sample}/mob_recon/contig_report.txt"
    #     params:
    #         outdir = "results/02.assembly/{sample}/mob_recon"
    #     log:
    #         "logs/02.assembly/mob_recon/{sample}.log"
    #     threads:
    #         config["threads"]
    #     conda:
    #         "../envs/mob-suite.yaml"
    #     shell:
    #         """
    #         mob_recon -i {input} -o {params.outdir} -n {threads} --force > {log} 2>&1
    #         """

    # 2.2 Separate assembly fasta into several fasta of chromosome and plasmid based on MOB-recon results
    rule separate_assembly:
        input:
            assembly = "results/02.assembly/{sample}/unicycler/assembly.fasta",
            viralverify = "results/02.assembly/{sample}/viralverify/assembly_result_table.csv"
            #mob_recon = "results/02.assembly/{sample}/mob_recon/contig_report.txt"
        output:
            assembly_all = "results/02.assembly/{sample}/assembly/{sample}.fasta",
            flag = "results/02.assembly/{sample}/assembly/separate_assembly.done"
        params:
            # Path to seqkit executable
            seqkit = config["seqkit"],
            outdir = "results/02.assembly/{sample}/assembly",
            # Python script to rename separated fasta files based on viralVerify results
            script = config["script_rename_fasta"]
        log:
            "logs/02.assembly/separate_assembly/{sample}.log"
        threads:
            config["threads"]
        shell:
            """
            cp -rva {input.assembly} {output.assembly_all} > {log} 2>&1
            
            # Split genome assembly fasta into several fasta files of each contig
            {params.seqkit} split -i -O {params.outdir} {output.assembly_all} >> {log} 2>&1
            
            # Rename separated fasta files based on viralVerify results
            python {params.script} {params.outdir} {wildcards.sample} {input.viralverify} >> {log} 2>&1
            
            # Remove separated fasta files which are identified as neither chromosome nor plasmid
            rm -rf {params.outdir}/{wildcards.sample}.part*.fasta 2>> {log}
            
            # Create a flag file when everything done
            touch {output.flag} 2>> {log}
            """

    # 2.3 Summarize genome assembly statistics to be presented in summary report
    rule assembly_stats:
        input:
            lambda wildcards: get_qualified_results("results/02.assembly/{sample}/assembly/{sample}.fasta", wildcards)
        output:
            "results/02.assembly/all_assembly_stats.xls"
        params:
            sample_list = SAMPLE_LIST_UPDATE,
            indir = "results/02.assembly",
            outdir = "results/02.assembly",
            # Python script to summarize genome assembly statistics
            script = config["script_assembly_stats"]
        log:
            "logs/02.assembly/assembly_stats.log"
        conda:
            # Conda env with biopython
            "../envs/icefinder.yaml"
        shell:
            """
            python {params.script} {params.sample_list} {params.indir} {params.outdir} hybrid > {log} 2>&1
            """
elif ASSEMBLY_STRATEGY == "short" or ASSEMBLY_STRATEGY == "long":
    # 2.2 No need to split fasta and rename contigs if short-ONLY or long-ONLY mode was selected
    rule separate_assembly:
        input:
            assembly = "results/02.assembly/{sample}/unicycler/assembly.fasta",
        output:
            assembly_all = "results/02.assembly/{sample}/assembly/{sample}.fasta",
            flag = "results/02.assembly/{sample}/assembly/separate_assembly.done"
        log:
            "logs/02.assembly/separate_assembly/{sample}.log"
        shell:
            """
            cp -rva {input.assembly} {output.assembly_all} > {log} 2>&1

            # Create a flag file when everything done
            touch {output.flag} 2>> {log}
            """
elif ASSEMBLY_STRATEGY == "skip":
    # 2.2 If assembly fasta files are input directly, just use them for the rest of analysis
    rule separate_assembly:
        input:
            assembly = config["assembly_input_dir"] + "/{sample}.fasta",
        output:
            assembly_all = "results/02.assembly/{sample}/assembly/{sample}.fasta",
            flag = "results/02.assembly/{sample}/assembly/separate_assembly.done"
        log:
            "logs/02.assembly/separate_assembly/{sample}.log"
        shell:
            """
            cp -rva {input.assembly} {output.assembly_all} > {log} 2>&1

            # Create a flag file when everything done
            touch {output.flag} 2>> {log}
            """

##### 3. Genome assembly QC #####
# 3.1 Quast: generate genome assembly quality control reports for each sample
if config["ref_genome"] and config["ref_genome_gff"]:
    rule quast:
        input:
            flag = "results/02.assembly/{sample}/assembly/separate_assembly.done",
            assembly = "results/02.assembly/{sample}/assembly/{sample}.fasta"
        output:
            "results/02.assembly/quast/{sample}/report.html"
        log:
            "logs/02.assembly/quast/{sample}.log"
        params:
            outdir = "results/02.assembly/quast/{sample}",
            ref_fna = config["ref_genome"],
            ref_gff = config["ref_genome_gff"]
        threads:
            config["threads"]
        conda:
            "../envs/qc.yaml"
        shell:
            """
            quast -o {params.outdir} -r {params.ref_fna} -g {params.ref_gff} -t {threads} {input.assembly} > {log} 2>&1
            """
else:
    rule quast:
        input:
            flag="results/02.assembly/{sample}/assembly/separate_assembly.done",
            assembly="results/02.assembly/{sample}/assembly/{sample}.fasta"
        output:
            "results/02.assembly/quast/{sample}/report.html"
        log:
            "logs/02.assembly/quast/{sample}.log"
        params:
            outdir="results/02.assembly/quast/{sample}"
        threads:
            config["threads"]
        conda:
            "../envs/qc.yaml"
        shell:
            """
            quast -o {params.outdir} -t {threads} {input.assembly} > {log} 2>&1
            """

if ASSEMBLY_STRATEGY == "short" or ASSEMBLY_STRATEGY == "long" or ASSEMBLY_STRATEGY == "skip":
    # Summarize genome assembly statistics to be presented in summary report
        rule assembly_stats:
            input:
                lambda wildcards: get_qualified_results("results/02.assembly/quast/{sample}/report.html", wildcards)
            output:
                "results/02.assembly/all_assembly_stats.xls"
            params:
                sample_list = SAMPLE_LIST_UPDATE,
                indir = "results/02.assembly",
                outdir = "results/02.assembly",
                # Python script to summarize genome assembly statistics
                script = config["script_assembly_stats"]
            log:
                "logs/02.assembly/assembly_stats.log"
            conda:
                # Conda env with biopython
                "../envs/icefinder.yaml"
            shell:
                """
                python {params.script} {params.sample_list} {params.indir} {params.outdir} short > {log} 2>&1
                """

# 3.2 MultiQC: merge all Quast reports
rule multiqc_quast:
    input:
        lambda wildcards: get_qualified_results("results/02.assembly/quast/{sample}/report.html", wildcards)
    output:
        "results/02.assembly/multiqc/multiqc_report.html"
    params:
        indir = "results/02.assembly/quast",
        outdir = "results/02.assembly/multiqc/"
    log:
        "logs/02.assembly/multiqc.log"
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc {params.indir} -m quast -f -o {params.outdir} > {log} 2>&1"



