#######################################################
# Step 01. Raw sequencing reads QC and pre-processing #
#######################################################

##### 1. Data cleaning #####
# 1.1 Fastp: Illumina fastq reads cleaning
## Run fastp
rule fastp:
    input:
        fq1 = config["illumina_data_dir"] + "%s/{sample}_%s1.%s" % (DATADIRINDEX2, FASTQREAD2, FASTQSUFFIX2),
        fq2 = config["illumina_data_dir"] + "%s/{sample}_%s2.%s" % (DATADIRINDEX2, FASTQREAD2, FASTQSUFFIX2)
    output:
        cfq1 = "results/01.data_clean/short_read/cleandata/{sample}_1.clean.fastq.gz",
        cfq2 = "results/01.data_clean/short_read/cleandata/{sample}_2.clean.fastq.gz",
        # Hopefully the u file is very small relative to the 1 and 2 files.If not,something might have gone wrong!
        cfqu = "results/01.data_clean/short_read/{sample}/fastp/{sample}_u.clean.fastq.gz",
        json = "results/01.data_clean/short_read/{sample}/fastp/{sample}_fastp.json",
        html = "results/01.data_clean/short_read/{sample}/fastp/{sample}_fastp.html"
    threads:
        config["threads"]
    log:
        "logs/01.data_clean/short_read/fastp/{sample}.log"
    conda:
        "../envs/data_clean.yaml"
    shell:
        """
        fastp --detect_adapter_for_pe -w {threads} -i {input.fq1} -I {input.fq2} \
        --unpaired1 {output.cfqu} --unpaired2 {output.cfqu} \
        -o {output.cfq1} -O {output.cfq2} -j {output.json} -h {output.html} > {log} 2>&1
        """

## Extract fastp QC statistics
rule fastp_stats:
    input:
        expand("results/01.data_clean/short_read/{sample}/fastp/{sample}_fastp.json", sample=SAMPLES)
    output:
        "results/01.data_clean/short_read/all_fastp_stats.xls"
    params:
        input_dir = "results/01.data_clean/short_read/",
        sample_list = config["sample_list"],
        script = config["script_extract_fastp_stats"]
    log:
        "logs/01.data_clean/short_read/fastp_stats.log"
    shell:
        """
        python {params.script} {params.input_dir} {params.sample_list} {output} > {log} 2>&1
        """

# 1.2 Filtlong: long fastq reads cleaning
rule filtlong:
    input:
        config["long_data_dir"] + "%s/{sample}.%s" % (DATADIRINDEX3, FASTQSUFFIX3)
    output:
        "results/01.data_clean/long_read/cleandata/{sample}.clean.fastq.gz",
    threads:
        config["threads"]
    log:
        "logs/01.data_clean/long_read/filtlong/{sample}.log"
    conda:
        "../envs/data_clean.yaml"
    shell:
        """
        (filtlong --min_length 1000 --keep_percent 95 {input} | gzip) > {output} 2> {log}
        """

##### 2. Species identification (Kraken2 + Bracken) #####
# 2.1 Generate a Kraken2 report
rule kraken2:
    input:
        fq1 = "results/01.data_clean/short_read/cleandata/{sample}_1.clean.fastq.gz",
        fq2 = "results/01.data_clean/short_read/cleandata/{sample}_2.clean.fastq.gz"
    output:
        report = "results/01.data_clean/short_read/{sample}/kraken2/{sample}.kraken2.report.xls"
    log:
        "logs/01.data_clean/short_read/kraken2/{sample}.log"
    conda:
        "../envs/kraken2.yaml"
    threads:
        config["threads"]
    params:
        kraken2_db = config["kraken2_db"]
    shell:
        """
        kraken2 --db {params.kraken2_db} --threads {threads} --report {output.report} \
        --paired --gzip-compressed {input.fq1} {input.fq2} --output - 2> {log}
        """

# 2.2 Estimate species-level abundance based on Kraken2 output by Bracken
rule bracken:
    input:
        "results/01.data_clean/short_read/{sample}/kraken2/{sample}.kraken2.report.xls"
    output:
        reads = "results/01.data_clean/short_read/{sample}/bracken/{sample}.bracken",
        report = "results/01.data_clean/short_read/{sample}/bracken/{sample}.bracken.report.xls"
    log:
        "logs/01.data_clean/short_read/bracken/{sample}.log"
    conda:
        "../envs/kraken2.yaml"
    threads:
        config["threads"]
    params:
        kraken2_db = config["kraken2_db"],
        # level to estimate abundance
        level = config["bracken_level"],
        # read length to get all classifications
        read_len = config["bracken_read_len"]
    shell:
        """
        bracken -d {params.kraken2_db} -i {input} -o {output.reads} -w {output.report} \
        -r {params.read_len} -l {params.level} > {log} 2>&1
        """

# 2.3 Extract species related info from Bracken report
rule bracken_filter:
    input:
        "results/01.data_clean/short_read/{sample}/bracken/{sample}.bracken.report.xls"
    output:
        "results/01.data_clean/short_read/{sample}/{sample}.kraken2.result.xls"
    log:
        "logs/01.data_clean/short_read/bracken_filter/{sample}.log"
    params:
        script =  config["script_extract_bracken"],
        genus = config["genus"],
        species = config["species"]
    shell:
        """
        echo -e "##The percentage of reads mapping to {params.genus} {params.species} and the main species detected by kraken2" > {output} 2> {log}
        echo -e "#sample\tpercentage_target\tmain_species_detected\tpercentage_main_species" >> {output} 2>> {log}
        python {params.script} {input} {output} {wildcards.sample} {params.genus} {params.species} >> {log} 2>&1
        """

# 2.4 Merge all sample's Kraken2 results
rule kraken2_qc:
    input:
        expand("results/01.data_clean/short_read/{sample}/{sample}.kraken2.result.xls", sample=SAMPLES)
    output:
        qc_result = "results/01.data_clean/short_read/all_kraken2_result.xls"
    params:
        qc_threshold = config["kraken2_threshold"],
        genus = config["genus"],
        species = config["species"]
    log:
        "logs/01.data_clean/short_read/kraken2_qc.log"
    shell:
        """
        # Generate kraken qc result output
        echo -e "##The percentage of reads mapping to {params.genus} {params.species} and the main species detected by kraken" > {output.qc_result} 2> {log}
        echo -e "#sample\tpercentage_target\tmain_species_detected\tpercentage_main_species" >> {output.qc_result} 2>> {log}
        cat {input} | grep -v "^#" >> {output.qc_result} 2>> {log}

        # extract qualified sample list
        TMP=`head -n 2 {output.qc_result} && grep -v "^#" {output.qc_result} | sort -k2n` 2>> {log}
        echo "$TMP" > {output.qc_result} 2>> {log}
        """

##### 3. FASTQ reads quality control #####
# 3.1 Illumina short reads QC (FastQC & MultiQC)
## Generate a quality control report of Illumina short fastq reads for each sample by FastQC
rule fastqc:
    input:
        fq1 = "results/01.data_clean/short_read/cleandata/{sample}_1.clean.fastq.gz",
        fq2 = "results/01.data_clean/short_read/cleandata/{sample}_2.clean.fastq.gz",
    output:
        fq1_html = "results/01.data_clean/short_read/fastqc/{sample}_1_fastqc.html",
        fq1_zip = "results/01.data_clean/short_read/fastqc/{sample}_1_fastqc.zip",
        fq2_html = "results/01.data_clean/short_read/fastqc/{sample}_2_fastqc.html",
        fq2_zip = "results/01.data_clean/short_read/fastqc/{sample}_2_fastqc.zip",
    log:
        "logs/01.data_clean/short_read/fastqc/{sample}.log"
    params:
        outdir = "results/01.data_clean/short_read/fastqc",
        fq1_prefix = "{sample}_1.clean",
        fq2_prefix = "{sample}_2.clean"
    threads:
        config["threads"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        # run fastqc
        fastqc -t {threads} -o {params.outdir} {input.fq1} {input.fq2} > {log} 2>&1
        # modify fastqc result file names
        mv {params.outdir}/{params.fq1_prefix}_fastqc.html {params.outdir}/{wildcards.sample}_1_fastqc.html 2>> {log}
        mv {params.outdir}/{params.fq1_prefix}_fastqc.zip {params.outdir}/{wildcards.sample}_1_fastqc.zip 2>> {log}
        mv {params.outdir}/{params.fq2_prefix}_fastqc.html {params.outdir}/{wildcards.sample}_2_fastqc.html 2>> {log}
        mv {params.outdir}/{params.fq2_prefix}_fastqc.zip {params.outdir}/{wildcards.sample}_2_fastqc.zip 2>> {log}
        """

## MultiQC: merge all fastqc reports
rule multiqc_fastqc:
    input:
        expand("results/01.data_clean/short_read/fastqc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=[1,2])
    output:
        "results/01.data_clean/short_read/multiqc/multiqc_report.html"
    log:
        "logs/01.data_clean/short_read/multiqc.log"
    params:
        indir = "results/01.data_clean/short_read/fastqc",
        outdir = "results/01.data_clean/short_read/multiqc/"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.indir} -m fastqc -f -o {params.outdir} > {log} 2>&1
        """

# 3.2 Nanopore long reads QC (NanoPlot & NanoStat)
## Generate a quality control report of Nanopore long fastq raw reads for each sample by NanoStat (supported by MultiQC)
rule nanostat_raw:
    input:
        config["long_data_dir"] + "%s/{sample}.%s" % (DATADIRINDEX3, FASTQSUFFIX3)
    output:
        "results/00.rawdata/long_read/nanostat_raw/{sample}_raw_NanoStats.txt"
    log:
        "logs/00.rawdata/long_read/nanostat_raw/{sample}.log"
    threads:
        config["threads"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        NanoStat --fastq {input} -n {output} -t {threads} > {log} 2>&1
        """

## Extract info from NanoStat results and merge
rule merge_nanostat_raw:
    input:
        expand("results/00.rawdata/long_read/nanostat_raw/{sample}_raw_NanoStats.txt", sample=SAMPLES)
    output:
        "results/00.rawdata/long_read/all_raw_NanoStats.xls"
    log:
        "logs/00.rawdata/long_read/nanostat_raw_merge.log"
    params:
        sample_list = config["sample_list"],
        indir = "results/00.rawdata/long_read/nanostat_raw",
        outdir = "results/00.rawdata/long_read",
        # Python script to extract info from NanoStat results and merge
        script = config["merge_nanostat_script"]
    shell:
        """
        python {params.script} {params.sample_list} {params.indir} {params.outdir} raw > {log} 2>&1
        """

## MultiQC: merge all NanoStat reports using
rule multiqc_nanostat_raw:
    input:
        expand("results/00.rawdata/long_read/nanostat_raw/{sample}_raw_NanoStats.txt", sample=SAMPLES)
    output:
        "results/00.rawdata/long_read/multiqc_raw/multiqc_report.html"
    log:
        "logs/00.rawdata/long_read/multiqc_raw.log"
    params:
        indir = "results/00.rawdata/long_read/nanostat_raw",
        outdir = "results/00.rawdata/long_read/multiqc_raw/"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.indir} -m nanostat -f -o {params.outdir} > {log} 2>&1
        """

## Generate a quality control report of Nanopore long fastq clean reads for each sample by NanoPlot
rule nanoplot_clean:
    input:
        "results/01.data_clean/long_read/cleandata/{sample}.clean.fastq.gz"
    output:
        "results/01.data_clean/long_read/{sample}/nanoplot_clean/{sample}_clean_NanoPlot-report.html"
    log:
        "logs/01.data_clean/long_read/nanoplot_clean/{sample}.log"
    params:
        outdir = "results/01.data_clean/long_read/{sample}/nanoplot_clean/",
        prefix = "{sample}_clean_"
    threads:
        config["threads"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        NanoPlot -t {threads} -o {params.outdir} -p {params.prefix} --tsv_stats --title {wildcards.sample} \
        -f png --fastq {input} --loglength --N50 --verbose > {log} 2>&1
        """

## Generate a quality control report of Nanopore long fastq clean reads for each sample by NanoStat (supported by MultiQC)
rule nanostat_clean:
    input:
        "results/01.data_clean/long_read/cleandata/{sample}.clean.fastq.gz"
    output:
        "results/01.data_clean/long_read/nanostat_clean/{sample}_clean_NanoStats.txt"
    log:
        "logs/01.data_clean/long_read/nanostat_clean/{sample}.log"
    threads:
        config["threads"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        NanoStat --fastq {input} -n {output} -t {threads} > {log} 2>&1
        """

## Extract info from NanoStat results and merge
rule merge_nanostat_clean:
    input:
        expand("results/01.data_clean/long_read/nanostat_clean/{sample}_clean_NanoStats.txt", sample=SAMPLES)
    output:
        "results/01.data_clean/long_read/all_clean_NanoStats.xls"
    log:
        "logs/01.data_clean/long_read/nanostat_clean_merge.log"
    params:
        sample_list = config["sample_list"],
        indir = "results/01.data_clean/long_read/nanostat_clean",
        outdir = "results/01.data_clean/long_read",
        # Python script to extract info from NanoStat results and merge
        script = config["merge_nanostat_script"]
    shell:
        """
        python {params.script} {params.sample_list} {params.indir} {params.outdir} clean > {log} 2>&1
        """

## MultiQC: merge all NanoStat reports using
rule multiqc_nanostat_clean:
    input:
        expand("results/01.data_clean/long_read/nanostat_clean/{sample}_clean_NanoStats.txt", sample=SAMPLES)
    output:
        "results/01.data_clean/long_read/multiqc_clean/multiqc_report.html"
    log:
        "logs/01.data_clean/long_read/multiqc_clean.log"
    params:
        indir = "results/01.data_clean/long_read/nanostat_clean",
        outdir = "results/01.data_clean/long_read/multiqc_clean/"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.indir} -m nanostat -f -o {params.outdir} > {log} 2>&1
        """