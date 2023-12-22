#################################
# Step 06. Genome Visualization #
#################################

##### 1. Chromosome Visualization #####
# MGCplotter: Microbial Genome Circular plotter
rule mgcplotter_chr:
    input:
        # prokka annotated Genbank file
        "results/03.genome_component/{sample}/prokka/{sample}.gbk"
    output:
        "results/06.genome_visualization/{sample}/mgcplotter/mgcplotter_chr.done"
    params:
        workdir = WORKDIR,
        outdir = "results/06.genome_visualization/{sample}/mgcplotter",
        # script to split genbank files
        splitgenbank_script = "workflow/scripts/splitgenbank.py",
    log:
        "logs/06.genome_visualization/mgcplotter_chr/{sample}.log"
    threads:
        64
    conda:
        "mgcplotter"
    shell:
        """
        # split genbank files
        mkdir -p {params.outdir}/gbk_input 
        cp {input} {params.outdir}/gbk_input/annotation.gbk
        cd {params.outdir}/gbk_input
        python {params.workdir}/{params.splitgenbank_script} annotation.gbk && rm annotation.gbk 2> {params.workdir}/{log}
        cd {params.workdir}

        # Circos plot for chromosome
        cd {params.outdir}
        shopt -s nullglob
        for f in gbk_input/chr*.gbk
        do
          chr_num=`basename $f .gbk`
          chr_len=`head -n1 $f | awk '{{print $3}}'`
           
          # Run MGCplotter with no COG assignment
          mkdir -p without_cog/{wildcards.sample}_${{chr_num}} 2>> {params.workdir}/{log}
          MGCplotter -r $f -o without_cog/{wildcards.sample}_${{chr_num}} \
          --gc_content_p_color orange --gc_content_n_color blue \
          --gc_skew_p_color pink --gc_skew_n_color green >> {params.workdir}/{log} 2>&1
          # Add some texts
          convert without_cog/{wildcards.sample}_${{chr_num}}/circos.png -font Ubuntu -gravity center \
          -pointsize 50 -fill black -annotate 0x0+0-20 "{wildcards.sample}_${{chr_num}}" \
          -annotate 0x0+0+50 "${{chr_len}} bp" \
          without_cog/{wildcards.sample}_${{chr_num}}.png >> {params.workdir}/{log} 2>&1

          # Run MGCplotter with COG assignment
          mkdir -p with_cog/{wildcards.sample}_${{chr_num}} 2>> {params.workdir}/{log}
          MGCplotter -r $f -o with_cog/{wildcards.sample}_${{chr_num}} \
          --assign_cog_color --gc_content_p_color orange \
          --gc_content_n_color blue --gc_skew_p_color pink --gc_skew_n_color green >> {params.workdir}/{log} 2>&1
          # Add some texts
          convert with_cog/{wildcards.sample}_${{chr_num}}/circos.png -font Ubuntu -gravity center \
          -pointsize 50 -fill black -annotate 0x0+0-20 "{wildcards.sample}_${{chr_num}}" \
          -annotate 0x0+0+50 "${{chr_len}} bp" \
          with_cog/{wildcards.sample}_${{chr_num}}.png >> {params.workdir}/{log} 2>&1
        done

        # Create a flag file when everything done
        touch {params.workdir}/{output} 2>> {params.workdir}/{log}
        cd {params.workdir}
        """

##### 2. Plasmid Visualization #####
# MGCplotter: Microbial Genome Circular plotter
rule mgcplotter_plas:
    input:
        # run after rule mgcplotter_chr
        "results/06.genome_visualization/{sample}/mgcplotter/mgcplotter_chr.done"
    output:
        "results/06.genome_visualization/{sample}/mgcplotter/mgcplotter_plas.done"
    params:
        workdir = WORKDIR,
        outdir = "results/06.genome_visualization/{sample}/mgcplotter",
        # not draw a plot for too short plasmid since MGCplotter might throw an error
        plas_length_threshold = config["plas_length_threshold"],
    log:
        "logs/06.genome_visualization/mgcplotter_plas/{sample}.log"
    threads:
        64
    conda:
        "mgcplotter"
    shell:
        """
        # Circos plot for plasmid
        cd {params.outdir}
        shopt -s nullglob
        for f in gbk_input/plas*.gbk
        do
          plas_num=`basename $f .gbk`
          plas_len=`head -n1 $f | awk '{{print $3}}'`
          
          if [ "${{plas_len}}" -gt "{params.plas_length_threshold}" ]; then
              # Run MGCplotter with no COG assignment
              mkdir -p without_cog/{wildcards.sample}_${{plas_num}} 2>> {params.workdir}/{log}
              MGCplotter -r $f -o without_cog/{wildcards.sample}_${{plas_num}} \
              --gc_content_p_color orange --gc_content_n_color blue \
              --gc_skew_p_color pink --gc_skew_n_color green >> {params.workdir}/{log} 2>&1
              # Add some texts
              convert without_cog/{wildcards.sample}_${{plas_num}}/circos.png -font Ubuntu -gravity center \
              -pointsize 50 -fill black -annotate 0x0+0-20 "{wildcards.sample}_${{plas_num}}" \
              -annotate 0x0+0+50 "${{plas_len}} bp" \
              without_cog/{wildcards.sample}_${{plas_num}}.png >> {params.workdir}/{log} 2>&1
    
              # Run MGCplotter with COG assignment
              mkdir -p with_cog/{wildcards.sample}_${{plas_num}} 2>> {params.workdir}/{log}
              MGCplotter -r $f -o with_cog/{wildcards.sample}_${{plas_num}} \
              --assign_cog_color --gc_content_p_color orange \
              --gc_content_n_color blue --gc_skew_p_color pink --gc_skew_n_color green >> {params.workdir}/{log} 2>&1
              # Add some texts
              convert with_cog/{wildcards.sample}_${{plas_num}}/circos.png -font Ubuntu -gravity center \
              -pointsize 50 -fill black -annotate 0x0+0-20 "{wildcards.sample}_${{plas_num}}" \
              -annotate 0x0+0+50 "${{plas_len}} bp" \
              with_cog/{wildcards.sample}_${{plas_num}}.png >> {params.workdir}/{log} 2>&1
          fi
        done

        # Create a flag file when everything done
        touch {params.workdir}/{output} 2>> {params.workdir}/{log}
        cd {params.workdir}
        """
