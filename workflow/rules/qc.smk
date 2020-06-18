#_____ RUN READ QC  __________________________________________________________#

rule run_pycoqc:
    input:
        "run_data/{run}/copy_finished"
    output:
        html = "qc/per_run/{run}.pycoQC.html",
        json = "qc/per_run/{run}.pycoQC.json"
    log:
        "logs/{run}_pycoqc.log"
    conda:
        "../env/pycoqc.yml"
    threads:
        1
    shell:
        """
        file=$(find run_data/{wildcards.run} -name 'sequencing_summary*')
        pycoQC \
            --summary_file $file\
            --html_outfile {output.html} \
            --json_outfile {output.json} \
            > {log} 2>&1
        """

rule run_multiqc:
    input:
        expand("qc/per_run/{run}.pycoQC.html", run = ID_runs)
    output:
        "qc/per_run/run_multiqc_report.html"
    log:
        "logs/run_multiqc.log"
    threads:
       1
    params:
        multiqc = config['apps']['multiqc']
    shell:
        """
        {params.multiqc} \
            --force \
            --config config/multiqc.yml \
            --outdir qc/per_run \
            --filename run_multiqc_report \
            qc/per_run/ \
            > {log} 2>> {log}
        """

#_____ SAMPLE READ QC  _________________________________________________________#

rule sample_pycoqc:
    input:
        unpack(get_summary_files),
        check_copy_finished
    output:
        html = "Sample_{sample}/{sample}.pycoQC.html",
        json = "Sample_{sample}/{sample}.pycoQC.json"
    conda:
        "../env/pycoqc.yml"
    log:
        "logs/{sample}_pycoqc.log"
    threads:
        1
    shell:
        """
        pycoQC \
            --summary_file {input.summary_files} \
            --html_outfile {output.html} \
            --json_outfile {output.json} \
            > {log} 2>&1
        """

#____ GENOME MAPPING QC __________________________________________________________#

rule qualimap:
    input:
        "Sample_{sample}/{sample}.bam"
    output:
        directory('qc/qualimap/{sample}_genome')
    log:
        "logs/{sample}_qualimap.log"
    threads:
        8
    params:
        qualimap = config['apps']['qualimap']
    shell:
        """
        {params.qualimap} bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            -nt {threads} \
            -outdir {output} \
            --java-mem-size=12G \
            > {log} 2>&1
        """

#_____ cDNA SPLICED MAPPING QC _________________________________________________#

rule rna_qualimap:
    input:
        "Sample_{sample}/{sample}.spliced.bam"
    output:
        directory('qc/qualimap/{sample}_rna')
    log:
        "logs/{sample}_qualimap_rna.log"
    threads:
        2
    params:
        gtf = config['ref']['annotation'],
        qualimap = config['apps']['qualimap']
    shell:
        """
        {params.qualimap} rnaseq \
            -bam {input} \
            -gtf {params.gtf} \
            -outdir {output}\
            --java-mem-size=12G \
            > {log} 2>&1
        """

#____ ASSEMBLY QC _____________________________________________________________#

rule quast:
    input:
        files = expand("Sample_{s}/{s}.asm.{asm}.fasta",
            s=ID_samples,
            asm=config['assembly']['methods'])
    output:
        "qc/quast_results/report.tsv"
    conda:
        "../env/quast.yml"
    log:
        "logs/quast.log"
    params:
        ref=config['ref']['genome']
    threads:
        8
    shell:
        """
        quast \
            --threads {threads} \
            --no-sv             \
            --reference {params.ref} \
            --output-dir qc/quast_results 
            {input}
        """

#_____ MULTI QC  _____________________________________________________________#

qc_out = {
    'mapping' : expand("qc/qualimap/{s}_genome", s = ID_samples),
    'assembly' : ["qc/quast_results/report.tsv"],
    'cDNA' : 
        expand("qc/qualimap/{s}_rna", s = ID_samples) + 
        expand("Sample_{s}/{s}.summary.tsv", s = ID_samples),
    'pinfish_annotation' : [],
    'qc' : ["qc/per_run/run_multiqc_report.html"],
}

qc_out_selected = [qc_out[step] for step in config['steps']]

rule multiqc:
    input:
        expand("Sample_{sample}/{sample}.pycoQC.json", sample = ID_samples),
        [y for x in qc_out_selected for y in x]
    output:
        "qc/multiqc_report.html"
    threads:
        1
    params:
        multiqc = config['apps']['multiqc']
    shell:
        """
        {params.multiqc} \
            --force \
            --config config/multiqc.yml \
            --outdir qc\
            --ignore-symlinks \
            {input}
        """
