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

#_____ cDNA ALIGNMENT QC _______________________________________________________#

rule rna_qualimap:
    input:
        "Sample_{sample}/{sample}.spliced.bam"
    output:
        report = "qc/qualimap/{sample}_rna/qualimap_report.pdf",
        stats = "qc/qualimap/{sample}_rna/rna_results.txt"
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
            -outdir qc/qualimap/{wildcards.sample}_rna/ \
            -outformat PDF:HTML \
            --java-mem-size=12G
        """

#_____ MULTI QC  _____________________________________________________________#

qc_out = {
    'mapping' : expand("qc/qualimap/{s}_genome/genome_results.txt", s = ID_samples),
    'assembly' : ["qc/quast_results/report.tsv"],
    'cDNA' : expand("qc/qualimap/{s}_rna/rna_results.txt", s = ID_samples),
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
