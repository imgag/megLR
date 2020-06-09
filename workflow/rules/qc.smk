#_____ RUN READ QC  __________________________________________________________#

rule run_pycoqc:
    input:
        "run_data/{run}/copy_finished"
    output:
        html = "qc/per_run/{run}.pycoQC.html",
        json = "qc/per_run/{run}.pycoQC.json"
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
            --json_outfile {output.json}
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
            --config config/multiqc.yaml \
            --outdir qc/per_run \
            --fileame run_multiqc_report \
            > {log} 2>> {log}
        """

#_____ SAMPLE READ QC  _________________________________________________________#

rule sample_pycoqc:
    input:
        unpack(get_summary_files),
        check_copy_finished
    output:
        html = "Sample_{sample}/{sample}.pycoQC.html",
        json = "Sample_{sample}/{sample}.pycoQC.json",
        qc_update = touch("qc/updates/pycoqc_{sample}")
    conda:
        "../env/pycoqc.yml"
    threads:
        1
    shell:
        """
        pycoQC \
            --summary_file {input.summary_files} \
            --html_outfile {output.html} \
            --json_outfile {output.json}
        """

#_____ MULTI QC  _____________________________________________________________#

rule multiqc:
    input:
        expand("Sample_{sample}/{sample}.pycoQC.json", sample = ID_samples)
    output:
        "multiqc_report.html"
    threads:
        1
    params:
        multiqc = config['apps']['multiqc']
    shell:
        """
        {params.multiqc} \
            --force \
            --config config/multiqc.yaml \
            --outdir qc\
            --ignore-symlinks Sample_*
        """