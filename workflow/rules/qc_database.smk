#_____ STORE QC VALUES IN CENTRAL LOCATION  ______________________________________________________#

rule copy_runqc:
    input:
        html = lambda wildcards: ["qc/pycoqc/per_run/" + r + ".pycoQC.html" for r in map_runs_folder[wildcards.run]],
        json = lambda wildcards: ["qc/pycoqc/per_run/" + r + ".pycoQC.json" for r in map_runs_folder[wildcards.run]],
    output:
        html = config['run_db_root'] + "/{run}/{run}.pycoQC.html",
        json = config['run_db_root'] + "/{run}/{run}.pycoQC.json"
    log:
        "logs/{run}_copy_runqc.log"
    shell:
        """
        cp -v {input.html} {output.html} > {log} 2>&1
        cp -v {input.json} {output.json} >> {log} 2>&1
        """

rule copy_run_report:
    input:
        pdf = lambda wildcards: [x for y in [glob(r + "/report_*.pdf") for r in map_runs_folder[wildcards.run]] for x in y],
        md = lambda wildcards:  [x for y in [glob(r + "/report_*.md") for r in map_runs_folder[wildcards.run]] for x in y]
    output:
        pdf = config['run_db_root'] + "/{run}/{run}.report.pdf",
        md = config['run_db_root'] + "/{run}/{run}.report.md"
    conda:
        "../env/poppler.yml"
    log:
        "logs/multiqc_copy_report_{run}.log"
    shell:
        """
        pdfunite {input.pdf} {output.pdf} > {log} 2>&1
        cat {input.md} > {output.md} 2> {log}
        """

rule copy_barcode_stats:
# TODO ptional copy the barcode stats (If available)
    input:
    
    output:
    
    shell:
        """

        """


rule all_multiqc:
    input:
        expand(config['run_db_root'] + "/{run}/{run}.pycoQC.html", run = ID_runs),
        expand(config['run_db_root'] + "/{run}/{run}.report.pdf", run = ID_runs)
    output:
        config['run_db_root'] + "/ont_runs_multiqc.html"
    log:
        "logs/multiqc_all_ont.log"
    params:
        multiqc = config['apps']['multiqc'],
        multiqc_config = srcdir('../../config/multiqc_config_db.yml'),
        multiqc_out = config['run_db_root'],
        multiqc_in = config['run_db_root']
    shell:
        """
        {params.multiqc} \
            --config  {params.multiqc_config} \
            --force \
            --outdir {params.multiqc_out} \
            --ignore-symlinks \
            {params.multiqc_in} > {log} 2>&1
        """


#rule copy_run_stats:
# TODO

#rule create_tables:
# TODO

#rule commit_new_files:
# TODO