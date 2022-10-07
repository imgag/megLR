#_____ STORE QC VALUES IN CENTRAL LOCATION  ______________________________________________________#

rule copy_runqc:
    input:
        html = "qc/pycoqc/per_run/{run}/{run}.pycoQC.html",
        json = "qc/pycoqc/per_run/{run}/{run}.pycoQC.json",
    output:
        html = config['run_db_root'] + "/runs/{run}/{run}.pycoQC.html",
        json = config['run_db_root'] + "/runs/{run}/{run}.pycoQC.json"
    log:
        "logs/{run}_copy_runqc.log"
    group:
        "qc_db"
    shell:
        """
        cp -v {input.html} {output.html} > {log} 2>&1
        cp -v {input.json} {output.json} >> {log} 2>&1
        """

rule copy_run_report:
    input:
        unpack(get_db_report)
    output:
        md = config['run_db_root'] + "/runs/{run}/{run}.report.md"
    log:
        "logs/multiqc_copy_report_{run}.log"
    group:
        "qc_db"
    shell:
        """
        cat {input.md} > {output.md} 2> {log}
        
        for path in {input.reports}
        do
            root="${{path%.*}}"
            ext="${{path#"$root"}}"
            cp ${{path}} $(dirname {output.md})/{wildcards.run}.report${{ext}}
        done
        """

rule copy_barcode_stats:
    input:
        ancient(get_db_barcode)
    output:
        tsv = config['run_db_root'] + "/runs/{run}/{run}.barcodes.tsv"
    log:
        "logs/{run}_copy_barcode_stats.log"
    group:
        "qc_db"
    shell:
        """
        cat  {input} > {output.tsv} 2>{log}       
        """

rule copy_mux_stats:
    input:
        csv = ancient(get_db_mux)
    output:
        mux = config['run_db_root'] + "/runs/{run}/{run}.mux.csv"
    log:
        "logs/{run}_copy_mux_stats.log"
    conda:
        "../env/R.yml"
    group:
        "qc_db"
    script: 
        "../scripts/db_mux_stats.R"

rule all_multiqc:
    input:
        expand(config['run_db_root'] + "/runs/{run}/{run}.pycoQC.json", run = ID_runs),
        expand(config['run_db_root'] + "/runs/{run}/{run}.report.md", run = ID_runs),
    output:
        config['run_db_root'] + "/ont_runs_multiqc.html"
    log:
        "logs/multiqc_all_ont.log"
    group:
        "qc_db"
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
            --filename ont_runs_multiqc.html \
            {params.multiqc_in} > {log} 2>&1
        """
        
rule update_table:
    input:
        expand(config['run_db_root'] + "/runs/{run}/{run}.pycoQC.json", run = ID_runs),
        expand(config['run_db_root'] + "/runs/{run}/{run}.report.md", run = ID_runs)
    output:
        csv = config['run_db_root'] + "/ont_runs.csv",
        xlsx = config['run_db_root'] + "/ont_runs.xlsx"
    conda:
        "../env/R.yml"
    group:
        "qc_db"
    params:
        db_root = config['run_db_root']
    script: 
        "../scripts/db_collect.R"

rule create_mux_plots:
    input:
        csv = config['run_db_root'] + "/runs/{run}/{run}.mux.csv"
    output:
        plot_scan = config['run_db_root'] + "/runs/{run}/{run}.pore_yield.png",
        plot_total = config['run_db_root'] + "/runs/{run}/{run}.pore_status.png"
    conda:
        "../env/R.yml"
    group:
        "qc_db"
    script: 
        "../scripts/db_mux_plot.R"

rule create_plots:
    input:
        csv = config['run_db_root'] + "/ont_runs.csv"
    output:
        timeline = config['run_db_root'] + "/plot/timeline.png",
        len_yield = config['run_db_root'] + "/plot/length_yield.png",
        device_pos = config['run_db_root'] + "/plot/device_position.png",
    conda:
        "../env/R.yml"
    group:
        "qc_db"
    script: 
        "../scripts/db_plot.R"
