# ____ VARIANT BENCHMARK TO REFERENCE ___________________________________________________#


rule benchmark_happy:
    input:
        truth=config["ref"]["vc_benchmark"],
        query=rules.process_vcf.output,
        ref=config["ref"]["genome"],
        target=config["ref"]["target_region"],
        conf_region=config["ref"]["vc_confidence_region"],
    output:
        "qc/happy/{sample}.summary.csv",
    log:
        "logs/{sample}_happy.log",
    conda:
        "../env/happy.yml"
    threads: 8
    shell:
        """
        hap.py \
            -r {input.ref} \
            -T {input.target} \
            -f {input.conf_region} \
            --threads {threads} \
            -o qc/happy/{wildcards.sample}\
            {input.truth} \
            {input.query} \
            2> {log}
        """
