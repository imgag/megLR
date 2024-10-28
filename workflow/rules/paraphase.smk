rule paraphase:
    input:
        bam="Sample_{sample}/{sample}.bam",
        ref=lambda wc: config['ref']['genome'],
    output:
        bam="Sample_{sample}/paraphase/{sample}.paraphase.bam",
    threads: 4
    conda:
        "../env/paraphase.yml"
    params:
        genome = config['paraphase']['genome_build'],
        region = lambda wc: "-g " + config['paraphase']['genes']
        if config['paraphase']['genes']
        else "",
    log:
        "logs/{sample}.paraphase.log",
    shell:
        """
        paraphase \
            -b {input.bam} \
            -o $(dirname {output.bam}) \
            -p {wildcards.sample} \
            -r {input.ref} \
            -t {threads} \
            --genome {params.genome} \
            {params.region} \
            >{log} 2>{log}
        """
