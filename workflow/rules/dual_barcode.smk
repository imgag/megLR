#____ cDNA READ ORIENTATION AND PRIMER TRIMMING _____________________________#

rule demux_minibar:
    input:
        fq = "Sample_{sample}/{sample}.fastq.gz",
        bc = config['demux']['bc_index']
    output:
        "Sample_{sample}/demux/{sample}_unk.fastq"
    threads:
        2
    log:
        "logs/{sample}_minibar.log"
    params:
        minibar = config['demux']['minibar'],
        dist = config['demux']['bc_distance']
    shell:
        """
        {params.minibar} \
            -p {params.dist} \
            -S -F -P Sample_{wildcards.sample}/demux/{wildcards.sample}_ \
            {input.bc} {input.fq} > {log} 2>&1
        """

rule demux_minibar_stats:
    input:
        fq = "Sample_{sample}/{sample}.fastq.gz",
        bc = config['demux']['bc_index']
    output:
        "Sample_{sample}/{sample}.demux.txt"
    log:
        "logs/{sample}_minibar_stats.log"
    threads:
        2
    params:
        minibar = config['demux']['minibar'],
        dist = config['demux']['bc_distance']
    shell:
        """
        {params.minibar} \
            -D -p {params.dist} \
            {input.bc} {input.fq} > {output} 2> {log}
        """