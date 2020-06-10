#____ ASSEMBLY _______________________________________________________________#

rule wtdbg2:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        "Sample_{sample}/{sample}.asm.wtdbg.fasta"
    conda:
        "env/wtdbg2.yml"
    log:
        "logs/{sample}_wtdbg.log"
    params:
        g=config['assembly']['genome_size']
    threads:
        config['sys']['max_threads']
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        mkdir -p $outf
        wtdbg \
            -x nanopore     \
            -g {params.g}   \
            -t {threads}    \
            -i {input}      \
            -fo $outf/wtdbg \
            > {log} 2>&1
        cp $outf/wtdbg.ctg.fa {output}
        """

rule flye:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        "Sample_{sample}/{sample}.asm.flye.fasta"
    conda:
        "env/flye.yml"
    log:
        "logs/{sample}_flye.log"
    params:
        g=config['assembly']['genome_size']
    threads:
        config['sys']['max_threads']
    shell:
        """
        flye \
            --nano-raw      {input}     \
            --genome-size   {params.g}  \
            --threads       {threads}   \
            --out-dir       assembly/{wildcards.sample}_flye/} \
            > {log} 2>&1
        cp assembly/{wildcards.sample}_flye/assembly.fasta {output}
        """

#____ ASSEMBLY QC _____________________________________________________________#

rule quast:
    input:
        files = expand("Sample_{s}/{s}.asm.{asm}.fasta",
            s=ID_samples,
            asm=config['assembly']['methods']) ,
        labels = expand("{s}_{asm}",
            s=ID_samples,
            asm=config['assembly']['methods'])
    output:
        "qc/quast_results/report.tsv"
    conda:
        "env/quast.yml"
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
            --output-dir qc/quast_results \
            --labels {input.labels}
            {input}
        """