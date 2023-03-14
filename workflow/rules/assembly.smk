#____ ASSEMBLY USING WTDBG2 ________________________________________________________#

rule wtdbg2:
    input:
        unpack(get_assembly_input)
    output:
        "assembly/{sample}_wtdbg/wtdbg.ctg.lay.gz"
    conda:
        "../env/wtdbg2.yml"
    log:
        "logs/{sample}_wtdbg.log"
    params:
        g=config['assembly']['genome_size']
    threads:
        config['max_threads']
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        mkdir -p $outf
        wtdbg2 \
            -x nanopore     \
            -g {params.g}   \
            -t {threads}    \
            -i {input.fq}      \
            -fo $outf/wtdbg \
            > {log} 2>&1
        """

rule wtdbg2_consensus:
    input:
        "assembly/{sample}_wtdbg/wtdbg.ctg.lay.gz"
    output:
        "assembly/{sample}_wtdbg/wtdbg.raw.fa"
    conda:
        "../env/wtdbg2.yml"
    log:
        "logs/{sample}_wtdbg_consensus.log"
    threads:
        config['max_threads']
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        wtpoa-cns \
            -t {threads}\
            -i {input}  \
            -fo {output}\
            > {log} 2>&1
        """

rule wtdbg2_polishing:
    input:
        asm = "assembly/{sample}_wtdbg/wtdbg.raw.fa",
        reads = "Sample_{sample}/{sample}.fastq.gz"
    output:
        "Sample_{sample}/{sample}.asm.wtdbg.fasta"
    conda:
        "../env/wtdbg2.yml"
    log:
        "logs/{sample}_wtdbg_consensus.log"
    threads:
        config['max_threads']
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        minimap2 -t{threads} -ax map-ont -r2k {input.asm} {input.reads} \
            | samtools sort -@4 > $outf/wtdbg.polishing.bam
        samtools view -F0x900 $outf/wtdbg.polishing.bam \
            | wtpoa-cns -t {threads} -d {input.asm} -i - -fo $outf/wtdbg.pol.fa \
            > {log} 2>&1
        cp $outf/wtdbg.pol.fa {output}
        rm $outf/wtdbg.polishing.bam
        """

#____ ASSEMBLY USING FLYE ________________________________________________________#

rule flye:
    input:
        unpack(get_assembly_input)
    output:
        "assembly/{sample}_flye/assembly.fasta"
    conda:
        "../env/flye.yml"
    log:
        "logs/{sample}_flye.log"
    params:
        g=config['assembly']['genome_size'],
        cov = "--asm-coverage "+str(config['assembly']['max_flye_cov']) if config['assembly']['max_flye_cov'] else ""
    threads:
        config['max_threads']
    shell:
        """
        flye \
            --nano-hq {input.fq} \
            --genome-size {params.g} \
            --threads {threads} \
            -o assembly/{wildcards.sample}_flye \
            {params.cov} > {log} 2>&1
        """

rule cp_flye:
    input:
        "assembly/{sample}_flye/assembly.fasta"
    output:
        "Sample_{sample}/{sample}.asm.flye.fasta"
    threads:
        1
    shell:
        """
        cp {input} {output}
        """

rule metaflye:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        "assembly/{sample}_metaflye/assembly.fasta"
    conda:
        "../env/flye.yml"
    log:
        "logs/{sample}_metaflye.log"
    threads:
        config['max_threads']
    shell:
        """
        flye --nano-raw {input} --meta --threads {threads} -o assembly/{wildcards.sample}_metaflye > {log} 2>&1
        """

rule cp_metaflye:
    input:
        "assembly/{sample}_metaflye/assembly.fasta"
    output:
        "Sample_{sample}/{sample}.asm.metaflye.fasta"
    threads:
        1
    shell:
        """
        cp {input} {output}
        """

#____ HIGHCOV GENOME ASSEMBLY ________________________________________________________#

rule filtlong:
    input:
        fq = "Sample_{sample}/{sample}.fastq.gz",
        qc = "qc/pycoqc/per_sample/{sample}.pycoQC.json"
    output:
        "assembly/{sample}_subsampled.fastq.gz"
    conda:
        "../env/filtlong.yml"
    log:
        "logs/{sample}_filtlong.log"
    threads:
        1
    params:
        n_bases = human2bytes(config['assembly']['genome_size'].upper())*config['assembly']['max_target_coverage']
    shell:
        """
        filtlong \
            --target_bases {params.n_bases} \
            --length_weight 2 \
            {input.fq} | gzip -c > {output} 2>{log}
        """
