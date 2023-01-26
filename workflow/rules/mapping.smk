#_____ MAPPING TO REFERENCE GENOME ___________________________________________#

rule map_genome_all:
    input: 
        genome = config['ref']['genome'],
        fq = "Sample_{sample}/{sample}.fastq.gz"
    output: 
        bam = "Sample_{sample}/{sample}.bam",
        bai= "Sample_{sample}/{sample}.bam.bai"
    conda:
        "../env/minimap2.yml"
    log:
        "logs/{sample}_minimap2.log"
    threads:
        20
    params:
        sample = "{sample}",
    shell:
        """
        minimap2 --MD -ax map-ont --eqx -t {threads} \
            -R "@RG\\tID:{params.sample}\\tSM:{params.sample}" \
            {input.genome} {input.fq} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """

rule map_genome_full_length:
    input:
        genome = config['ref']['genome'],
        fq = "Sample_{sample}/{sample}.full_length.fastq"
    output: 
        bam = "Sample_{sample}/{sample}.flts.bam",
        bai= "Sample_{sample}/{sample}.flts.bam.bai"
    conda:
        "../env/minimap2.yml"
    log:
        "logs/{sample}_minimap2.log"
    threads:
        20
    params:
        sample = "{sample}",
        ref = config['ref']['genome']
    shell:
        """
        minimap2 --MD -ax map-ont --eqx -t {threads} \
            -R "@RG\\tID:{params.sample}\\tSM:{params.sample}" \
            {input.genome} {input.fq} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """

#_____ SPLICE AWARE MAPPING WITH MINIMAP2  ___________________________________#

#TODO Adjust parameters;
# https://github.com/lh3/minimap2#map-long-mrnacdna-reads

rule map_genome_splice:
    input:
        genome = config['ref']['genome'],
        fq = "Sample_{sample}/{sample}.full_length.fastq"
    output:
        bam = "Sample_{sample}/{sample}.spliced.bam",
        bai= "Sample_{sample}/{sample}.spliced.bai"
    log:
        "logs/{sample}_minimap2_splice.log"
    conda:
        "../env/minimap2.yml"
    threads:
        config['max_threads']
    params:    
        min_mq = config['mapping']['min_qual']
    shell:
        """
        minimap2 -ax splice --eqx -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            {input.genome} {input.fq} 2> {log} \
             | samtools view -q {params.min_mq} -F 2304 -Sb \
             | samtools sort -@ 4 -m 4G  - -o {output.bam} >{log} 2>&1
        samtools index {output}
        """

#_____ MAPPING TO TRANSCRIPTOME __________________________________________#

rule map_to_transcriptome:
    input:
        trs = config['ref']['cDNA'],
        fq = "Sample_{sample}/{sample}.full_length.fastq"
    output:
        bam = "Sample_{sample}/{sample}.transcripts.bam"
    log:
        "logs/{sample}_minimap_transcripts.log"
    conda:
        "../env/minimap2.yml"
    threads:
        16
    params:
        opts = config['transcript']["minimap2_opts"],
        msec = config['transcript']["maximum_secondary"],
        psec = config['transcript']["secondary_score_ratio"]
    shell:
        """
        minimap2 --MD -ax map-ont --eqx -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            -p {params.psec} -N {params.msec} {params.opts} \
            {input.trs} {input.fq} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """