#_____ MAPPING TO REFERENCE GENOME ___________________________________________#

rule map_genome_full:
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
        ref = config['ref']['genome']
    shell:
        """
        minimap2 --MD -ax map-ont -t {threads} \
            -R "@RG\\tID:{params.sample}\\tSM:{params.sample}" \
            {input.genome} {input.ref} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """

#_____ SPLICE AWARE MAPPING WITH MINIMAP2  ___________________________________#

rule map_genome_splice:
    input:
        genome = config['ref']['genome'],
        fq = "Sample_{sample}/{sample}.fastq.gz"
        #TODO Better to use full length transcripts only?
    output:
        bam = "Sample_{sample}/{sample}.spliced.bam",
        bai= "Sample_{sample}/{sample}.bam.bai"
    log:
        "logs/{sample}_minimap2_splice.log"
    conda:
        "../env/minimap2.yml"
    threads:
        config['sys']['max_threads']
    shell:
        """
        minimap2 -ax splice -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            {input.genome} {input.fq} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output} - >{log} 2>&1
        samtools index {output}
        """

#_____ MAPPING TO TRANSCRIPTOME __________________________________________#

rule map_to_transcriptome:
    input:
        trs = config['ref']['cDNA'],
        fq = "Sample_{sample}/{sample}.fastq.gz"
        #rules.pychopper.output.fq TODO use pychopped output (or both?)
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
        minimap2 --MD -ax map-ont -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            -p {params.psec} -N {params.msec} {params.opts} \
            {input.trs} {input.fq} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """

#____ MAPPING QC _____________________________________________________________#

rule qualimap:
    input:
        "Sample_{sample}/{sample}.bam"
    output:
        report = "qc/qualimap/{sample}_genome/qualimap_report.pdf",
        stats = "qc/qualimap/{sample}_genome/genome_results.txt"
    log:
        "logs/{sample}_qualimap.log"
    threads:
        8
    params:
        qualimap = config['apps']['qualimap']
    shell:
        """
        {params.qualimap} bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            -nt {threads} \
            -outdir qc/qualimap/{wildcards.sample}_genome \
            -outfile qualimap_report \
            > {log} 2>&1
        """
