#____ cDNA READ ORIENTATION AND PRIMER TRIMMING _____________________________#

rule pychopper:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        fq = "Sample_{sample}/{sample}.full_length.fastq",
        unclass = "Sample_{sample}/{sample}.unclassified.fastq",
        rescued = "Sample_{sample}/{sample}.rescued.fastq",
        report = "qc/porechopper/{sample}_report.pdf",
        stats = "qc/porechopper/{sample}_stats.txt"
    log:
        "logs/{sample}_pychopper.log"
    conda:
        "../env/pychopper.yml"
    threads:
        8
    shell:
        """
        fq=$(mktemp)
        unpigz -d -c {input} > $fq 2> {log}
        cdna_classifier.py \
            -r {output.report} \
            -t {threads} \
            -S {output.stats} \
            -u {output.unclass} \
            -w {output.rescued} \
            $fq \
            {output.fq} \
            >> {log} 2>&1
        rm $fq
        """

rule map_to_transcriptome:
    input:
        "Sample_{sample}/{sample}.full_length.fastq"
    output:
        bam = "Sample_{sample}/{sample}.transcripts.bam"
    log:
        "logs/{sample}_minimap_transcripts.log"
    conda:
        "../env/minimap2.yml"
    threads:
        16
    params:
        transcriptome = config['ref']['cDNA'],
        opts = config['transcript']["minimap2_opts"],
        msec = config['transcript']["maximum_secondary"],
        psec = config['transcript']["secondary_score_ratio"]
    shell:
        """
        minimap2 --MD -ax map-ont -t {threads} \
            -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
            -p {params.psec} -N {params.msec} {params.opts} \
            {params.transcriptome} {input} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """
