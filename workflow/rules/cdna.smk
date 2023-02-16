#____ cDNA READ ORIENTATION AND PRIMER TRIMMING _____________________________#

rule pychopper:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        fq = "Sample_{sample}/{sample}.full_length.fastq",
        unclass = "Sample_{sample}/{sample}.unclassified.fastq",
        rescued = "Sample_{sample}/{sample}.rescued.fastq",
        report = "qc/pychopper/{sample}_report.pdf",
        stats = "qc/pychopper/{sample}_stats.txt"
    log:
        "logs/{sample}_pychopper.log"
    conda:
        "../env/pychopper.yml"
    threads:
        24
    params:
        detect_umi = "-U -y" if config['cdna']['with_umi'] else "",
        kit = config['cdna']['primer_kit']
    shell:
        """
        fq=$(mktemp)
        unpigz -d -c {input} > $fq 2> {log}
        pychopper \
            -r {output.report} \
            -t {threads} \
            -S {output.stats} \
            -u {output.unclass} \
            -w {output.rescued} \
            -k {params.kit} \
            {params.detect_umi} $fq \
            {output.fq} \
            >> {log} 2>&1
        rm $fq
        """

#_____ DEDUPLICATE READS WITH UMIS __________________________________________#

rule deduplicate_umitools:
    input:
        bam = "Sample_{sample}/{sample}.spliced.bam"
    output:
        bam = "Sample_{sample}/{sample}.spliced.dedup.bam",
        bai= "Sample_{sample}/{sample}.spliced.dedup.bai",
        stats= "qc/umitools_dedup/{sample}_stats_per_umi_per.tsv"
    log:
        "logs/{sample}_dedup_umitools.log"
    conda:
        "../env/umi_tools.yml"
    threads:
        1
    shell:
        """
        umi_tools dedup \
            --stdin {input.bam} \
            --stdout {output.bam} \
            --output-stats qc/dedup_stats/{wildcards.sample} \
            --extract-umi-method tag \
            --umi-tag RX \
            --log {log} \
            --error {log} \

        samtools index {output.bam}
        """