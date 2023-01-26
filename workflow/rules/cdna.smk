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
        detect_umi = "-U" if config['cdna']['with_umi'] else "",
        kit = config['cdna']['primer_kit']
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
        "../env/umitools.yml"
    threads:
        1
    shell:
        """
        umitools \
            --stdin {input.bam} \
            --stout {output.bam} \
            --output-stats Sample_{wildcards.sample}/dedup_stats/{wildcards.sample} \
            --log {log} \
            --error {log} \

        samtools index {output.bam}
        """