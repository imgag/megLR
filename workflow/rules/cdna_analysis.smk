#_____ SPLICE AWARE MAPPING WITH MINIMAP2  ___________________________________#

rule splice_mapping:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        "Sample_{sample}/{sample}.spliced.bam"
    log:
        "logs/{sample}_minimap2_splice.log"
    conda:
        "../env/minimap2.yml"
    threads:
        config['sys']['max_threads']
    params:
        sample = "{sample}",
        ref = config['ref']['genome']
    shell:
        """
        minimap2 -ax splice -t {threads} \
            -R "@RG\\tID:{params.sample}\\tSM:{params.sample}" \
            {params.ref} {input} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output} - >{log} 2>&1
        samtools index {output}
        """

#_____ FEATURE COUNT ANALYSIS _________________________________________________#

rule featureCount:
    input:
        "Sample_{sample}/{sample}.spliced.bam"
    output:
        "Sample_{sample}/{sample}"

#_____ TRANSCRIPT RECONSTRUCTION ______________________________________________#

rule stringtie:
    input:
        "Sample_{sample}/{sample}.spliced.bam"
    output:
        "Sample_{sample}/{sample}.stringtie.gtf"
    log:
        "logs/{sample}_stringtie.log"
    threads:
        2
    params:
        annot = config['ref']['annotation']
    conda:
        "../env/stringtie.yml"
    shell:
        """
        stringtie -L -G {params.annot} -A -B -o {output} {input} >{log} 2>&1
        """

#____ ISOFORM ANNOTATION ____________________________________________________#

rule sqanti:
    input:
        "Sample_{sample}/{sample}.stringtie.gtf"
    output:
        "Sample_{sample}/{sample}"
