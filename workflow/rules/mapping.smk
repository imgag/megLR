#_____ MAPPING TO REFERENCE GENOME ___________________________________________#

rule ref_alignment:
    input: 
        "Sample_{sample}/{sample}.fastq.gz"
    output: 
        bam = "Sample_{sample}/{sample}.bam",
        bai= "Sample_{sample}/{sample}.bam.bai"
    conda:
        "env/minimap2.yaml"
    log:
        "logs/{sample}_minimap2.log"
    threads:
        config['sys']['max_threads']
    params:
        sample = "{sample}",
        ref = config['ref']['genome']
    shell:
        """
        minimap2 --MD -ax map-ont -t {threads} \
            -R "@RG\\tID:{params.sample}\\tSM:{params.sample}" \
            {params.ref} {input} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output} - >{log} 2>&1
        samtools index {output}
        """

#____ MAPPING QC _____________________________________________________________#

rule samtools_stats:
    input:
        "Sample_{sample}/{sample}.bam"
    output:
        stats = "Sample_{sample}/{sample}.bam.stats",
        flagstats = "Sample_{sample}/{sample}.bam.flagstats",
    log:
        "logs/{sample}_samtools_stats.log"
    threads:
        4
    shell:
        """
        samtools stats -@{threads} > {output.stats}
        samtoools flagstat -@{threads} > {output.flagstats}
        """ 

#____ COVERAGE DEPTH ANALYSIS ________________________________________________#
