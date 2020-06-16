#_____ MAPPING TO REFERENCE GENOME ___________________________________________#

rule ref_alignment:
    input: 
        "Sample_{sample}/{sample}.fastq.gz"
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
            {params.ref} {input} 2> {log} \
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
