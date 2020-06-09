rule splice_mapping:
    input: 
        "Sample_{sample}/{sample}.fastq.gz"
    output: 
        "Sample_{sample}/{sample}.spliced.bam"
    log:
        "logs/{sample}_minimap2.log"
    threads:
        12
    params:
        sample = "{sample}",
        ref = "/mnt/share/data/genomes/GRCh38.fa"
    shell:
        """
        minimap2 -ax splice -t {threads} \
            -R "@RG\\tID:{params.sample}\\tSM:{params.sample}" \
            {params.ref} {input} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output} - 2>> {log}
        samtools index {output}
        """

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
        annot = "Homo_sapiens.GRCh38.85.gtf"
    conda:
        "env/stringtie.yml"
    shell:
        """
        stringtie -L -G -A -B {params.annot} -o {output} {input}
        """
