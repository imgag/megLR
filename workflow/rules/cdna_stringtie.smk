#_____ TRANSCRIPTOME ASSEMBLY ______________________________________________#

rule stringtie:
    input:
        annot = config['ref']['annotation'],
        bam = "Sample_{sample}/{sample}.spliced.bam"
    output:
        gtf = "Sample_{sample}/{sample}.stringtie.gtf",
        abundance = "Sample_{sample}/{sample}.stringtie.abundance.tsv"

    log:
        "logs/{sample}_stringtie.log"
    threads:
        2
    params:
        stringtie = config['apps']['stringtie']
    conda:
        "../env/stringtie.yml"
    shell:
        """
        {params.stringtie} -L -R -v -B\
        -A {output.abundance} \
        -G {input.annot} \
        -o {output.gtf} \
        {input.bam} \
        >{log} 2>&1
        """
