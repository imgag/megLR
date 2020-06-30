#_____ TRANSCRIPTOME ASSEMBLY ______________________________________________#

rule stringtie:
    input:
        annot = config['ref']['annotation'],
        bam = "Sample_{sample}/{sample}.spliced.bam"
    output:
        "Sample_{sample}/{sample}.stringtie.gtf"
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
        {params.stringtie} -L -R -A -B -v\
        -G {input.annot} \
        -o {output} \
        {input.bam} \
        >{log} 2>&1
        """
