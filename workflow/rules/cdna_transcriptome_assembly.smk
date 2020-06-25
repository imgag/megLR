#_____ TRANSCRIPT ASSEMBLY ______________________________________________#

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
        {params.stringtie} -L -G {input.annot} -R -A -B -o {output} {input.bam} >{log} 2>&1
        """

#____ ISOFORM ANNOTATION ____________________________________________________#

#rule sqanti:
#    input:
#        "Sample_{sample}/{sample}.stringtie.gtf"
#    output:
#        "Sample_{sample}/{sample}"
