#_____ TRANSCRIPTOME ASSEMBLY ______________________________________________#
# We run Stringtie in stranded mode since the bam file after
# Pychopper contains correct strand informations

rule stringtie:
    input:
        bam = get_cdna_bam  
    output:
        gtf = "Sample_{sample}/{sample}.stringtie.gtf",
        abundance = "Sample_{sample}/{sample}.stringtie.abundance.tsv"
    log:
        "logs/{sample}_stringtie.log"
    threads:
        2
    conda:
        "../env/stringtie.yml"
    params:
        annot = config['ref']['annotation']
    # Stringtie Parameters:
    # -L : Input consists of long reads
    # -B : Output coverage data (.ctab) in Ballgown Format
    # -G : Reference guided (requires input .GTF)
    # -A : Output gene Abundances
    shell:
        """
        stringtie -L -R -v -B \
        --rf \
        -A {output.abundance} \
        -G {params.annot} \
        -o {output.gtf} \
        {input.bam} \
        >{log} 2>&1
        """
