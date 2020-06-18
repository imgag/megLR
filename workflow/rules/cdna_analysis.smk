#_____ FEATURE COUNT (TRANSCRIPT LEVEL) ____________________________________________#

rule quant_salmon:
    input:
        bam = rules.map_to_transcriptome.output.bam,
        trs = config['ref']['cDNA']
    output:
        counts = "counts/{sample}/quant.sf",
        stats = "counts/{sample}/aux_info/meta_info.json"
    conda:
        "../env/salmon.yml"
    log: 
        "logs/{sample}_salmon.log"
    threads:
        10
    params:
        library = "U", # U = unspecific, S = strand-specific
        args = "--noErrorModel" #Ignore indels/mismatches in aln (because ONT)
    shell:
        """
        salmon quant {params.args} \
            --threads {threads} \
            --libType {params.library} \
            --targets {input.trs} \
            --alignments {input.bam} \
            --output counts/{wildcards.sample}/ 
        """

#_____ FEATURE COUNT (GENE LEVEL) _____________________________________________#

rule quant_genes:
    input:
        bam = "Sample_{sample}/{sample}.spliced.bam",
        ann = config['ref']['annotation']
    output:
        counts = "Sample_{sample}/{sample}.tsv",
        stats = "Sample_{sample}/{sample}.summary.tsv"
    conda:
        "../env/subread.yml"
    log:
        "logs/{sample}_featureCounts.log"
    threads:
        1
    shell:
        """
        featureCounts \
            -a {input.ann}  \
            -o Sample_{wildcards.sample}/{wildcards.sample} \ 
            -L > {log} 2>&1
        """

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
        annot = config['ref']['annotation'],
        stringtie = config['apps']['stringtie']
    conda:
        "../env/stringtie.yml"
    shell:
        """
        {params.stringtie} -L -G {params.annot} -A -B -o {output} {input} >{log} 2>&1
        """

#____ ISOFORM ANNOTATION ____________________________________________________#

rule sqanti:
    input:
        "Sample_{sample}/{sample}.stringtie.gtf"
    output:
        "Sample_{sample}/{sample}"
