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
        library = "S", # U = unspecific, S = strand-specific
        args = "--noErrorModel" #Ignore indels/mismatches in aln (because ONT)
    shell:
        """
        salmon quant {params.args} \
            --threads {threads} \
            --libType {params.library} \
            --targets {input.trs} \
            --alignments {input.bam} \
            --output counts/{wildcards.sample}/ \
            > {log} 2>&1 
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
            -a {input.ann} \
            -o Sample_{wildcards.sample}/{wildcards.sample} \
            -L {input.bam} \
            > {log} 2>&1
        """