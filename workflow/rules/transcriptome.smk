#_____ STRINGTIE RANSCRIPTOME ASSEMBLY ______________________________________________#
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

#_____ FLAIR ISOFORM ANALYSIS ______________________________________________#

rule bamtobed12:
    input:
        bam = get_cdna_bam
    output:
        bed = "Sample_{sample}/{sample}.spliced.bed12"
    conda:
        "../env/flair.yml"
    threads:
        1
    log:
        "logs/{sample}_flair_bamtobed12.log"
    shell:
        """
        b2b=$(dirname $(which flair.py))/bin/bam2Bed12.py
        $b2b \
            --input_bam {input.bam} \
            > {output.bed} 2> {log}
        """

rule flairCorrect:
    input:
        bed = rules.bamtobed12.output.bed,
        genome = config['ref']['genome'],
        annot = config['ref']['annotation'],
        chromsize = config['ref']['target_region']
    output:
        bed_corrected = "Sample_{sample}/flair/{sample}_all_corrected.bed",
        psl_corrected = "Sample_{sample}/flair/{sample}_all_corrected.psl"
    conda:
        "../env/flair.yml"
    threads:
        8
    log:
        "logs/{sample}_flair_correct.log"
    shell:
        """
        flair.py correct \
            -q {input.bed} \
            -f {input.annot} \
            -g {input.genome} \
            -o Sample_{wildcards.sample}/flair/{wildcards.sample} \
            --nvrna \
            --threads {threads} \
            >{log} 2>&1

        bed2psl=$(dirname $(which flair.py))/bin/bed_to_psl.py
        $bed2psl {input.chromsize} {output.bed_corrected} {output.psl_corrected}
        """

rule flairCollapse:
    input:
        psl = rules.flairCorrect.output.psl_corrected,
        fastq = rules.pychopper.output.fq,
        genome = config['ref']['genome'],
        annot = config['ref']['annotation']
    output:
        gtf = "Sample_{sample}/flair/{sample}.isoforms.gtf",
        gtf_main = "Sample_{sample}/{sample}.flair.gtf",
        fasta = "Sample_{sample}/flair/{sample}.isoforms.fa",
    conda:
        "../env/flair.yml"
    threads:
        8
    log:
        "logs/{sample}_flair_collapse.log"
    params:
        trust_ends = "--trust_ends" if config['flair']['trust_ends'] else ""
    shell:
        """
        flair.py collapse \
            --genome {input.genome} \
            --query {input.psl} \
            --reads {input.fastq} \
            --threads {threads} \
            --gtf {input.annot} \
            {params.trust_ends} \
            --generate_map \
            --output Sample_{wildcards.sample}/flair/{wildcards.sample} \
            >{log} 2>&1
        cp {output.gtf} {output.gtf_main}
        """

rule flairQuantify:
    input:
        fa = rules.flairCollapse.output.fasta,
        manifest = config['flair']['manifest']
    output:
        mat = "Sample_{sample}/flair/{sample}.counts_matrix.tsv"
    conda:
        "../env/flair.yml"
    threads:
        20
    log:
        "logs/{sample}_flair_quantify.log"
    params:
        trust_ends = "--trust_ends" if config['flair']['trust_ends'] else "",
        qual = config['flair']['min_qual']
    shell:
        """ 
        flair.py quantify \
            --quality {params.qual} \
            {params.trust_ends} \
            --output {output} \
            --tpm \
            --reads_manifest {input.manifest} \
            --isoforms {input.fa} \
            >{log} 2>&1
        """

# Other unimplemented rules:
#rule flairDiffExp:
#rule flairDiffSpice:
