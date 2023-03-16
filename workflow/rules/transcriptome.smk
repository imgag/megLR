#_____ STRINGTIE RANSCRIPTOME ASSEMBLY ______________________________________________#
# We run Stringtie in stranded mode since the bam file after
# Pychopper contains correct strand informations

rule stringtie:
    input:
        bam = get_cdna_bam  
    output:
        gtf = "Sample_{sample}/stringtie/{sample}.stringtie.gtf",
        abundance = "Sample_{sample}/stringtie/{sample}.stringtie.abundance.tsv"
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

rule copy_stringtie:
    input:
        "Sample_{sample}/stringtie/{sample}.stringtie.gtf"
    output:
        "Sample_{sample}/{sample}.stringtie.gtf"
    threads:
        1
    shell:
        """
        cp {input} {output}
        """

rule create_stringtie_transcriptome:
    input:
        gtf = "Sample_{sample}/stringtie/{sample}.stringtie.gtf",
        ref = config['ref']['genome']
    output:
        "Sample_{sample}/{sample}.stringtie.fasta"
    log:
        "logs/{sample}_create_stringtie_transcriptome.log"
    threads:
        1
    conda:
        "../env/gffread.yml"
    shell:
        """
        gffread -w {output} -g {input.ref} {input.gtf} > {log} 2>&1
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
        """

rule copy_flair_results:
    input:
        gtf = "Sample_{sample}/flair/{sample}.isoforms.gtf",
        fa = "Sample_{sample}/flair/{sample}.isoforms.fa"
    output:
        gtf = "Sample_{sample}/{sample}.flair.gtf",
        fa = "Sample_{sample}/{sample}.flair.fasta"
    shell:
        """
        cp {input.gtf} {output.gtf}
        cp {input.fa} {output.fa}
        """

# Other unimplemented rules:
#rule flairDiffExp:
#rule flairDiffSpice:
