rule bamtobed12:
    input:
        bam = rules.map_genome_splice.output.bam
    output:
        bed = "Sample_{sample}/{sample}.splice.bed12"
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
        bed_corrected = "flair/{sample}_all_corrected.bed",
        psl_corrected = "flair/{sample}_all_corrected.psl"
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
            -o flair/{wildcards.sample} \
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
        annot = config['ref']['annotation']
    output:
        gtf = "flair/{sample}.isoforms.gtf",
        gtf_main = "Sample_{sample}/{sample}.flair.gtf",
        fasta = "flair/{sample}.isoforms.fa",
    conda:
        "../env/flair.yml"
    threads:
        8
    log:
        "logs/{sample}_flair_correct.log"
    params:
        trust_ends = "--trust_ends" if config['flair']['trust_ends'] else ""
    shell:
        """
        flair.py collapse \
            --gtf {input.annot} \
            --reads {input.fastq} \
            --query P
            --threads {threads} \
            {params.trust_ends} \
            --generate_map \
            --output flair/{wildcards.sample} \
            >{log} 2>&1
        cp {output.gtf} {output.gtf_main}
        """

rule flairQuantify:
    input:
        fa = rules.flairCollapse.output.fasta,
        manifest = config['flair']['manifest']
    output:
        mat = "flair/{sample}.counts_matrix.tsv"
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

#rule flairDiffExp:

#rule flairDiffSpice:
