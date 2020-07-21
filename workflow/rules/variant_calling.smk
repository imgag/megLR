#____ VARIANT CALLING WITH MEDAKA _______________________________________________________#

rule medaka_variants:
    input:
        bam = rules.map_genome_all.output.bam,
        ref = config['ref']['genome']
    output:
        vcf = "variant_calling/{sample}/{chr}/round_1.vcf"
    conda:
        "../env/medaka.yml"
    threads:
        2
    log:
        "logs/{sample}_{chr}_medaka.log"
    params:
        model_snp = config['vc']['model_initial'],
        model_final = config['vc']['model_final']
    shell:
        """
        export OMP_NUM_THREADS={threads}
        medaka_variant -d \
            -f {input.ref} \
            -i {input.bam} \
            -r {wildcards.chr} \
            -p  \
            -o variant_calling/{wildcards.sample}/{wildcards.chr} \
            -t {threads} \
            -s {params.model_snp} \
            -m {params.model_final} \
            >{log} 2>&1
        """

rule combine_vcf:
    input:
        expand("variant_calling/{{sample}}/{C}/round_1.vcf",  C = get_chromosomes())
    output:
        "variant_calling/{sample}/{sample}_var.vcf"
    conda:
        "../env/bcftools.yml"
    threads:
        1
    log:
        "logs/{sample}_combinevcf.log"
    conda:
        "env/bcftools.yml"
    shell:
        """
        bcftools concat $(echo '{input}' | sort -V) | bcftools annotate > {output} 2> {log}
        """

rule process_vcf:
    input:
        rules.combine_vcf.output
    output:
        "Sample_{sample}/{sample}_var.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/{sample}_processvcf.log"
    shell:
        """
        bgzip -c {input} > {output} 2> {log}
        tabix {output} 2>>{log}
        """

#____ VARIANT BENCHMARK TO REFERENCE ___________________________________________________#

rule benchmark_happy:
    input:
        truth = config['ref']['vc_benchmark'],
        query = rules.process_vcf.output,
        ref = config['ref']['genome'],
        target = config['ref']['target_region'],
        conf_region = config['ref']['vc_confidence_region']
    output:
        "qc/happy/{sample}.summary.csv"
    log:
        "logs/{sample}_happy.log"
    conda:
        "../env/happy.yml"
    threads:
        8
    shell:
        """
        hap.py \
            -r {input.ref} \
            -T {input.target} \
            -f {input.conf_region} \
            --threads {threads} \
            -o qc/happy/{wildcards.sample}\
            {input.truth} \
            {input.query} \
            2> {log}
        """
