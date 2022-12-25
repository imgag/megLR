rule trgt:
    input:
        bam = "Sample_{sample}/{sample}.bam"
    output:
        vcf = "trgt/{sample}/{sample}.vcf.gz",
        bam = "trgt/{sample}/{sample}.spanning.bam"
    params:
        trgt = config['repeat_expansion']['trgt'],
        ref = config['ref']['genome'],
        repeats = config['repeat_expansion']['loci']
    log:
        "logs/{sample}_trgt.log"
    shell:
        """
        {params.trgt} \
            --genome {params.ref} \
            --repeats {params.repeats} \
            --reads {input.bam} \
            --output-prefix trgt/{wildcards.sample}/{wildcards.sample} \
            >{log} 2>&1
        """

rule copy_vcf:
    input:
        rules.trgt.output.vcf 
    output:
        "Sample_{sample}/{sample}.repeats.vcf.gz"
    shell:
        "cp {input} {output}"

rule index_spanning_reads:
    input:
        bam = rules.trgt.output.bam
    output:
        bai = "trgt/{sample}/{sample}.spanning.bam.bai"
    conda:
        "../env/samtools.yml"
    shell:
        """
        samtools index {input}
        """

rule trvz:
    input:
        bam  = rules.trgt.output.bam,
        bai = rules.index_spanning_reads.output.bai,
        vcf = rules.trgt.output.vcf
    output:
        image = "trgt/{sample}/plots/{sample}_{loci}.png"
    params:
        trvz = config['repeat_expansion']['trvz'],
        ref = config['ref']['genome'],
        repeats = config['repeat_expansion']['loci']
    log:
        "logs/{sample}_trvz_{loci}.log"
    shell:
        """
        {params.trvz} \
            --genome {params.ref} \
            --repeats {params.loci} \
            --vcf {input.vcf} \
            --repeat_id {wildcards.loci} \
            --spanning-reads {input.bam} \
            --image {output} \
            >{log} 2>&1
        """