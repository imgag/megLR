checkpoint paraphase:
    input:
        bam="Sample_{sample}/{sample}.bam",
        ref=lambda wc: config['ref']['genome'],
    output:
        bam="Sample_{sample}/paraphase/{sample}.paraphase.bam"
    threads: 4
    conda:
        "../env/paraphase.yml"
    params:
        genome = config['paraphase']['genome_build'],
        region = lambda wc: "-g " + config['paraphase']['genes']
        if config['paraphase']['genes']
        else "",
    log:
        "logs/{sample}.paraphase.log",
    shell:
        """
        paraphase \
            -b {input.bam} \
            -o $(dirname {output.bam}) \
            -p {wildcards.sample} \
            -r {input.ref} \
            -t {threads} \
            --genome {params.genome} \
            {params.region} \
            >{log} 2>{log}
        """

rule paraphase_add_sample_header:
    input:
        vcf="Sample_{sample}/paraphase/{sample}_paraphase_vcfs/{sample}_{gene}.vcf",
    output:
        vcf="Sample_{sample}/paraphase/{sample}_paraphase_vcfs/{sample}_{gene}.reheader.vcf"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/{sample}.{gene}.paraphase_add_sample_header.log"
    shell:
        """
        samples=$(bcftools query -l {input.vcf} | tr '\n' ',' | head -c -1) 2>{log}
        bcftools annotate -h "##SAMPLE=<${{samples}}>" > {output.vcf} 2>>{log}
        """

rule paraphase_annotate_vcf:
    input:
        vcf="Sample_{sample}/paraphase/{sample}_paraphase_vcfs/{sample}_{gene}.reheader.vcf"
    output:
        vcf="Sample_{sample}/paraphase/{sample}_paraphase_vcfs/{sample}_{gene}_var_annotated.vcf"
    log:
        "logs/{sample}.{gene}.paraphase_annotate_vcf.log"
    params:
        megsap = config['apps']['megsap'],
        system = config['paraphase']['annotation_system']
    shell:
        """
        {params.megsap}/src/Pipelines/annotate.php \
            -vcf {input} \
            -multi \
            -out_folder $dirname({output.vcf}) \
            -out_bane {wildcards.sample}_{wildcards.gene} \
            -system {params.system} > {log} 2>{log}
        """

def find_paraphase_vcf(wc):
    paraphase_bam = checkpoints.paraphase.get(sample = wc.sample).output.bam
    paraphase_folder = os.path.dirname(paraphase_bam)
    id_genes, = glob_wildcards(f"{paraphase_folder}/{wc.sample}_paraphase_vcfs/{wc.sample}_{{gene}}.vcf")
    vcfs = expand("Sample_{sample}/paraphase/{sample}_paraphase_vcfs/{sample}_{gene}_var_annotated.vcf", sample = wc.sample, gene = id_genes)
    return(vcfs)

rule collect_paraphase_vcfs:
    input:
        find_paraphase_vcf
    output:
        "Sample_{sample}/paraphase/annotation.done"
    shell:
        """
        touch {output}
        """