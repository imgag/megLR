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
    wildcard_constraints:
        gene="[^.]*"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/{sample}.{gene}.paraphase_add_sample_header.log"
    shell:
        """
        samples=$(bcftools query -l {input.vcf} | sed 's/.*/##SAMPLE=<ID=&>/') 2>{log}
        bcftools annotate -h <(echo "${{samples}}") {input.vcf} > {output.vcf} 2>>{log}
        """

rule paraphase_annotate_vcf:
    input:
        vcf="Sample_{sample}/paraphase/{sample}_paraphase_vcfs/{sample}_{gene}.reheader.vcf"
    output:
        vcf="Sample_{sample}/paraphase/{sample}_paraphase_vcfs/{sample}_{gene}_var_annotated.vcf.gz"
    wildcard_constraints:
        gene="[^.]*"
    log:
        "logs/{sample}.{gene}.paraphase_annotate_vcf.log"
    params:
        megsap = config['apps']['megsap'],
        system = config['paraphase']['annotation_system']
    shell:
        """
        php {params.megsap}/src/Pipelines/annotate.php \
            -vcf {input} \
            -multi \
            -out_folder $(dirname {output.vcf}) \
            -out_name {wildcards.sample}_{wildcards.gene} \
            -system {params.system} > {log} 2>{log}
        """

rule paraphase_vcf_to_tsv:
    input:
        "Sample_{sample}/paraphase/{sample}_paraphase_vcfs/{sample}_{gene}_var_annotated.vcf.gz"
    output:
        "Sample_{sample}/paraphase/{sample}_{gene}_var_annotated.tsv"
    wildcard_constraints:
        gene="[^.]*"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/{sample}.{gene}.paraphase_vcf_to_tsv.log"
    shell:
        """
        gt_samples=$(bcftools query -l {input} | sed -e 's/^/GT_/' | tr '\\n' '\\t') 2>{log}
        dp_samples=$(bcftools query -l {input} | sed -e 's/^/DP_/' | tr '\\n' '\\t') 2>>{log}
        echo "CHR\tPOS\tREF\tALT\t${{gt_samples}}${{dp_samples}}CSQ\tGNOMAD_AF\tPHYLOP\tCADD_SNV" > {output} 2>>{log}
        bcftools query \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT\\t][%DP\\t]%INFO/CSQ/\\t%INFO/gnomADg_AF\\t%INFO/PHYLOP\\t%INFO/CADD_SNV\\n' \
            {input} >> {output} 2>> {log}
        """

def find_paraphase_outs(wc):
    paraphase_bam = checkpoints.paraphase.get(sample = wc.sample).output.bam
    paraphase_folder = os.path.dirname(paraphase_bam)
    id_genes, = glob_wildcards(f"{paraphase_folder}/{wc.sample}_paraphase_vcfs/{wc.sample}_{{gene}}.vcf")
    id_genes=[x for x in id_genes if "reheader" not in x]
    outs = expand("Sample_{sample}/paraphase/{sample}_{gene}_var_annotated.tsv", sample = wc.sample, gene = id_genes)
    return(outs)

rule collect_paraphase_outs:
    input:
        find_paraphase_outs
    output:
        "Sample_{sample}/paraphase/annotation.done"
    shell:
        """
        touch {output}
        """