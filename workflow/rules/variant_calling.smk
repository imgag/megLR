#____ VARIANT CALLING WITH DEEPVARIANT _______________________________________________________#

#  Inference is sped up massively using GPU support.
rule pepper_marging_deepvariant:
    input:
        bam = use_bam,
        ref = config['ref']['genome']
    output:
        vcf = "variant_calling/{sample}_pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz"
    log:
        "logs/{sample}_deepvariant_pepper.log"
    params:
        model = "--" + config['vc_pepper']['model'],
        gpu_id = config['gpu_id']['id'],
        target_region = "--region "+ config['vc_pepper']['target_region'] if config['vc_pepper']['target_region'] else "",
        phased_output = "--phased_output" if config['vc']['phased_output'] else "",
        keep_supp = '--pepper_include_supplementary' if config['vc']['keep_supplementary'] else ""
    threads: 20
    run:
        if config['use_gpu']:
            shell(
                """
                docker run \
                -v "$(dirname $(realpath {input.bam}))":"/mnt/input_bam" \
                -v "$(dirname $(realpath {input.ref}))":"/mnt/input_ref" \
                -v "$(dirname $(realpath {output.vcf}))":"/mnt/output" \
                --user $(id -u):$(id -g) \
                --gpus device={params.gpu_id} \
                kishwars/pepper_deepvariant:r0.8-gpu \
                run_pepper_margin_deepvariant call_variant \
                --bam "/mnt/input_bam/$(basename {input.bam})" \
                --fasta "/mnt/input_ref/$(basename {input.ref})" \
                --threads 8 \
                --gpu \
                --output_dir "/mnt/output" \
                {params.target_region} {params.keep_supp} {params.phased_output} {params.model} \
                >{log} 2>&1
                """
            )
        else:
            shell(
                """
                docker run \
                -v "$(dirname $(realpath {input.bam}))":"/mnt/input_bam" \
                -v "$(dirname $(realpath {input.ref}))":"/mnt/input_ref" \
                -v "$(dirname $(realpath {output.vcf}))":"/mnt/output" \
                --user $(id -u):$(id -g) \
                kishwars/pepper_deepvariant:r0.8 \
                run_pepper_margin_deepvariant call_variant \
                --bam "/mnt/input_bam/$(basename {input.bam})" \
                --fasta "/mnt/input_ref/$(basename {input.ref})" \
                --threads {threads} \
                --output_dir "/mnt/output" \
                {params.target_region} {params.keep_supp} {params.phased_output} {params.model} \
                >{log} 2>&1
                """
            )

rule copy_vcf_pepper:
    input:
        vcf = "variant_calling/{sample}_pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz"
    output:
        vcf = "Sample_{sample}/{sample}.pepper_margin_dv.vcf.gz"
    run:
        if config['vc']['phased_output']:
            shell(
                """
                cp {input.vcf} {output.vcf}
                cp {input.vcf}.tbi {output.vcf}.tbi
                cp variant_calling/{wildcards.sample}_pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.phased.vcf.gz Sample_{wildcards.sample}/{wildcards.sample}.pepper_margin_dv.phased.vcf.gz
                cp variant_calling/{wildcards.sample}_pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.phased.vcf.gz.tbi Sample_{wildcards.sample}/{wildcards.sample}.pepper_margin_dv.phased.vcf.gz.tbi
                """)
        else: 
            shell(
                """
                cp {input.vcf} {output.vcf}
                cp {input.vcf} {output.vcf}
                """)

rule output_haplotagged_bam:
    input:
        vcf = "variant_calling/{sample}_pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz"
    output:
        bam = "Sample_{sample}/{sample}.haplotagged.bam"
    threads:
        4
    shell:
        """
        mv variant_calling/{wildcards.sample}_pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.haplotagged.bam {output}
        samtools index -@{threads} {output}
        """

#____ VARIANT CALLING WITH CLAIR3 ______________________________________________________#

rule clair3_variants:
    input:
        bam = use_bam,
        ref = config['ref']['genome']
    output:
        "variant_calling/{sample}_clair3/merge_output.vcf.gz"
    conda:
        "../env/clair3.yml" 
    params:
        platform = "ont",
        model = config['vc_clair3']['model'],
        phased_output = "--enable_phasing" if config['vc']['phased_output'] else ""
    threads:
        20
    log:
        "logs/{sample}_clair3.log"
    shell:
        """
        run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.ref} \
            --threads={threads} \
            --platform={params.platform} \
            --model_path="${{CONDA_PREFIX}}/bin/models/{params.model}" \
            --output=$(dirname {output}) \
            --sample_name={wildcards.sample}  {params.phased_output} \
            >{log} 2>&1        
        """

rule copy_vcf_clair:
    input:
        vcf = "variant_calling/{sample}_clair3/merge_output.vcf.gz"
    output:
        vcf = "Sample_{sample}/{sample}.clair3.vcf.gz"
    log:
        "logs/{sample}_copy_vcf_clair3.log"
    run:
        if config['vc']['phased_output']:
            shell(
                """
                tabix 
                cp {input.vcf} {output.vcf} >{log} 2>&1
                cp {input.vcf}.tbi {output.vcf}.tbi >{log} 2>&1
                cp variant_calling/{wildcards.sample}_clair3/phased_merge_output.vcf.gz Sample_{wildcards.sample}/{wildcards.sample}.clair3.phased.vcf.gz > {log} 2>&1
                cp variant_calling/{wildcards.sample}_clair3/phased_merge_output.vcf.gz.tbi Sample_{wildcards.sample}/{wildcards.sample}.clair3.phased.vcf.gz.tbi > {log} 2>&1
                """)
        else: 
            shell(
                """
                cp {input.vcf} {output.vcf} >{log} 2>&1
                cp {input.vcf} {output.vcf} >{log} 2>&1
                """)

#____ VARIANT CALLING WITH MEDAKA (DEPRECATED BY ONT) ________________________________________________#

rule medaka_variants:
    input:
        bam = use_bam,
        ref = config['ref']['genome']
    output:
        vcf = "variant_calling/{sample}_medaka/{chr}/round_1.vcf"
    conda:
        "../env/medaka.yml"
    threads:
        3
    log:
        "logs/{sample}_{chr}_medaka.log"
    params:
        model_snp = config['vc_medaka']['model_initial'],
        model_final = config['vc_medaka']['model_final']
    shell:
        """
        export OMP_NUM_THREADS={threads}
        medaka_variant -d \
            -f {input.ref} \
            -i {input.bam} \
            -r {wildcards.chr} \
            -p  \
            -o variant_calling/{wildcards.sample}_medaka/{wildcards.chr} \
            -t {threads} \
            -s {params.model_snp} \
            -m {params.model_final} \
            >{log} 2>&1
        """

rule combine_vcf:
    input:
        expand("variant_calling/{{sample}}_medaka/{C}/round_1.vcf",  C = get_chromosomes())
    output:
        "variant_calling/{sample}_medaka/{sample}_var.vcf"
    conda:
        "../env/bcftools.yml"
    threads:
        1
    log:
        "logs/{sample}_combinevcf.log"
    group:
        "medaka"
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
        "Sample_{sample}/{sample}.medaka.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/{sample}_processvcf.log"
    group:
        "medaka"
    shell:
        """
        bgzip -c {input} > {output} 2> {log}
        tabix {output} 2>>{log}
        """
