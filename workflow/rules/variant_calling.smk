#____ VARIANT CALLING WITH DEEPVARIANT _______________________________________________________#

#  Inference is speed up massively using GPU support.
rule pepper_marging_deepvariant:
    input:
        bam = rules.map_genome_all.output.bam,
        ref = config['ref']['genome']
    output:
        vcf = "variant_calling/{sample}_pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz"
    log:
        "logs/{sample}_deepvariant_pepper.log"
    params:
        model = "--" + config['vc_pepper']['model'],
        gpu_id = config['gpu_id'],
        target_region = "--region "+ config['vc_pepper']['target_region'] if config['vc_pepper']['target_region'] else "",
        output_phased = "--phased_output" if config['vc']['phased_output'] else "",
        keep_supp = '--pepper_include_supplementary' if config['vc']['include_supplementary'] else ""
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
                {params.target_region} {params.keep_supp} {params.output_phased} {params.model} \
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
                {params.target_region} {params.keep_supp} {params.output_phased} {params.model} \
                >{log} 2>&1
                """
            )

rule copy_vcf:
    input:
        vcf = "variant_calling/{sample}_pepper/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz"
    output:
        vcf = "Sample_{sample}/{sample}.pepper_margin_dv.vcf.gz"
    run:
        if config['output_phased']:
            shell(
                """
                cp {input.vcf} {output.vcf}
                cp {input.vcf}.tbi {output.vcf}.tbi
                cp variant_calling/{sample}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.phased.vcf.gz Sample_{sample}/{sample}.pepper_margin_dv.phased.vcf.gz
                cp variant_calling/{sample}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.phased.vcf.gz.tbi Sample_{sample}/{sample}.pepper_margin_dv.phased.vcf.gz.tbi
                """)
        else: 
            shell(
                """
                cp {input.vcf} {output.vcf}
                cp {input.vcf} {output.vcf}
                """)

#____ VARIANT CALLING WITH CLAIR3 ______________________________________________________#

rule clair3_variants:
    input:
        bam = rules.map_genome_all.output.bam,
        ref = config['ref']['genome']
    output:
        directory("variant_calling/{sample}_clair3")
    conda:
        "../env/clair3.yml" 
    params:
        platform = "ont",
        model = config['vc_clair3']['model'],
        output_phased = "--enable_phasing" if config['vc']['phased_output'] else ""
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
            --output={output} \
            --sample_name={wildcards.sample}  {params.output_phased} \
            >{log} 2>&1        
        """


#____ VARIANT CALLING WITH MEDAKA (DEPRECATED BY ONT) ________________________________________________#

rule medaka_variants:
    input:
        bam = rules.map_genome_all.output.bam,
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
