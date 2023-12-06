#____ VARIANT CALLING WITH DEEPVARIANT _______________________________________________________#
# Works only for R10 Flowcells !
# Todo implement downstream phasing using longphase

if config['use_gpu']:
    ruleorder:
        deepvariant_gpu > deepvariant
else:
    ruleorder:
        deepvariant > deepvariant_gpu

rule deepvariant:
    input:
        bam = use_bam,
        ref = config['ref']['genome']
    output:
        vcf = "variant_calling/{sample}_deepvariant/{sample}.dv.vcf.gz",
        gvcf =  "variant_calling/{sample}_deepvariant/{sample}.dv.gvcf.gz"
    log:
        "logs/{sample}_deepvariant.log"
    params:
        version = "1.6.0"
    threads: 20
    shell:
        """
        docker run \
        -v "$(dirname $(realpath {input.bam}))":"/mnt/input_bam" \
        -v "$(dirname $(realpath {input.ref}))":"/mnt/input_ref" \
        -v "$(dirname $(realpath {output.vcf}))":"/mnt/output" \
        -v "/tmp":"/tmp" \
        --user $(id -u):$(id -g) \
        google/deepvariant:{params.version} \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type="ONT_R104" \
        --ref="/mnt/input_ref/$(basename {input.ref})" \
        --reads="/mnt/input_bam/$(basename {input.bam})" \
        --output_vcf="/mnt/output/$(basename {output.vcf})" \
        --output_gvcf="/mnt/output/$(basename {output.gvcf})" \
        --num_shards={threads} \
        >{log} 2>&1
        """

rule deepvariant_gpu:
    input:
        bam = use_bam,
        ref = config['ref']['genome']
    output:
        vcf = "variant_calling/{sample}_deepvariant/{sample}.dv.vcf.gz",
        gvcf =  "variant_calling/{sample}_deepvariant/{sample}.dv.gvcf.gz"
    log:
        "logs/{sample}_deepvariant.log"
    params:
        version = "1.6.0",
        gpu_id = config['gpu_id']['id'],
    threads: 1
    resources:
        queue="gpu_srv019"
    shell:
        """
        GPU_OCCUPIED=$(nvidia-smi --query-compute-apps=gpu_uuid --format=csv,noheader | head -n1)
        if [ -z $GPU_OCCUPIED ] 
        then
            GPU_ACTIVE="0"
        else 
            GPU_ACTIVE=$(nvidia-smi --query-gpu=index,gpu_uuid --format=csv,noheader \
                | grep -v $GPU_OCCUPIED \
                | cut -f1 -d\,)
        fi
        docker run \
        -v "$(dirname $(realpath {input.bam}))":"/mnt/input_bam" \
        -v "$(dirname $(realpath {input.ref}))":"/mnt/input_ref" \
        -v "$(dirname $(realpath {output.vcf}))":"/mnt/output" \
        -v "/tmp":"/tmp" \
        -e CUDA_LAUNCH_BLOCKING=1 \
        --user $(id -u):$(id -g) \
        --gpus device="$GPU_ACTIVE" \
        google/deepvariant:{params.version}-gpu \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type="ONT_R104" \
        --ref="/mnt/input_ref/$(basename {input.ref})" \
        --reads="/mnt/input_bam/$(basename {input.bam})" \
        --output_vcf="/mnt/output/$(basename {output.vcf})" \
        --output_gvcf="/mnt/output/$(basename {output.gvcf})" \
        --num_shards=20 \
        >{log} 2>&1
        """

rule copy_vcf_deepvariant:
    input:
        vcf = "variant_calling/{sample}_deepvariant/{sample}.dv.vcf.gz"
    output:
        vcf = "Sample_{sample}/{sample}.dv.vcf.gz"
    log:
        "logs/{sample}_copy_vcf_deepvariant.log"
    shell:
        """
        cp {input.vcf} {output.vcf} >{log} 2>&1
        """
#____ VARIANT CALLING WITH PEPPER_MARGIN_DEEPVARIANT _________________________________________#

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
                GPU_OCCUPIED=$(nvidia-smi --query-compute-apps=gpu_uuid --format=csv,noheader | head -n1)
                if [ -z $GPU_OCCUPIED ] 
                then
                    GPU_ACTIVE="0"
                else 
                    GPU_ACTIVE=$(nvidia-smi --query-gpu=index,gpu_uuid --format=csv,noheader \
                        | grep -v $GPU_OCCUPIED \
                        | cut -f1 -d\,)
                fi
                docker run \
                -v "$(dirname $(realpath {input.bam}))":"/mnt/input_bam" \
                -v "$(dirname $(realpath {input.ref}))":"/mnt/input_ref" \
                -v "$(dirname $(realpath {output.vcf}))":"/mnt/output" \
                -v "/tmp":"/tmp" \
                -e CUDA_LAUNCH_BLOCKING=1
                --user $(id -u):$(id -g) \
                --gpus device="cuda:$GPU_ACTIVE" \
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
                mkdir -p $(dirname $(realpath {output.vcf}))
                singularity run \
                -B "$(dirname $(realpath {input.bam}))":"/mnt/input_bam" \
                -B "$(dirname $(realpath {input.ref}))":"/mnt/input_ref" \
                -B "$(dirname $(realpath {output.vcf}))":"/mnt/output" \
                -B "/tmp":"/tmp" \
                docker://kishwars/pepper_deepvariant:r0.8 \
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

rule pepper_haplotagged_bam:
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
            --model_path={params.model} \
            --output=$(dirname {output}) \
            --sample_name={wildcards.sample} \
            >{log} 2>&1        
        """

rule copy_vcf_clair:
    input:
        vcf = "variant_calling/{sample}_clair3/merge_output.vcf.gz"
    output:
        vcf = "Sample_{sample}/{sample}.clair3.vcf.gz"
    log:
        "logs/{sample}_copy_vcf_clair3.log"
    shell:
        """
        cp {input.vcf} {output.vcf} >{log} 2>&1
        """

#____ PHASING USING LONGPHASE ______________________________________________________#

rule longphase_phase:
    input:
        vcf_snp = "Sample_{sample}/{sample}.{vc}.vcf.gz",
        vcf_sv = "Sample_{{sample}}/{{sample}}.{sv}.vcf.gz".format(sv = config['phasing']['sv_vcf']),
        bam = "Sample_{sample}/{sample}.bam",
        ref = config['ref']['genome']
    output:
        vcf = "Sample_{sample}/{sample}.phased.{vc}.vcf.gz",
        vcf_sv = "Sample_{sample}/{sample}.phased.{vc}_SV.vcf.gz"
    threads:
        4
    log:
        "logs/{sample}_longphase_phase_{vc}.log"
    params:
        longphase = config['apps']['longphase']
    shell:
        """
        {params.longphase} phase \
        --snp-file {input.vcf_snp} \
        --sv-file {input.vcf_sv} \
        --bam-file {input.bam} \
        --ont \
        --reference {input.ref} \
        --threads {threads} \
        --out-prefix Sample_{wildcards.sample}/{wildcards.sample}.phased.{wildcards.vc} \
        > {log} 2>&1

        bgzip Sample_{wildcards.sample}/{wildcards.sample}.phased.{wildcards.vc}.vcf
        bgzip Sample_{wildcards.sample}/{wildcards.sample}.phased.{wildcards.vc}_SV.vcf
        tabix {output.vcf}
        tabix {output.vcf_sv}
        """

rule longphase_haplotag:
    input:
        vcf_snp = "Sample_{{sample}}/{{sample}}.phased.{vc}.vcf.gz".format(vc = config['phasing']['haplotagging_input']),
        bam = "Sample_{sample}/{sample}.bam",
    output:
        "Sample_{sample}/{sample}.haplotagged.bam"
    threads:
        4
    log:
        "logs/{sample}_longphase_haplotagging.log"
    params:
        longphase = config['apps']['longphase']
    shell:
        """
        {params.longphase} haplotag \
        --snp-file {input.vcf_snp} \
        --bam-file {input.bam} \
        --threads {threads} \
        --out-prefix Sample_{wildcards.sample}/{wildcards.sample}.haplotagged \
        > {log} 2>&1

        samtools index {output}
        """

        