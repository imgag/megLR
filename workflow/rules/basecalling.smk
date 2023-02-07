#_____ BASECALLING METHYLATION MAPPING (BONITO) __________________________________________#

rule bonito:
    input:
        unpack(get_input_folders_fast5)
    output:
        "basecalling_bonito/{sample}.mod.bonito.bam"
    conda:
        "../env/bonito.yml"
    log:
        "logs/{sample}_bonito.log"
    params:
        model = config['bonito']['model'],
        modbases = config['bonito']['modified_bases'],
        ref = config['ref']['genome']
    threads:
        8
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
        export CUDA_LAUNCH_BLOCKING=1
        
        bonito basecaller \
            {params.model} \
            {input} \
            --modified-bases {params.modbases} \
            --recursive \
            --reference {params.ref} \
            --alignment-threads {threads} \
            --device "cuda:$GPU_ACTIVE" \
            > {output} 2> {log}
         """

rule process_bonito_bam:
    input:
        "basecalling_bonito/{sample}.mod.bonito.bam"
    output:
        "Sample_{sample}/{sample}.mod.bam"
    conda:
        "../env/samtools.yml"
    log:
        "logs/{sample}_bonito_process.log"
    threads: 4
    shell:
        """
        samtools sort -@ 4 -m 4G {input} > {output}
        samtools index {output}
        """

#______ BASECALLING METHYLATION GUPPY  _________________________________________________________#

rule create_fast5_filelist:
    input:
        unpack(get_input_folders_fast5)
    output:
        "basecalling_guppy/{wildcards.sample}.fast5files.txt"
    shell:
        "find {input} -name '*.fast5' > {output}"

rule guppy_basecalling_with_mapping:
    input:
        "basecalling_guppy/{wildcards.sample}.fast5files.txt"
    output:
        directory("basecalling_guppy/{sample}")
    log:
        "logs/{sample}_guppy_methylation.log"
    threads:
        8
    params:
        guppy = config['apps']['guppy'],
        model = config['guppy']['model'],
        ref = config['ref']['genome'],
        min_q_score = config['guppy']['min_qscore'],
        gpu = config['gpu_id']['cuda']
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
        export CUDA_LAUNCH_BLOCKING=1
        {params.guppy}  \
            --input_file_list {input}\
            --config {params.model} \
            --align_ref 
            --bam_out \
            --index \
            --minimap_opt_string --eqx \
            --num_alignment_threads 4 \
            --min_qscore {params.qscore} \
            --save_path {output} \
            --device "cuda:$GPU_ACTIVE" \
        """

rule guppy_basecalling:
    input:
        "basecalling_guppy/{wildcards.sample}.fast5files.txt"
    output:
        directory("basecalling_guppy/{sample}")
    log:
        "logs/{sample}_guppy_methylation.log"
    threads:
        8
    params:
        guppy = config['apps']['guppy'],
        model = config['guppy']['model'],
        ref = config['ref']['genome'],
        gpu = config['gpu_id']['cuda']
    shell:
        """
        {params.guppy}  \
            --config {params.model} \
            --save_path {output} \
            --recursive \
            --device {params.gpu} \
        """
