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
        ref = config['ref']['genome'],
        gpu = config['gpu_id']['cuda']
    threads:
        8
    shell:
        """
        bonito basecaller \
            {params.model} \
            {input} \
            --modified-bases {params.modbases} \
            --recursive \
            --reference {params.ref} \
            --alignment-threads {threads} \
            --device '{params.gpu}' \
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
            --device {params.gpu} \
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
