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

#______ BASECALLING METHYLATION GUPPY  _________________________________________________________#

rule create_fast5_filelist:
    input:
        unpack(get_input_folders_fast5)
    output:
        "basecalling_guppy/{sample}.fast5files.txt"
    shell:
        "find {input} -name '*.fast5' > {output}"

rule guppy_basecalling_mod:
    input:
        "basecalling_guppy/{sample}.fast5files.txt"
    output:
        directory("basecalling_guppy/guppy_{sample}/")
    log:
        "logs/{sample}_guppy_methylation.log"
    threads:
        8
    params:
        guppy = config['apps']['guppy'],
        model = config['guppy']['model'],
        ref = "--align_ref {}".format(config['ref']['genome']) if config['map_during_basecalling'] else "",
        min_qscore = config['guppy']['min_qscore']
    shell:
        """
        GPU_OCCUPIED=$(nvidia-smi --query-compute-apps=gpu_uuid --format=csv,noheader | head -n1)
        if [ -z $GPU_OCCUPIED ] 
        then
            if test -f "/var/lock/gpu0"
            then
                GPU_ACTIVE="1"
            else
                GPU_ACTIVE="0"
                touch /var/lock/gpu0
            fi
        else 
            GPU_ACTIVE=$(nvidia-smi --query-gpu=index,gpu_uuid --format=csv,noheader \
                | grep -v $GPU_OCCUPIED \
                | cut -f1 -d\,)
        fi
        {params.guppy}  \
            --input_file_list {input}\
            --config {params.model} \
            --align_ref {params.ref}\
            --bam_out \
            --index \
            --min_qscore {params.min_qscore} \
            --save_path {output} \
            --device "cuda:$GPU_ACTIVE" \
            >{log} 2>&1
        
        if [$GPU_ACTIVE=="0"]
        then 
            rm /var/lock/gpu0
        fi
        """
    
rule guppy_merge_bams:
    input:
        "basecalling_guppy/guppy_{sample}"
    output:
        "basecalling_guppy/{sample}.mod.guppy.bam"
    conda:
        "../env/samtools.yml"
    shell:
        """
        find {input} -name '.bam' | exec - samtools merge {output} {{}}+
        """

#_____ BASECALLING METHYLATION MAPPING (DORADO) __________________________________________#




#______ PROCESS BASECALLED MOD _____________________________________________________#

rule process_mod_bam:
    input:
        lambda wc: "basecalling_{t}/{s}.mod.{t}.bam".format(s = wc.sample, t = config['basecaller'] )
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