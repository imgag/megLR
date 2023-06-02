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
        1
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
        bonito basecaller \
            {params.model} \
            {input} \
            --modified-bases {params.modbases} \
            --recursive \
            --reference {params.ref} \
            --alignment-threads 16 \
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

#_____ BASECALLING METHYLATION (DORADO) __________________________________________#


rule fast5_to_pod5:
    input:
        unpack(get_input_folders_fast5)
    output:
        pod5="pod5/{sample}/{sample}.pod5",
    log:
        "logs/{sample}_pod5.log"
    params:
        findpath="fast5_*" if config["use_failed_reads"] else "fast5_pass"
    conda:
        "../env/pod5.yml"
    threads:
        8
    shell:
        """
        pod5 convert fast5 \
            --output {output.pod5} \
            --threads {threads} \
            $(find {input} -name '{params.findpath}' -type d -printf '%p ') \
            >{log} 2>&1
        """

rule dorado_basecalling_mod:
    input:
        pod5="pod5/{sample}/{sample}.pod5",
    output:
        bam="Sample_{sample}/{sample}.mod.unmapped.bam",
    log:
        "logs/{sample}_dorado.log"
    threads:
        2
    resources:
        queue="gpu_srv019"
    params:
        dorado = config['apps']['dorado'],
        model = config['dorado']['model'],
        min_qscore = config['dorado']['min_qscore'],
        modbases = f"--modified-bases {config['dorado']['modified_bases']}" if config["dorado"]["modified_bases"] else ""
    shell:
        """
        {params.dorado} basecaller \
            --device "cuda:all" \
            --min-qscore {params.min_qscore} \
            --verbose \
            {params.model} \
            $(dirname {input.pod5}) \
            {params.modbases} 2>{log} | \
            samtools view -b -o {output.bam} \
            >>{log} 2>&1
        """

rule mod_unmapped_fastq:
    input:
        bam="Sample_{sample}/{sample}.mod.unmapped.bam",
    output:
        fastq="Sample_{sample}/{sample}.fastq.gz"
    threads:
        1
    conda:
        "../env/samtools.yml"
    shell:
        """
        samtools fastq -0 {output.fastq} {input.bam}
        """

rule dorado_basecalling_duplex:
    input:
        pod5="pod5/{sample}/{sample}.pod5",
    output:
        bam="Sample_{sample}/{sample}.duplex.unmapped.bam",
    log:
        "logs/{sample}_dorado-duplex.log"
    threads:
        2
    resources:
        queue="gpu_srv019"
    params:
        dorado = config['apps']['dorado'],
        model = config['dorado']['model'],
        min_qscore = config['dorado']['min_qscore'],
    #     modbases = f"--modified-bases {config['dorado']['modified_bases']}" if config["dorado"]["modified_bases"] else ""
    shell:
        """
        {params.dorado} duplex \
            --device "cuda:all" \
            --min-qscore {params.min_qscore} \
            --verbose \
            {params.model} \
            $(dirname {input.pod5}) \
            2>{log} \
            | \
            samtools view -b -o {output.bam} \
            >>{log} 2>&1
        """

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
    threads:
        4
    shell:
        """
        samtools sort -@ 4 -m 4G {input} > {output}
        samtools index {output}
        """
