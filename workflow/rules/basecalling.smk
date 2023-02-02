#_____ BASECALLING METHYLATION MAPPING (BONITO) __________________________________________#

rule bonito:
    input:
        unpack(get_input_folders_fast5)
    output:
        "bonito_basecalled/{sample}.mod.bonito.bam"
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
        export CUDA_LAUNCH_BLOCKING=1
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
        "bonito_basecalled/{sample}.mod.bonito.bam"
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

#______ BASECALLING WITH GUPPY _________________________________________________________#

