rule process_trgt_loci:
    input:
        lambda wc: workflow.source_path("../../resources/{repeat_target}.bed".format(repeat_target = wc.repeat_target))
    output:
        "{repeat_target}.nanorepeat.bed"
    log:
        "logs/{repeat_target}_process_trgt_bed.log"
    shell:
        """
        cat {input} \
            | awk -F '\\t' 'OFS = "\\t" {{split($4, a, ";"); split(a[2], b, "="); match(b[2], /([ACTGN]*)/, c);  print $1, $2, $3, c[1]}}' \
            > {output} 2>{log}        
        """

rule nanorepeat:
    input:
        bam = use_bam,
        ref = config['ref']['genome'],
        bed = "{repeat_target}.nanorepeat.bed".format(repeat_target = config['nanorepeat']['repeat_loci'])
    output:
        directory("Sample_{sample}/repeat_expansions/nanorepeat")
    conda:
        "../env/nanorepeat.yml"
    log:
        "logs/{sample}_nanorepeat.log"
    params:
        data_type = config['nanorepeat']['datatype']
    threads:
        4
    shell:
        """
        nanoRepeat.py \
            -i {input.bam} \
            -t bam \
            -d {params.data_type} \
            -r {input.ref} \
            -b {input.bed} \
            -c {threads} \
            -o {output}/{wildcards.sample} >{log} 2>&1
        """

rule straglr:
    input:
        bam = use_bam,
        ref = config['ref']['genome'],
        bed = workflow.source_path("../../resources/{repeat_target}.nanorepeat.bed".format(repeat_target = config['nanorepeat']['repeat_loci']))
    output:
        "Sample_{sample}/repeat_expansions/straglr.bed"
    conda:
        "../env/straglr.yml"
    log:
        "logs/{sample}_straglr.log"
    params:
        data_type = config['nanorepeat']['datatype']
    threads:
        1
    shell:
        """
        straglr.py \
            {input.bam} \
            {input.ref} \
            Sample_{wildcards.sample}/repeat_expansions/straglr \
            --loci {input.bed} \
            --nprocs {threads} \
        >{log} 2>&1
        """    