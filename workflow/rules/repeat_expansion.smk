rule process_trgt_loci:
    input:
        srcdir("../../resources/{repeat_target}.bed")
    output:
        srcdir("../../resources/{repeat_target}.nanorepeat.bed")
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
        bed = srcdir("../../resources/{repeat_target}.nanorepeat.bed".format(repeat_target = config['nanorepeat']['repeat_loci']))
    output:
        directory("Sample_{sample}/repeat_expansions")
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