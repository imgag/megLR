#_____ COMPRESS TO SINGLE .FASTQ.GZ __________________________________________#

rule join_fastq:
    input:
        unpack(get_input_folders)
    output:
        "Sample_{sample}/{sample}.fastq.gz"
    log:
        "logs/{sample}_fetch_reads.log"
    threads:
        8
    conda:
        "../env/pigz.yml"
    params:
        exclude_failed = "-not -path '*fail*' -a" if not config['use_failed_reads'] else ""
    shell:
        """
        find {input.folders} -type f  \
        {params.exclude_failed} '(' \
          -name '*.fastq' -o \
          -name '*.fastq.gz' -o \
          -name '*.fq' -o \
          -name '*.fq.gz' \
          ')'  -exec zcat -f {{}} + \
            | pigz -p {threads} -c > {output}
        
        echo $(find {input.folders} -type f  \
        {params.exclude_failed} '(' \
          -name '*.fastq' -o \
          -name '*.fastq.gz' -o \
          -name '*.fq' -o \
          -name '*.fq.gz' \
          ')') > {log}
        """


rule join_bam:
    input:
        unpack(get_input_folders_bam)
    output:
        bam="Sample_{sample}/{sample}.mod.unmapped.bam",
    log:
        "logs/{sample}_join_bam.log"
    conda:
        "../env/samtools.yml"
    threads:
        2
    params:
        exclude_failed = "-not -path '*fail*' -a" if not config['use_failed_reads'] else ""
    shell:
        """
        find {input} -type f {params.exclude_failed} -name '*.bam' | \
        samtools cat \
            --threads {threads} \
            -o {output.bam} \
            -b - \
            >{log} 2>&1
        """
