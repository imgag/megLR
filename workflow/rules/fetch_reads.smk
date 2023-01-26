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
