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
    shell:
        """
        zcat -f {input.fastqs} | pigz -p {threads} -c > {output}
        echo {input.fastqs} > {log}
        """


rule bam_to_fastq:
    input:
        get_bam_to_fastq_input
    output:
        "Sample_{sample}/{sample}.fastq.gz"
    log:
        "logs/{sample}_bam_to_fastq.log"
    threads:
        8
    conda:
        "../env/samtools.yml"
    shell:
        """
        # Unmapped BAM/CRAM files contain single-end reads; -0 writes all reads to one output file
        samtools cat --threads {threads} {input} \
            | samtools fastq -@ {threads} -0 {output} - 2> {log}
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
