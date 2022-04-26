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

#____ cDNA READ ORIENTATION AND PRIMER TRIMMING _____________________________#

rule pychopper:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        fq = "Sample_{sample}/{sample}.full_length.fastq",
        unclass = "Sample_{sample}/{sample}.unclassified.fastq",
        rescued = "Sample_{sample}/{sample}.rescued.fastq",
        report = "qc/pychopper/{sample}_report.pdf",
        stats = "qc/pychopper/{sample}_stats.txt"
    log:
        "logs/{sample}_pychopper.log"
    conda:
        "../env/pychopper.yml"
    threads:
        24
    shell:
        """
        fq=$(mktemp)
        unpigz -d -c {input} > $fq 2> {log}
        cdna_classifier.py \
            -r {output.report} \
            -t {threads} \
            -S {output.stats} \
            -u {output.unclass} \
            -w {output.rescued} \
            $fq \
            {output.fq} \
            >> {log} 2>&1
        rm $fq
        """
