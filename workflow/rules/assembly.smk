#____ ASSEMBLY USING WTDBG2 ________________________________________________________#

rule wtdbg2:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        "assembly/{sample}_wtdbg/wtdbg.ctg.lay.gz"
    conda:
        "../env/wtdbg2.yml"
    log:
        "logs/{sample}_wtdbg.log"
    params:
        g=config['assembly']['genome_size']
    threads:
        config['sys']['max_threads']
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        mkdir -p $outf
        wtdbg2 \
            -x nanopore     \
            -g {params.g}   \
            -t {threads}    \
            -i {input}      \
            -fo $outf/wtdbg \
            > {log} 2>&1
        """

rule wtdbg2_consensus:
    input:
        "assembly/{sample}_wtdbg/wtdbg.ctg.lay.gz"
    output:
        "assembly/{sample}_wtdbg/wtdbg.raw.fa"
    conda:
        "../env/wtdbg2.yml"
    log:
        "logs/{sample}_wtdbg_consensus.log"
    threads:
        config['sys']['max_threads']
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        wtpoa-cns \
            -t {threads}\
            -i {input}  \
            -fo {output}\
            > {log} 2>&1
        """

rule wtdbg2_polishing:
    input:
        asm = "assembly/{sample}_wtdbg/wtdbg.raw.fa",
        reads = "Sample_{sample}/{sample}.fastq.gz"
    output:
        "Sample_{sample}/{sample}.asm.wtdbg.fasta"
    conda:
        "../env/wtdbg2.yml"
    log:
        "logs/{sample}_wtdbg_consensus.log"
    threads:
        config['sys']['max_threads']
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        minimap2 -t{threads} -ax map-ont -r2k {input.asm} {input.reads} \
            | samtools sort -@4 > $outf/wtdbg.polishing.bam
        samtools view -F0x900 $outf/wtdbg.polishing.bam \
            | wtpoa-cns -t {threads} -d {input.asm} -i - -fo $outf/wtdbg.pol.fa \
            > {log} 2>&1
        cp $outf/wtdbg.pol.fa {output}
        rm $outf/wtdbg.polishing.bam
        """

#____ ASSEMBLY USING FLYE ________________________________________________________#

rule flye:
    input:
        "Sample_{sample}/{sample}.fastq.gz"
    output:
        "Sample_{sample}/{sample}.asm.flye.fasta"
    conda:
        "../env/flye.yml"
    log:
        "logs/{sample}_flye.log"
    params:
        g=config['assembly']['genome_size']
    threads:
        config['sys']['max_threads']
    shell:
        """
        flye \
            --nano-raw      {input}     \
            --genome-size   {params.g}  \
            --threads       {threads}   \
            --out-dir       assembly/{wildcards.sample}_flye/ \
            > {log} 2>&1
        cp assembly/{wildcards.sample}_flye/assembly.fasta {output}
        """
        
