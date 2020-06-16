#____ ASSEMBLY _______________________________________________________________#

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
        cp $outf/wtdbg.ctg.lay.gz {output}
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
        12
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        wtpoa-cns \
            -t 16 \
            -i {input} \
            -fo {output}
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
        16
    shell:
        """
        outf=assembly/{wildcards.sample}_wtdbg
        minimap2 -t{threads} -ax map-ont -r2k {input.asm} {input.reads} \
            | samtools sort -@4 > $outf/wtdbg.polishing.bam
        samtools view -F0x900 $outf/wtdbg.polishing.bam \
            | wtpoa-cns -t {threads} -d {input.asm} -i - -fo $outf/wtdbg.pol.fa
        cp $outf/wtdbg.pol.fa {output}
        rm $outf/wtdbg.polishing.bam
        """

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

#____ ASSEMBLY QC _____________________________________________________________#

rule quast:
    input:
        files = expand("Sample_{s}/{s}.asm.{asm}.fasta",
            s=ID_samples,
            asm=config['assembly']['methods'])
    output:
        "qc/quast_results/report.tsv"
    conda:
        "../env/quast.yml"
    log:
        "logs/quast.log"
    params:
        ref=config['ref']['genome']
    threads:
        8
    shell:
        """
        quast \
            --threads {threads} \
            --no-sv             \
            --reference {params.ref} \
            --output-dir qc/quast_results 
            {input}
        """
        
