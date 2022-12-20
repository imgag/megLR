#____ EXTRACT READS FROM TARGET REGION ___________________________________________#
rule create_target_beds:
    output:
        "Sample_{sample}/local_assembly/{target}.bed"
    params:
        t = lambda wildcards: wildcards['target']
    run:
        with open(output[0], "w") as out:
            bed_line = params['t'].replace(':', '\t').replace('-','\t')
            out.write(bed_line)
    
rule extract_read_ids:
    input:
        bam = "Sample_{sample}/{sample}.bam",
        region = "Sample_{sample}/local_assembly/{target}.bed"
    output:
        "Sample_{sample}/local_assembly/{target}.read_ids.txt"
    conda:
        "../env/samtools.yml"
    shell:
        """
        samtools view {input.bam} --region-file {input.region} | cut -f 1 | sort | uniq > {output}
        """

rule extract_reads:
    input:
        fastq = "Sample_{sample}/{sample}.fastq.gz",
        read_ids = "Sample_{sample}/local_assembly/{target}.read_ids.txt"
    output:
        fastq = "Sample_{sample}/local_assembly/{target}.fastq"
    conda:
        "../env/seqtk.yml"
    shell:
        """
        seqtk subseq {input.fastq} {input.read_ids} > {output.fastq}
        """

#____ ASSEMBLY USING RAVEN ________________________________________________________#

rule raven_assmbly:
    input:
        fastq = "Sample_{sample}/local_assembly/{target}.fastq"
    output:
        fa = "Sample_{sample}/local_assembly/{target}.fasta",
        gfa = "Sample_{sample}/local_assembly/{target}.gfa"
    conda:
        "../env/raven.yml"
    log:
        "logs/{sample}_{target}_raven.log"
    params:
        polishing = 2
    threads:
        8
    shell:
        """
        raven \
            --threads {threads} \
            --graphical-fragment-assembly {output.gfa} \
            --polishing-rounds {params.polishing} \
            {input.fastq} > {output.fa} 2>{log}
        """

#____ SMASHPP PAIRWISE GENOME COMP ________________________________________________________#
rule extract_ref_seq:
    input:
        ref = config['ref']
    output:
        "Sample_{sample}/local_assembly/{target}.ref.fasta"
    shell:
        """
        samtools faidx {input.ref} {wildcards.target} > {output}
        """

rule smashpp:
    input:
        asm = "Sample_{sample}/local_assembly/{target}.fasta",
        ref = config['ref']
    output:
        json = "{sample}.{target}.json"    
    conda:
        "../env/smashpp.yml"
    log:
        "logs/{sample}_{target}_smashpp.log"
    threads:
        4
    shell:
        """
        smashpp \
            --reference {input.ref} \
            --target {input.asm} \
            --verbose \
            --format json \
            --num-threads {threads} \
            >{log} 2>&1
        """

rule smashpp_viz:
    input:
        json = rules.smashpp.output.json,
    output:
        "Sample_{sample}/local_assembly/{target}.synteny.svg"
    conda:
        "../env/smashpp.yml"
    log:
        "logs/{sample}_{target}_smashpp_viz.log"
    shell:
        """
        smashpp viz\
            --output {output} \
            {input.json} \
            >{log} 2>&1
        """
