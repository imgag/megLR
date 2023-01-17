#____ EXTRACT READS FROM TARGET REGION ___________________________________________#
rule create_target_beds:
    output:
        "local_assembly/{target}.bed"
    params:
        t = lambda wildcards: wildcards['target']
    run:
        with open(output[0], "w") as out:
            bed_line = params['t'].replace(':', '\t').replace('-','\t')
            out.write(bed_line)
    
rule extract_read_ids:
    input:
        bam = "Sample_{sample}/{sample}.bam",
        region = "local_assembly/{target}.bed"
    output:
        "local_assembly/Sample_{sample}/{target}.read_ids.txt"
    conda:
        "../env/samtools.yml"
    shell:
        """
        samtools view {input.bam} --region-file {input.region} | cut -f 1 | sort | uniq > {output}
        """

rule extract_reads:
    input:
        fastq = "Sample_{sample}/{sample}.fastq.gz",
        read_ids = "local_assembly/Sample_{sample}/{target}.read_ids.txt"
    output:
        fastq = "local_assembly/Sample_{sample}/{target}.fastq"
    conda:
        "../env/seqtk.yml"
    shell:
        """
        seqtk subseq {input.fastq} {input.read_ids} > {output.fastq}
        """

#____ ASSEMBLY USING RAVEN ________________________________________________________#

rule raven_assmbly:
    input:
        fastq = "local_assembly/Sample_{sample}/{target}.fastq"
    output:
        fa = "local_assembly/Sample_{sample}/{target}.fasta",
        gfa = "local_assembly/Sample_{sample}/{target}.gfa"
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
        ref = config['ref']['genome']
    output:
        "local_assembly/{target}.ref.fasta"
    shell:
        """
        samtools faidx {input.ref} {wildcards.target} > {output}
        """

rule smashpp:
    input:
        asm = "local_assembly/Sample_{sample}/{target}.fasta",
        ref = "local_assembly/{target}.ref.fasta"
    output:
        json = "Sample_{sample}/local_assembly/{sample}.{target}.json"    
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
        "local_assembly/Sample_{sample}/{target}.synteny.svg"
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


#____ MINIMAP2 PAIRWISE GENOME ALN ________________________________________________________#