#____ EXTRACT READS FROM TARGET REGION ___________________________________________#
rule create_target_beds:
    output:
        "local_assembly/{target}.bed"
    params:
        t = lambda wildcards: wildcards['target']
    run:
        with open(output[0], "w") as out:
            bed_line = params['t'].replace(':', '\t').replace('-','\t').replace('_', '\t')
            bed_line = bed_line.replace('.hp1', '').replace('.hp2', '')
            out.write(bed_line)
    
rule extract_read_ids:
    input:
        bam = "Sample_{sample}/{sample}.bam",
        region = "local_assembly/{target}.bed"
    output:
        "local_assembly/Sample_{sample}/{target,[^.]+}.read_ids.txt" #Regex: Does not contain dot
    conda:
        "../env/samtools.yml"
    shell:
        """
        samtools view {input.bam} --region-file {input.region} | cut -f 1 | sort | uniq > {output}
        """

rule extract_phased_read_ids:
    input:
        bam = "Sample_{sample}/{sample}.haplotagged.bam",
        region = "local_assembly/{target}.bed"
    output:
        "local_assembly/Sample_{sample}/{target}.hp{phase}.read_ids.txt"
    conda:
        "../env/samtools.yml"
    shell:
        """
        samtools view {input.bam} \
            --region-file {input.region} \
            --tag HP:{wildcards.phase} \
            | cut -f 1 \
            | sort \
            | uniq > {output}
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
        fa = "local_assembly/Sample_{sample}/{target}.asm.fasta",
        gfa = "local_assembly/Sample_{sample}/{target}.asm.gfa"
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
    params:
        region = lambda wc: wc.target.replace('_',':')
    shell:
        """    
        samtools faidx {input.ref} {params.region} > {output}
        """

rule syri:
    input:
        sam = "local_assembly/Sample_{sample}/{target}.aln_region.sam",
        asm = "local_assembly/Sample_{sample}/{target}.asm.fasta",
        ref = lambda wc: "local_assembly/"+ wc.target.replace('.hp1',"").replace('.hp2',"") +".ref.fasta"
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
            -reference {input.ref} \
            -t {input.asm} \
            --verbose \
            --format json \
            --num-threads {threads} \
            >{log} 2>&1
        """

rule smashpp_viz:
    input:
        json="Sample_{sample}/local_assembly/{sample}.{target}.json"    
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

rule local_assembly_mapping_genome:
    input:
        asm = "local_assembly/Sample_{sample}/{target}.asm.fasta",
        ref = config['ref']['genome']
    output:
        bam = "local_assembly/Sample_{sample}/{target}.aln_ref.bam"
    conda:
        "../env/minimap2.yml"
    log:
        "logs/{sample}_{target}_minimap2.log"
    shell:
        """
        minimap2 -ax asm5 --eqx {input.ref} {input.asm} \
            | samtools sort -O BAM - > {output.bam}
        samtools index {output}
        """

rule local_assembly_mapping_region:
    input:
        asm = "local_assembly/Sample_{sample}/{target}.asm.fasta",
        ref = lambda wc: "local_assembly/"+ wc.target.replace('.hp1',"").replace('.hp2',"") +".ref.fasta"
    output:
        sam = "local_assembly/Sample_{sample}/{target}.aln_region.sam"
    conda:
        "../env/minimap2.yml"
    log:
        "logs/{sample}_{target}_minimap2.log"
    shell:
        """
        minimap2 -ax asm5 --eqx {input.ref} {input.asm} > {output.sam}
        """