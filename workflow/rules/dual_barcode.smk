#____ cDNA READ ORIENTATION AND PRIMER TRIMMING _____________________________#

rule demux_minibar:
    input:
        fq = "Sample_{sample}/{sample}.fastq.gz",
        bc = config['demux']['bc_index']
    output:
        folder = directory("Sample_{sample}/demux")
    threads:
        2
    log:
        "logs/{sample}_minibar.log"
    params:
        minibar = config['demux']['minibar'],
        dist = config['demux']['bc_distance']
    shell:
        """
        mkdir -p {output}
        {params.minibar} \
            -p {params.dist} \
            -T -F -P Sample_{wildcards.sample}/demux/{wildcards.sample}_ \
            {input.bc} {input.fq} > {log} 2>&1
        """

rule demux_minibar_stats:
    input:
        fq = "Sample_{sample}/{sample}.fastq.gz",
        bc = config['demux']['bc_index']
    output:
        "Sample_{sample}/{sample}.demux.txt"
    log:
        "logs/{sample}_minibar_stats.log"
    threads:
        2
    params:
        minibar = config['demux']['minibar'],
        dist = config['demux']['bc_distance']
    shell:
        """
        {params.minibar} \
            -D -p {params.dist} \
            {input.bc} {input.fq} > {output} 2> {log}
        """

rule demux_qc:
    input:
        demux_stats = rules.demux_minibar_stats.output
    output:
        pdf = "Sample_{sample}/{sample}_barcode_distribution.pdf",
        png = "Sample_{sample}/{sample}_barcode_distribution.png",
        mult_files = "Sample_{sample}/multi_fastq_files.fofn",
        demux_samples = "Sample_{sample}/demuxed_samples.txt"
    log:
        "logs/{sample}_demux_qc.log"
    threads:
        1
    conda:
        "../env/R.yml"
    script:
        "../scripts/extract_barcode_info.R"
       
checkpoint mv_mult_files:
    input:
        folder = rules.demux_minibar.output.folder,
        fofn = rules.demux_qc.output.mult_files
    output:
        bc_pass = directory("Sample_{sample}/demux_pass"),
        bc_fail = directory("Sample_{sample}/demux_fail")
    log:
        "logs/{sample}_mv_mult.log"
    threads:
        1
    shell:
        """
        set +e
        set +o pipefail
        mkdir -v -p {output.bc_fail} > {log}
        mkdir -v -p {output.bc_pass} >> {log}
        for file in $(cat {input.fofn})
        do
            mv -v "$file" {output.bc_fail} >>{log} 2>&1 || true
        done
        cp -r -v {input.folder}/* {output.bc_pass} >> {log}
        """

rule bc_map:
    input:
        fq = "Sample_{sample}/demux_pass/{sample}_{bc}.fastq",
        genome = config['ref']['genome']
    output:
        bam = "Sample_{sample}/demux_analysis/{bc}/{sample}_{bc}.bam",
        bai = "Sample_{sample}/demux_analysis/{bc}/{sample}_{bc}.bam.bai"
    conda:
        "../env/minimap2.yml"
    log:
        "logs/{sample}_{bc}_map.yml"
    threads:
        8
    shell:
        """
         minimap2 --MD -ax map-ont -t 4 \
            -R "@RG\\tID:{wildcards.bc}\\tSM:{wildcards.bc}" \
            {input.genome} {input.fq} 2> {log} \
            | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """

rule bc_bam_qc:
    input:
        "Sample_{sample}/demux_analysis/{bc}/{sample}_{bc}.bam"
    output:
        "Sample_{sample}/demux_analysis/{bc}/genome_results.txt"
    conda:
        "../env/ngs-bits.yml"
    log:
        "logs/{sample}_{bc}qualimap.log"
    threads:
        4
    params:
        qualimap = config['apps']['qualimap'],
        roi = config['demux']['target_region']
    shell:
        """
        {params.qualimap} bamqc \
            -bam {input} \
            --feature-file {params.roi} \
            --paint-chromosome-limits \
            -nt {threads} \
            -outdir Sample_{wildcards.sample}/demux_analysis/{wildcards.bc} \
            --java-mem-size=3G \
            >{log} 2>&1
        """

def get_demuxed_samples(wildcards):
        pass_folder = checkpoints.mv_mult_files.get(sample=wildcards.sample).output.bc_pass
        sample_files = expand("Sample_{sample}/demux_analysis/{bc_id}/genome_results.txt",
            sample = wildcards.sample,
            bc_id = glob_wildcards(os.path.join(pass_folder, wildcards.sample+"_{bc_id}.fastq" )).bc_id)
        #print(sample_files)
        return sample_files

rule coverage_stats:
    input:
        get_demuxed_samples
    output:
        "Sample_{sample}/coverage_stats.txt"
    shell:
        """
        echo {input} > {output}
        """