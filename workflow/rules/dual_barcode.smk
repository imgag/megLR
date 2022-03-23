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
        bam = "Sample_{sample}/demux_analysis/{sample}_{bc}/{sample}_{bc}.bam",
        bai = "Sample_{sample}/demux_analysis/{sample}_{bc}/{sample}_{bc}.bam.bai"
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
        "Sample_{sample}/demux_analysis/{sample}_{bc}/{sample}_{bc}.bam"
    output:
        "Sample_{sample}/demux_analysis/{sample}_{bc}/genome_results.txt"
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
            --java-mem-size=12G \
            >{log} 2>&1
        """

rule bc_var:
    input:
        bam = rules.bc_map.output.bam,
        ref = config['ref']['genome']
    output:
        vcf = "Sample_{sample}/demux_analysis/{sample}_{bc}/round_1.vcf"
    conda:
        "../env/medaka.yml"
    threads:
        2
    log:
        "logs/{sample}_{bc}_medaka.log"
    params:
        model_snp = config['vc_medaka']['model_initial'],
        roi = config['demux']['target_region_str']
    shell:
        """
        export OMP_NUM_THREADS={threads}
        medaka_variant -d \
            -f {input.ref} \
            -i {input.bam} \
            -r {params.roi} \
            -p  \
            -o Sample_{wildcards.sample}/demux_analysis/{wildcards.sample}_{wildcards.bc}/ \
            -t {threads} \
            -s {params.model_snp} \
            -S \
            >{log} 2>&1
        """

rule vcf2tsv:
    input:
        rules.bc_var.output.vcf
    output:
        "Sample_{sample}/demux_analysis/{sample}_{bc}/round_1.tsv"
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t%GQ]\n' {input} > {output}
        """

def get_demuxed_samples(wildcards):
        pass_folder = checkpoints.mv_mult_files.get(sample=wildcards.sample).output.bc_pass
        sample_files = expand("Sample_{sample}/demux_analysis/{sample}_{bc_id}/genome_results.txt",
            sample = wildcards.sample,
            bc_id = glob_wildcards(os.path.join(pass_folder, wildcards.sample+"_{bc_id}.fastq" )).bc_id)
        #print(sample_files)
        return sample_files

def get_demuxed_variants_mm(wildcards):
        pass_folder = checkpoints.mv_mult_files.get(sample=wildcards.sample).output.bc_pass
        sample_files = expand("Sample_{sample}/demux_analysis/{sample}_{bc_id}/round_1.tsv",
            sample = wildcards.sample,
            bc_id = glob_wildcards(os.path.join(pass_folder, wildcards.sample+"_{bc_id}.fastq" )).bc_id)
        #print(sample_files)
        return sample_files

def get_demuxed_variants_dv(wildcards):
        pass_folder = checkpoints.mv_mult_files.get(sample=wildcards.sample).output.bc_pass
        sample_files = expand("Sample_{sample}/demux_analysis/{sample}_{bc_id}/{sample}_{bc_id}.dv.vcf.gz",
            sample = wildcards.sample,
            bc_id = glob_wildcards(os.path.join(pass_folder, wildcards.sample+"_{bc_id}.fastq" )).bc_id)
        #print(sample_files)
        return sample_files

rule all_demux:
    input:
        get_demuxed_samples,
        get_demuxed_variants_mm
    output:
        "Sample_{sample}/coverage_stats.txt"
    shell:
        """
        echo {input} > {output}
        """