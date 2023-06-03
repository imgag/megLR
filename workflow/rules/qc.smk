#_____ RUN READ QC  __________________________________________________________#

rule folder_pycoqc:
    input:
        summarystats = lambda wildcards: [x for y in [glob(r + "/**/sequencing_summary*.txt") for r in map_runs_folder[wildcards.run]] for x in y]
    output:
        html = "qc/pycoqc/per_run/{run}/{run}.pycoQC.html",
        json = "qc/pycoqc/per_run/{run}/{run}.pycoQC.json"
    log:
        "logs/{run}_pycoqc.log"
    conda:
        "../env/pycoqc.yml"
    threads:
        1
    shell:
        """
        pycoQC \
            --summary_file {input.summarystats}\
            --html_outfile {output.html} \
            --json_outfile {output.json} \
            >{log} 2>&1
        """

rule run_multiqc:
    input:
        expand("qc/pycoqc/per_run/{run}/{run}.pycoQC.html", run = ID_runs)
    output:
        "qc/pycoqc/per_run/run_multiqc_report.html"
    log:
        "logs/run_multiqc.log"
    conda:
        "../env/multiqc.yml"
    threads:
       1
    params:
        multiqc_config = srcdir('../../config/multiqc_config.yml')
    shell:
        """
        multiqc \
            --force \
            --outdir qc/pycoqc/per_run \
            --config  {params.multiqc_config} \
            --filename run_multiqc_report \
            qc/pycoqc/per_run/ \
            >{log} 2>> {log}
        """

#_____ SAMPLE READ QC  _________________________________________________________#


checkpoint split_summary_perbarcode:
    input:
        get_summary_files_to_split
    output:
        directory("qc/pycoqc/split_barcodes_{run}")
    log:
        "logs/pycoqc_split_{run}.log"
    conda:
        "../env/pycoqc.yml"
    shell:
        """
        mkdir -p {output}
        Barcode_split \
            --summary_file {input} \
            --output_dir {output} \
            --output_unclassified \
            --verbose >{log} 2>&1
        """

rule rename_split_summary_files:
    input:
        lookup_split_summary_file
    output:
        "Sample_{sample}/sequencing_summary_bc_{sample}.txt"
    wildcard_constraints:
        folder = "^\w+$"
    shell:
        "cp {input} {output}"

rule sample_pycoqc:
    input:
        unpack(get_summary_files_sample),
    output:
        html = "qc/pycoqc/per_sample/{sample}.pycoQC.html",
        json = "qc/pycoqc/per_sample/{sample}.pycoQC.json"
    conda:
        "../env/pycoqc.yml"
    log:
        "logs/{sample}_pycoqc.log"
    threads:
        1
    shell:
        """
        pycoQC \
            --summary_file {input.summary_files} \
            --html_outfile {output.html} \
            --json_outfile {output.json} \
            >{log} 2>&1
        """

#____ GENOME MAPPING QC __________________________________________________________#

rule qualimap:
    input:
        "Sample_{sample}/{sample}.bam"
    output:
        directory('qc/qualimap/{sample}_genome')
    log:
        "logs/{sample}_qualimap.log"
    threads:
        8
    conda:
        "../env/qualimap.yml"
    params:
        target = "--outside-stats --feature-file " + config['qualimap']['target'] if config['qualimap']['target'] else ""

    shell:
        """
        unset DISPLAY
        qualimap bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            -nt {threads} {params.target}\
            -outdir qc/qualimap/{wildcards.sample}_genome \
            --java-mem-size=24G \
            >{log} 2>&1
        """
    

rule qualimap_mod:
    input:
        "Sample_{sample}/{sample}.mod.bam"
    output:
        directory('qc/qualimap/{sample}_modbases')
    log:
        "logs/{sample}_qualimap.log"
    threads:
        8
    conda:
        "../env/qualimap.yml"
    shell:
        """
        unset DISPLAY
        qualimap bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            -nt {threads} \
            -outdir qc/qualimap/{wildcards.sample}_modbases \
            --java-mem-size=24G \
            >{log} 2>&1
        """

#_____ cDNA SPLICED MAPPING QC _________________________________________________#

rule rna_qualimap:
    input:
        "Sample_{sample}/{sample}.spliced.bam"
    output:
        "qc/qualimap/{sample}_rna/rnaseq_qc_results.txt"
    log:
        "logs/{sample}_qualimap_rna.log"
    threads:
        2
    params:
        gtf = config['ref']['annotation'],
    conda:
        "../env/qualimap.yml"
    shell:
        """
        unset DISPLAY
        qualimap rnaseq \
            -bam {input} \
            -gtf {params.gtf} \
            -outdir qc/qualimap/{wildcards.sample}_rna \
            --java-mem-size=24G \
            >{log} 2>&1
        """

rule rseqc_read_distribution:
    input:
        bam = "Sample_{sample}/{sample}.spliced.bam",
        annot = config['ref']['annotation_bed']
    output:
        "qc/rseqc/{sample}.read_distribution.txt"
    log:
        "logs/{sample}_rseqc_distribution.log"
    threads:
        2
    conda:
        "../env/rseqc.yml"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.annot} > {output} 2> {log}
        """

rule rseqc_geneBody_coverage:
    input:
        bam = "Sample_{sample}/{sample}.spliced.bam",
        annot = config['ref']['annotation_bed']
    output:
        "qc/rseqc/{sample}.geneBodyCoverage.txt"
    log:
        "logs/{sample}_rseqc_coverage.log"
    threads:
        2
    conda:
        "../env/rseqc.yml"
    shell:
        """
        geneBody_coverage.py -i {input.bam} -r {input.annot} -o qc/rseqc/{wildcards.sample} >{log} 2>&1
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
            --output-dir qc/quast_results \
            {input} \
            >{log} 2>&1
        """

# #____ VARIANT QC _____________________________________________________________#
# rule bcftools_stats:
#     input:
#         vcf = rules.combine_vcf.output
#     output:
#         "qc/variants/{sample}.stats"
#     log:
#         "logs/{sample}_bcftools.log"
#     conda:
#         "../env/bcftools.yml"
#     shell:
#         """
#         bcftools stats {input} > {output} 2> {log}
#         """

#_____ TRANSCRIPTOME ANNOTATION QC ____________________________________________#
# Exons with ambigous strand informations are filtered out to avoid sqanti errors
# This happens only with stringtie annotation

rule sqanti:
    input:
        gtf = "Sample_{sample}/{sample}.{tool}.gtf",
        anno_ref = config['ref']['annotation'],
        genome = config['ref']['genome']
    output:
        "qc/sqanti/{sample}_{tool}/{sample}_{tool}_classification.txt"
    log:
        "logs/{sample}_{tool}_sqanti.log"
    conda:
        "../env/sqanti.yml"
    params:
        sqanti = config['apps']['sqanti']
    shell:
        """
        set +u;
        tmp=$(mktemp)
        awk '$7!="." {{print $0}}' {input.gtf} > $tmp
        {params.sqanti} \
            -d qc/sqanti/{wildcards.sample}_{wildcards.tool} \
            -o {wildcards.sample}_{wildcards.tool} \
            --report both \
            $tmp {input.anno_ref} {input.genome} \
            >{log} 2>&1
        """

#_______ GFFCOMPARE ___________________________________________________________#

rule gffcompare:
    input:
        gtf = "Sample_{sample}/{sample}.{tool}.gtf",
        ref = config['ref']['annotation'],
        genome = config['ref']['genome']
    output:
        "qc/gffcompare/{sample}_{tool}/{sample}_{tool}.stats"
    log:
        "logs/{sample}_{tool}.log"
    conda:
        "../env/gffcompare.yml"
    params:
        opts = config['gffcompare']['options']
    shell:
        """
        gffcompare \
            {params.opts} \
            -r {input.ref} \
            -o qc/gffcompare/{wildcards.sample}_{wildcards.tool}/{wildcards.sample}_{wildcards.tool} \
            {input.gtf} \
            > {log} 2>&1
        mv Sample_{wildcards.sample}/{wildcards.sample}_{wildcards.tool}.{wildcards.sample}.{wildcards.tool}.gtf.refmap qc/gffcompare/{wildcards.sample}_{wildcards.tool}/{wildcards.sample}.{wildcards.tool}.gtf.refmap
        mv Sample_{wildcards.sample}/{wildcards.sample}_{wildcards.tool}.{wildcards.sample}.{wildcards.tool}.gtf.tmap qc/gffcompare/{wildcards.sample}_{wildcards.tool}/{wildcards.sample}.{wildcards.tool}.gtf.tmap
        """

#_____ MULTI QC  _____________________________________________________________#

qc_out = {
    'fastq' : expand("qc/pycoqc/per_sample/{s}.pycoQC.json", s = ID_samples),
    'mapping' : expand("qc/qualimap/{s}_genome", s = ID_samples),
    'assembly' : ["qc/quast_results/report.tsv"],
    'modbases' : expand("qc/qualimap/{s}_modbases", s = ID_samples),
    'modification' : [],
    'duplex' : [],
    'variant_calling':[],
    'structural_variant_calling' : [],
    'cdna' : 
        expand("qc/pychopper/{s}_stats.txt", s = ID_samples) + 
        expand("qc/rseqc/{s}.read_distribution.txt", s = ID_samples) + 
        expand("qc/rseqc/{s}.geneBodyCoverage.txt", s = ID_samples),
    'transcriptome' : expand("qc/gffcompare/{s}_{t}/{s}_{t}.stats", 
        s = ID_samples, t = config['transcriptome']['methods']),
    'expression' : 
        expand("Sample_{s}/{s}.counts.tsv.summary", s = ID_samples),
    'dual_demux' : [],
    'de_analysis' : [],
    'cnv': [],
    'local_assembly': [],
    'repeat_expansion': [],
    'qc' : ["qc/pycoqc/per_run/run_multiqc_report.html"],
    'qc_db': []
}

if not config['disable_sampleqc']:
    qc_out['qc'] += expand("qc/pycoqc/per_sample/{s}.pycoQC.json", s = ID_samples)


if config['cdna']['with_umi'] and config['cdna']['deduplicate']:
    qc_out['cdna'] += expand("qc/umitools_dedup/{s}_stats_per_umi_per.tsv", s = ID_samples)

qc_out_selected = [qc_out[step] for step in config['steps']]

rule multiqc:
    input:
        [y for x in qc_out_selected for y in x]
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "../env/multiqc.yml"
    threads:
        1
    params:
        multiqc_config = srcdir('../../config/multiqc_config.yml')
    shell:
        """
        multiqc \
            --config  {params.multiqc_config} \
            --force \
            --outdir qc\
            --ignore-symlinks \
            {input} > {log} 2>&1
        """
