#_____ RUN READ QC  __________________________________________________________#

rule folder_pycoqc:
    input:
        lambda wildcards: glob(wildcards.folder + "/**/sequencing_summary*")
    output:
        html = "qc/pycoqc/per_run/{folder}.pycoQC.html",
        json = "qc/pycoqc/per_run/{folder}.pycoQC.json"
    log:
        "logs/{folder}_pycoqc.log"
    conda:
        "../env/pycoqc.yml"
    threads:
        1
    shell:
        """
        pycoQC \
            --summary_file {input}\
            --html_outfile {output.html} \
            --json_outfile {output.json} \
            >{log} 2>&1
        """

rule run_multiqc:
    input:
        expand("qc/pycoqc/per_run/{folder}.pycoQC.html", folder = ID_folders)
    output:
        "qc/pycoqc/per_run/run_multiqc_report.html"
    log:
        "logs/run_multiqc.log"
    threads:
       1
    params:
        multiqc = config['apps']['multiqc'],
        multiqc_config = srcdir('../../config/multiqc_config.yml')
    shell:
        """
        {params.multiqc} \
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
        lambda wildcards: glob(wildcards.folder + "/sequencing_summary*")
    output:
        directory("qc/pycoqc/split_{folder}")
    log:
        "logs/{folder}_pycoqc_split.log"
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
    output:
        "Sample_{sample}/sequencing_summary_bc_{sample}.txt"
    params:
        ssfile = lookup_split_summary_file
    shell:
        "cp {params.ssfile} {output}"

rule sample_pycoqc:
    input:
        unpack(get_summary_files),
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
        'qc/qualimap/{sample}_genome/genome_results.txt'
    log:
        "logs/{sample}_qualimap.log"
    threads:
        8
    params:
        qualimap = config['apps']['qualimap']
    shell:
        """
        {params.qualimap} bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            -nt {threads} \
            -outdir qc/qualimap/{wildcards.sample}_genome \
            --java-mem-size=12G \
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
        qualimap = config['apps']['qualimap']
    shell:
        """
        {params.qualimap} rnaseq \
            -bam {input} \
            -gtf {params.gtf} \
            -outdir qc/qualimap/{wildcards.sample}_rna \
            --java-mem-size=12G \
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

#____ VARIANT QC _____________________________________________________________#
rule bcftools_stats:
    input:
        vcf = rules.combine_vcf.output
    output:
        "qc/variants/{sample}.stats"
    log:
        "logs/{sample}_bcftools.log"
    conda:
        "../env/bcftools.yml"
    shell:
        """
        bcftools stats {input} > {output} 2> {log}
        """

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
        sqanti = config['apps']['sqanti'],
        cdna_cupcake = config['apps']['cdna_cupcake']
    shell:
        """
        set +u;
        export PYTHONPATH=$PYTHONPATH:{params.cdna_cupcake}sequence/
        export PYTHONPATH=$PYTHONPATH:{params.cdna_cupcake}
        tmp=$(mktemp)
        awk '$7!="." {{print $0}}' {input.gtf} > $tmp
        {params.sqanti} \
            -d qc/sqanti/{wildcards.sample}_{wildcards.tool} \
            -o {wildcards.sample}_{wildcards.tool} \
            --gtf $tmp \
            {input.anno_ref} {input.genome} \
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
    'mapping' : expand("qc/qualimap/{s}_genome/genome_results.txt", s = ID_samples),
    'assembly' : ["qc/quast_results/report.tsv"],
    'variant_calling':[expand("qc/variants/{s}.stats", s = ID_samples)],
    'structural_variant_calling' : [],
    'cDNA_stringtie' : expand("qc/gffcompare/{s}_stringtie/{s}_stringtie.stats", s = ID_samples) +
        expand("qc/pychopper/{s}_stats.txt", s = ID_samples), 
    'cDNA_flair': 
        expand("qc/rseqc/{s}.read_distribution.txt", s = ID_samples) + 
        expand("qc/rseqc/{s}.geneBodyCoverage.txt", s = ID_samples) +    
        expand("qc/gffcompare/{s}_flair/{s}_flair.stats", s = ID_samples),
    'cDNA_expression' : 
        #expand("qc/qualimap/{s}_rna/rnaseq_qc_results.txt", s = ID_samples) + 
        expand("qc/rseqc/{s}.read_distribution.txt", s = ID_samples) + 
        expand("qc/rseqc/{s}.geneBodyCoverage.txt", s = ID_samples) + 
        expand("qc/pychopper/{s}_stats.txt", s = ID_samples) +
        expand("Sample_{s}/{s}.counts.tsv.summary", s = ID_samples),
    'cDNA_pinfish' : [],
    'dual_demux' : [],
    'de_analysis' : [],
    'qc' : ["qc/pycoqc/per_run/run_multiqc_report.html",
        expand("qc/pycoqc/per_sample/{s}.pycoQC.json", s = ID_samples)],
}

# Additional output options
if config['vc']['create_benchmark']:
    qc_out['variant_calling'] +=  expand("qc/happy/{s}.summary.csv", s=ID_samples)

if map_samples_barcode: 
    qc_out += aggregate_sample_pycoqc

qc_out_selected = [qc_out[step] for step in config['steps']]

rule multiqc:
    input:
        [y for x in qc_out_selected for y in x]
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    threads:
        1
    params:
        multiqc = config['apps']['multiqc'],
        multiqc_config = srcdir('../../config/multiqc_config.yml')
    shell:
        """
        {params.multiqc} \
            --config  {params.multiqc_config} \
            --force \
            --outdir qc\
            --ignore-symlinks \
            {input} > {log} 2>&1
        """