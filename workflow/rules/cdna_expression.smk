
configfile: srcdir('../../config/de_analysis_config.yml')

#_____ QUANTIFICATION (SALMON) ____________________________________________#

# Use Transcriptome alignment as input
rule quant_salmon:
    input:
        bam = rules.map_to_transcriptome.output.bam,
        trs = config['ref']['cDNA']
    output:
        counts = "Sample_{sample}/counts/quant.sf",
        stats = "Sample_{sample}/counts/aux_info/meta_info.json"
    conda:
        "../env/de_analysis.yml"
    log: 
        "logs/{sample}_salmon.log"
    threads:
        10
    params:
        library = config['libtype'], # U = unspecific, S = strand-specific
        args = "--noErrorModel" if config['ignore_error'] else ""
    shell:
        """
        salmon quant {params.args} \
            --threads {threads} \
            --libType {params.library} \
            --targets {input.trs} \
            --alignments {input.bam} \
            --output Sample_{wildcards.sample}/counts \
            > {log} 2>&1 
        """

# Merge Salmon counts
rule merge_counts:
    input:
        count_tsvs = expand("Sample_{sample}/counts/quant.sf", sample=ID_samples)
    output:
        tsv ="de_analysis/all_counts.tsv"
    conda:
        "../env/de_analysis.yml"
    params:
        script = srcdir("../scripts/merge_counts_tsv.py")
    shell:
        """
        python {params.script} -z -o {output.tsv} {input.count_tsvs}
        """

rule write_coldata:
    input:
    output:
        coldata = "de_analysis/coldata.tsv"
    run:
        samples, conditions, types = [], [], []
        for sample in config['control_samples'].keys():
            samples.append(config['control_samples'][sample])
            conditions.append("untreated")
            types.append("single-read")
        for sample in config['treated_samples'].keys():
            samples.append(config['treated_samples'][sample])
            conditions.append("treated")
            types.append("single-read")

        df = pd.DataFrame(OrderedDict([('sample', samples),('condition', conditions),('type', types)]))
        df.to_csv(output.coldata, sep="\t", index=False)

rule write_de_params:
    input:
    output:
        de_params = "de_analysis/de_params.tsv"
    run:
        d = OrderedDict()
        d["Annotation"] = [config['ref']["annotation"]]
        d["min_samps_gene_expr"] = [config["min_samps_gene_expr"]]
        d["min_samps_feature_expr"] = [config["min_samps_feature_expr"]]
        d["min_gene_expr"] = [config["min_gene_expr"]]
        d["min_feature_expr"] = [config["min_feature_expr"]]
        df = pd.DataFrame(d)
        df.to_csv(output.de_params, sep="\t", index=False)


rule de_analysis:
    input:
        de_params = rules.write_de_params.output.de_params,
        coldata = rules.write_coldata.output.coldata,
        tsv = rules.merge_counts.output.tsv,
    output:
        res_dge = "de_analysis/results_dge.tsv",
        pdf_dge = "de_analysis/results_dge.pdf",
        res_dtu_gene = "de_analysis/results_dtu_gene.tsv",
        res_dtu_trs = "de_analysis/results_dtu_transcript.tsv",
        res_dtu_stager = "de_analysis/results_dtu_stageR.tsv",
        flt_counts = "de_analysis/all_counts_filtered.tsv",
        flt_counts_gens = "de_analysis/all_gene_counts.tsv",
    conda:
        "../env/de_analysis.yml"
    params:
        script = srcdir("../scripts/de_analysis.R")
    shell:
        """
        Rscript {params.script}
        """

rule plot_dtu_res:
    input:
        res_dtu_stager = "de_analysis/results_dtu_stageR.tsv",
        flt_counts = "de_analysis/all_counts_filtered.tsv",
    output:
        dtu_pdf = "de_analysis/dtu_plots.pdf",
    conda: "../env/de_analysis.yml"
    params:
        script = srcdir("scripts/plot_dtu_results.R")
    shell:
        """
        {params.script}
        """
#_____ FEATURE COUNT (GENE LEVEL) _____________________________________________#

rule quant_genes:
    input:
        bam = "Sample_{sample}/{sample}.spliced.bam",
        ann = config['ref']['annotation']
    output:
        counts = "Sample_{sample}/{sample}.counts.tsv",
        stats = "Sample_{sample}/{sample}.counts.tsv.summary"
    conda:
        "../env/subread.yml"
    log:
        "logs/{sample}_featureCounts.log"
    threads:
        1
    params:
        level = config['expression']['level'],
        min_qual = config['expression']['min_qual'],
        strandedness = config['expression']['strandedness']
    shell:
        """
        featureCounts \
            -t {params.level} \
            -Q {params.min_qual} \
            -s {params.strandedness} \
            -L {input.bam} \
            -a {input.ann} \
            -o {output.counts} \
            > {log} 2>&1
        """