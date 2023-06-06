configfile: srcdir('../../config/de_analysis_config.yml')

rule stringtie_merge:
    input:
        gtf =  expand("Sample_{sample}/stringtie/{sample}.stringtie.gtf", sample=ID_samples)
    output:
        gtf = "isoform_counts/all_samples_stringtie.isoforms.gtf"
    log:
        "logs/stringtie_merge.log"
    threads:
        2
    conda:
        "../env/stringtie.yml"
    params:
        annot = config['ref']['annotation']
    # Stringtie Merge Parameters:
    # -m : minimum transcript length (default: 50)
    # -c : minimum coverage (default: 0)
    # -F : minimum fpkm (default:  0)
    # -T : minimum tpm (default: 0)
    # -f : minimum isoform fraction (default: 0.01)
    # -i : keep merged transcripts with retained introns (default: false)
    # -G : Reference guided (requires input .GTF)
    shell:
        """
        stringtie --merge \
        -G {params.annot} \
        -o {output.gtf} \
        {input.gtf} \
        >{log} 2>&1
        """

rule create_stringtie_merged_transcriptome:
    input:
        gtf = "isoform_counts/all_samples_stringtie.isoforms.gtf",
        ref = config['ref']['genome']
    output:
        "isoform_counts/all_samples_stringtie.isoforms.fasta"
    log:
        "logs/convert_stringtie_transcriptome.log"
    threads:
        1
    conda:
        "../env/gffread.yml"
    params:
        rename_tab = srcdir("../../resources/chr_rename_table.tsv")
    shell:
        """
        gffread \
            -w {output} \
            -g {input.ref} \
            -m {params.rename_tab} \
            {input.gtf} \
            >{log} 2>&1
        """


#_____ QUANTIFICATION (SALMON) ____________________________________________#

# Use Transcriptome alignment as input
rule quant_salmon:
    input:
        bam = rules.map_to_transcriptome.output.bam,
        trs = get_fasta_annotation_for_transcriptome_alignment
    output:
        counts = "isoform_counts/{sample}_{method}_counts/quant.sf",
        stats = "isoform_counts/{sample}_{method}_counts/aux_info/meta_info.json"
    conda:
        "../env/de_analysis.yml"
    log: 
        "logs/{sample}_{method}_salmon.log"
    threads:
        10
    params:
        library = config['expression']['salmon_libtype'],
        args = "--noErrorModel" if config['ignore_error'] else ""
    shell:
        """
        salmon quant {params.args} \
            --threads {threads} \
            --libType {params.library} \
            --targets {input.trs} \
            --alignments {input.bam} \
            --output isoform_counts/{wildcards.sample}_{wildcards.method}_counts \
            > {log} 2>&1 
        """

# Merge Salmon counts
rule merge_counts:
    input:
        count_tsvs = expand("isoform_counts/{s}_stringtie-consensus_counts/quant.sf", s=ID_samples)
    output:
        tsv ="isoform_counts/merged_counts_stringtie.tsv"
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
        bam = get_cdna_bam
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
        strandedness = config['expression']['strandedness'],
        ann = config['ref']['annotation']
    shell:
        """
        featureCounts \
            -t {params.level} \
            -Q {params.min_qual} \
            -s {params.strandedness} \
            -L {input.bam} \
            -a {params.ann} \
            -o {output.counts} \
            > {log} 2>&1
        """
    
rule flairQuantify:
    input:
        fa = rules.flairCollapse.output.fasta,
        manifest = config['flair']['manifest']
    output:
        mat = "Sample_{sample}/flair/{sample}.counts_matrix.tsv"
    conda:
        "../env/flair.yml"
    threads:
        20
    log:
        "logs/{sample}_flair_quantify.log"
    params:
        trust_ends = "--trust_ends" if config['flair']['trust_ends'] else "",
        qual = config['flair']['min_qual']
    shell:
        """ 
        flair.py quantify \
            --quality {params.qual} \
            {params.trust_ends} \
            --output {output} \
            --tpm \
            --reads_manifest {input.manifest} \
            --isoforms {input.fa} \
            --threads {threads} \
            >{log} 2>&1
        """

rule merge_stringtie_counts:
    input:
        abund = "Sample_{sample}/stringtie/{sample}.stringtie.abundance.tsv",
        classif = "qc/sqanti/{sample}_stringtie/{sample}_stringtie_classification.txt"
    output:
        tsv = "Sample_{sample}/{sample}.stringtie.annotated_expression.tsv"
    conda:
        "../env/R.yml"
    log:
        "logs/{sample}_stringtie_annotation.log"
    script:
        "../scripts/annotate_stringtie_abundance.R"
