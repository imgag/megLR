#____ PINFISH PIPELINE (ONT) ______________________________________#

rule convert_bam: ## convert BAM to GFF
    input:
        bam = rules.map_genome_splice.output.bam
    output:
        raw_gff = "Sample_{sample}/pinfish/raw_transcripts.gff"
    conda: 
        "../env/pinfish.yml"
    params:
        opts = config['pinfish']["spliced_bam2gff_opts"]
    threads: 12
    shell:
        """
        spliced_bam2gff {params.opts} -t {threads} -M {input.bam} > {output.raw_gff}
        """

rule cluster_gff: ## cluster transcripts in GFF
    input:
        raw_gff = rules.convert_bam.output.raw_gff
    output:
        cls_gff = "Sample_{sample}/pinfish/clustered_transcripts.gff",
        cls_tab = "Sample_{sample}/pinfish/cluster_memberships.tsv",
    params:
        c = config['pinfish']["minimum_cluster_size"],
        d = config['pinfish']["exon_boundary_tolerance"],
        e = config['pinfish']["terminal_exon_boundary_tolerance"],
        min_iso_frac = config['pinfish']["minimum_isoform_percent"],
    conda: 
        "../env/pinfish.yml"
    # This step does not seem to benefit from multithreading --> very slow
    threads: 3
    shell:
        """
        cluster_gff \
            -p {params.min_iso_frac} \
            -t {threads} \
            -c {params.c} \
            -d {params.d} \
            -e {params.e} \
            -a {output.cls_tab} \
            {input.raw_gff} > {output.cls_gff}
        """

rule collapse_clustered: ## collapse clustered read artifacts
    input:
        cls_gff = rules.cluster_gff.output.cls_gff
    output:
        cls_gff_col = "Sample_{sample}/pinfish/clustered_transcripts_collapsed.gff"
    params:
        d = config['pinfish']["collapse_internal_tol"],
        e = config['pinfish']["collapse_three_tol"],
        f = config['pinfish']["collapse_five_tol"],
    conda: 
        "../env/pinfish.yml"
    shell:
        """
        collapse_partials \
            -d {params.d} \
            -e {params.e} \
            -f {params.f} \
            {input.cls_gff} > {output.cls_gff_col}
        """

rule polish_clusters: ## polish read clusters
    input:
        cls_gff = rules.cluster_gff.output.cls_gff,
        cls_tab = rules.cluster_gff.output.cls_tab,
        bam = rules.map_genome_splice.output.bam
    output:
        pol_trs = "Sample_{sample}/pinfish/polished_transcripts.fas",
    params:
        c = config['pinfish']["minimum_cluster_size"],
    conda: 
        "../env/pinfish.yml"
    threads:
        12
    shell:
        """
        polish_clusters \
            -t {threads} \
            -a {input.cls_tab} \
            -c {params.c} \
            -o {output.pol_trs} \
            {input.bam}
        """

rule map_polished: ## map polished transcripts to genome
    input:
       fasta = rules.polish_clusters.output.pol_trs,
    output:
       pol_bam = "Sample_{sample}/pinfish/polished_reads_aln_sorted.bam"
    params:
        extra = config['pinfish']["minimap2_opts_polished"],
        genome = config['ref']['genome']
    conda: 
        "../env/pinfish.yml"
    threads: 12
    shell:
        """
        minimap2 -t {threads} {params.extra} -ax splice {params.genome} {input.fasta}\
        | samtools view -Sb -F 2304 | samtools sort -@ {threads} - -o {output.pol_bam};
        samtools index {output.pol_bam}
        """

rule convert_polished: ## convert BAM of polished transcripts to GFF
    input:
        bam = rules.map_polished.output.pol_bam
    output:
        pol_gff = "Sample_{sample}/pinfish/polished_transcripts.gff"
    params:
        extra = config['pinfish']["spliced_bam2gff_opts_pol"]
    conda: 
        "../env/pinfish.yml"
    threads: 12
    shell:
        """
        spliced_bam2gff {params.extra} \
            -t {threads} \
            -M {input.bam} > {output.pol_gff}
        """

rule collapse_polished: ## collapse polished read artifacts
    input:
        pol_gff = rules.convert_polished.output.pol_gff
    output:
        pol_gff_col = "Sample_{sample}/pinfish/polished_transcripts_collapsed.gff"
    params:
        d = config['pinfish']["collapse_internal_tol"],
        e = config['pinfish']["collapse_three_tol"],
        f = config['pinfish']["collapse_five_tol"],
    conda: 
        "../env/pinfish.yml"
    shell:
        """
        collapse_partials \
            -d {params.d} \
            -e {params.e} \
            -f {params.f} \
            {input.pol_gff} > {output.pol_gff_col}
        """

rule gen_corr_trs: ## Generate corrected transcriptome.
    input:
        genome = config["ref"]['genome'],
        gff = rules.collapse_polished.output.pol_gff_col,
    output:
        fasta = "Sample_{sample}/pinfish/corrected_transcriptome_polished_collapsed.fasta"
    conda: 
        "../env/pinfish.yml"
    shell:
        """
        gffread -g {input.genome} -w {output.fasta} {input.gff}
        """

