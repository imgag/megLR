#_____ POREC ANALYSIS ____________________________________________________________#
#
# PoreC analysis pipeline adapted from imgag/T2T_ONT.
# Steps:
#   1. Falign mapping (build repeat regions, map reads)
#   2. Contact generation  (pairtools parse2, sort/flip, stats, HTML report)
#   3. HiC generation      (cooler, mcool, balance, juicer .hic)
#   4. TAD calling         (hicFindTADs via HiCExplorer)
#   5. QC & region plots   (hicCorrectMatrix, hicPlotDistVsCounts, hicPlotMatrix)
#
# Config keys used from config['porec']:
#   falign               - path to Falign binary
#   juicer_tools         - path to juicer_tools JAR (set to "" to skip .hic creation)
#   assembly_name        - genome assembly name embedded in pairs headers
#   min_bin_width        - finest cooler resolution (bp)
#   cooler_resolutions   - comma-separated list for mcool zoom levels
#   juicer_resolutions   - comma-separated list for .hic file
#   porec_resolutions    - list of resolutions for QC plots and TAD calling inputs
#   tad_resolutions      - list of resolutions for hicFindTADs
#   tad_min_depth        - hicFindTADs minDepth
#   tad_max_depth        - hicFindTADs maxDepth
#   tad_step             - hicFindTADs step
#   tad_fdr_threshold    - hicFindTADs FDR threshold
#   tad_delta            - hicFindTADs delta
#   tad_correction_threshold - hicFindTADs correction factor threshold
#   plot_regions         - list of genomic regions for contact map plots,
#                          e.g. ["chr11:1000000-3000000"]  (empty = skip)
#   plot_resolution      - resolution used for hicPlotMatrix target region plots

# ---------------------------------------------------------------------------
# Helper: decide whether to generate .hic (juicer_tools path must be set)
# ---------------------------------------------------------------------------
def _porec_juicer_enabled():
    return bool(config['porec'].get('juicer_tools', ''))


# ---------------------------------------------------------------------------
# Chromosome sizes from reference FAI
# ---------------------------------------------------------------------------
rule porec_create_chromsizes:
    """Derive chromosome size file from the reference FASTA index."""
    input:
        fai = ancient(config['ref']['genome'] + ".fai")
    output:
        "porec/ref.chrom.sizes"
    run:
        with open(input.fai) as fin, open(output[0], 'w') as fout:
            for line in fin:
                parts = line.split('\t')
                fout.write(f"{parts[0]}\t{parts[1]}\n")

# ---------------------------------------------------------------------------
# Step 0: Build repeat regions for Falign (once per reference)
# ---------------------------------------------------------------------------
rule porec_build_repeat_regions:
    """Build repeat-region BED used by Falign to avoid spurious contacts."""
    input:
        ref = ancient(config['ref']['genome'])
    output:
        bed = "porec/ref.repeat_regions.bed"
    params:
        falign = config['porec']['falign']
    log:
        "logs/porec/build_repeat_regions.log"
    shell:
        """
        {params.falign} build-repeat \
            {input.ref} \
            {output.bed} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 1: Map PoreC reads with Falign
# ---------------------------------------------------------------------------
rule porec_falign_map:
    """Map PoreC reads to reference with Falign in fragment-BAM mode."""
    input:
        fq  = "Sample_{sample}/{sample}.fastq.gz",
        fa  = ancient(config['ref']['genome']),
        bed = "porec/ref.repeat_regions.bed"
    output:
        bam = "porec/{sample}/1-falign/{sample}.fragments.bam"
    params:
        falign = config['porec']['falign']
    threads: 30
    log:
        "logs/porec/falign_map.{sample}.log"
    shell:
        """
        {params.falign} \
            -repeat_bed {input.bed} \
            -num_threads {threads} \
            -outfmt frag-bam \
            -out {output.bam} \
            {input.fa} \
            {input.fq} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 2a: Parse fragment BAM → pairs with pairtools parse2
# ---------------------------------------------------------------------------
rule porec_pairtools_parse2:
    """Convert Falign fragment BAM to pairtools pairs format."""
    input:
        bam       = "porec/{sample}/1-falign/{sample}.fragments.bam",
        chromsize = "porec/ref.chrom.sizes"
    output:
        pairs = temp("porec/{sample}/pairs/{sample}.raw.pairs.gz")
    params:
        assembly    = config['porec'].get('assembly_name', 'genome'),
        orientation = "pair",
        position    = "junction"
    log:
        "logs/porec/pairtools_parse2.{sample}.log"
    conda:
        "../env/pairtools.yml"
    shell:
        """
        pairtools parse2 \
            --chroms-path {input.chromsize} \
            --assembly {params.assembly} \
            --report-position {params.position} \
            --report-orientation {params.orientation} \
            --add-pair-index \
            --single-end \
            --expand \
            --flip \
            --readid-transform 'readID.split(":")[0]' \
            --drop-seq \
            --drop-sam \
            --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches \
            --output {output.pairs} \
            {input.bam} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 2b: Sort and flip pairs
# ---------------------------------------------------------------------------
rule porec_sort_pairs:
    """Flip (canonicalise strand orientation) and sort the pairs file."""
    input:
        pairs     = "porec/{sample}/pairs/{sample}.raw.pairs.gz",
        chromsize = "porec/ref.chrom.sizes"
    output:
        "porec/{sample}/pairs/{sample}.pairs.gz"
    threads: 8
    resources:
        mem_gb = 30
    log:
        "logs/porec/sort_pairs.{sample}.log"
    conda:
        "../env/pairtools.yml"
    shell:
        """
        pairtools flip --chroms-path {input.chromsize} {input.pairs} | \
        pairtools sort \
            --nproc {threads} \
            --memory 25G \
            --output {output} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 2c: Pairs statistics
# ---------------------------------------------------------------------------
rule porec_pairtools_stats:
    """Compute summary statistics for the pairs file."""
    input:
        "porec/{sample}/pairs/{sample}.pairs.gz"
    output:
        "porec/{sample}/pairs/{sample}.pairs.stats.txt"
    log:
        "logs/porec/pairtools_stats.{sample}.log"
    conda:
        "../env/pairtools.yml"
    shell:
        """
        pairtools stats \
            --output {output} \
            {input} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 2d: HTML pairs statistics report
# ---------------------------------------------------------------------------
rule porec_pairs_stats_report:
    """Generate an interactive HTML report from the pairs statistics."""
    input:
        "porec/{sample}/pairs/{sample}.pairs.stats.txt"
    output:
        "porec/{sample}/pairs/{sample}.pairs.stats.html"
    params:
        report_script = workflow.source_path("../scripts/create_pairs_report.py")
    log:
        "logs/porec/pairs_stats_report.{sample}.log"
    conda:
        "../env/py_report.yml"
    shell:
        """
        python {params.report_script} \
            {input} {output} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 3a: Convert pairs to cooler at the finest resolution
# ---------------------------------------------------------------------------
rule porec_pairs_to_cooler:
    """Build a single-resolution cooler file from the pairs."""
    input:
        fai   = ancient(config['ref']['genome'] + ".fai"),
        pairs = "porec/{sample}/pairs/{sample}.pairs.gz"
    output:
        cool = "porec/{sample}/cooler/{sample}_{resolution}.cool"
    wildcard_constraints:
        resolution = r"\d+"
    threads: 2
    resources:
        mem_gb = 64
    log:
        "logs/porec/pairs_to_cooler.{sample}_{resolution}.log"
    conda:
        "../env/cooler.yml"
    shell:
        """
        cooler cload pairs \
            -c1 2 -p1 3 -c2 4 -p2 5 \
            {input.fai}:{wildcards.resolution} \
            {input.pairs} \
            {output.cool} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 3b: Zoomify to multi-resolution .mcool
# ---------------------------------------------------------------------------
rule porec_cooler_zoomify:
    """Create a multi-resolution .mcool from the base-resolution cooler."""
    input:
        expand(
            "porec/{{sample}}/cooler/{{sample}}_{resolution}.cool",
            resolution=config['porec'].get('min_bin_width', 1000)
        )
    output:
        mcool = "porec/{sample}/cooler/{sample}.mcool"
    params:
        resolutions = config['porec'].get('cooler_resolutions', '1000,5000,10000,50000,100000')
    threads: 2
    resources:
        mem_gb = 64
    log:
        "logs/porec/cooler_zoomify.{sample}.log"
    conda:
        "../env/cooler.yml"
    shell:
        """
        cooler zoomify \
            -r {params.resolutions} \
            -o {output.mcool} \
            {input} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 3c: Balance each resolution cooler (ICE normalisation)
# ---------------------------------------------------------------------------
rule porec_cooler_balance:
    """ICE-balance a single-resolution cooler (in-place copy)."""
    input:
        "porec/{sample}/cooler/{sample}_{resolution}.cool"
    output:
        "porec/{sample}/cooler/{sample}_{resolution}_balanced.cool"
    threads: 16
    resources:
        mem_gb = 64
    log:
        "logs/porec/cooler_balance.{sample}_{resolution}.log"
    conda:
        "../env/cooler.yml"
    shell:
        """
        cp {input} {output}
        cooler balance \
            --cis-only \
            --ignore-diags 2 \
            --mad-max 5 \
            --min-nnz 10 \
            --tol 1e-05 \
            --max-iters 500 \
            --nproc {threads} \
            {output} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 3d: Convert pairs to .hic (Juicebox) – optional
# ---------------------------------------------------------------------------
rule porec_clean_pairs_for_hic:
    """Reformat pairs to the MEDIUM format expected by juicer_tools pre."""
    input:
        "porec/{sample}/pairs/{sample}.pairs.gz"
    output:
        temp("porec/{sample}/pairs/{sample}.pairs.for_juice")
    resources:
        mem_gb = 64
    shell:
        """
        zcat {input} | \
        awk '!/^#/ {{OFS="\\t"; \
        strand1 = ($6=="+")?0:1; \
        strand2 = ($7=="+")?0:1; \
        print $1,strand1,$2,$3,0,strand2,$4,$5,1,$11,$12}}' \
        > {output}
        """

rule porec_pairs_to_hic:
    """Generate a .hic contact file for visualisation in Juicebox."""
    input:
        pairs     = "porec/{sample}/pairs/{sample}.pairs.for_juice",
        chromsize = "porec/ref.chrom.sizes"
    output:
        hic = "porec/{sample}/hic/{sample}.hic"
    params:
        juicer_tools = config['porec']['juicer_tools'],
        resolutions  = config['porec'].get('juicer_resolutions', '5000,10000,25000,50000,100000,500000')
    log:
        "logs/porec/pairs_to_hic.{sample}.log"
    resources:
        mem_gb = 128
    shell:
        """
        java -Xmx{resources.mem_gb}G -jar {params.juicer_tools} \
            pre {input.pairs} {output.hic} {input.chromsize} \
            -r {params.resolutions} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 4: TAD calling with HiCExplorer
# ---------------------------------------------------------------------------
rule porec_hic_find_tads:
    """Call TADs using hicFindTADs on a balanced cooler."""
    input:
        cool = "porec/{sample}/cooler/{sample}_{resolution}_balanced.cool"
    output:
        tads       = "porec/{sample}/tad/{sample}_{resolution}_domains.bed",
        boundaries = "porec/{sample}/tad/{sample}_{resolution}_boundaries.bed",
        score      = "porec/{sample}/tad/{sample}_{resolution}_score.bedgraph"
    params:
        min_depth  = lambda wc: max(
            config['porec'].get('tad_min_depth', 20000),
            3 * int(wc.resolution)
        ),
        max_depth  = config['porec'].get('tad_max_depth', 200000),
        step       = lambda wc: max(
            config['porec'].get('tad_step', 10000),
            int(wc.resolution)
        ),
        fdr        = config['porec'].get('tad_fdr_threshold', 0.05),
        delta      = config['porec'].get('tad_delta', 0.01),
        correction = config['porec'].get('tad_correction_threshold', 1.5)
    threads: 4
    log:
        "logs/porec/hic_find_tads.{sample}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicFindTADs \
            --matrix {input.cool} \
            --outPrefix porec/{wildcards.sample}/tad/{wildcards.sample}_{wildcards.resolution} \
            --minDepth {params.min_depth} \
            --maxDepth {params.max_depth} \
            --step {params.step} \
            --thresholdComparisons {params.fdr} \
            --delta {params.delta} \
            --correctForMultipleTesting fdr \
            --numberOfProcessors {threads} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 5a: Diagnostic plot (matrix quality check)
# ---------------------------------------------------------------------------
rule porec_hic_diagnostic_plot:
    """Generate a diagnostic plot to evaluate matrix quality before balancing."""
    input:
        cool = "porec/{sample}/cooler/{sample}_{resolution}_balanced.cool"
    output:
        "porec/{sample}/qc/{sample}_{resolution}_diagnostic.png"
    resources:
        mem_gb = 64
    log:
        "logs/porec/hic_diagnostic.{sample}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicCorrectMatrix diagnostic_plot \
            --matrix {input.cool} \
            -o {output} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 5b: Distance-vs-counts plot
# ---------------------------------------------------------------------------
rule porec_hic_plot_dist_vs_counts:
    """Plot contact frequency as a function of genomic distance."""
    input:
        cool = "porec/{sample}/cooler/{sample}_{resolution}_balanced.cool"
    output:
        "porec/{sample}/qc/plot_vs_counts_{resolution}.png"
    resources:
        mem_gb = 64
    log:
        "logs/porec/plot_dist_vs_counts.{sample}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicPlotDistVsCounts \
            --matrices {input.cool} \
            -o {output} \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Step 5c: Contact map plots for user-defined target regions
# ---------------------------------------------------------------------------
rule porec_hic_plot_matrix:
    """Plot a contact matrix for a specific genomic target region."""
    input:
        cool = "porec/{sample}/cooler/{sample}_{resolution}_balanced.cool"
    output:
        "porec/{sample}/plots/{sample}_{resolution}_{region}.png"
    wildcard_constraints:
        region = r"[^/]+"
    params:
        region_fmt = lambda wc: wc.region.replace("_", ":")
    log:
        "logs/porec/hic_plot_matrix.{sample}.{resolution}.{region}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicPlotMatrix \
            --matrix {input.cool} \
            --outFileName {output} \
            --region {params.region_fmt} \
            --log1p \
            >{log} 2>&1
        """

# ---------------------------------------------------------------------------
# Collect: aggregate all PoreC outputs for a sample
# ---------------------------------------------------------------------------
def _porec_outputs(sample):
    """Return the list of all expected PoreC output files for *sample*."""
    s = sample
    porec_res   = config['porec'].get('porec_resolutions', ['5000', '25000'])
    tad_res     = config['porec'].get('tad_resolutions',   ['25000'])
    plot_res    = config['porec'].get('plot_resolution',   '25000')
    plot_regions = config['porec'].get('plot_regions', [])

    outs = []

    # Pairs + stats
    outs += [f"porec/{s}/pairs/{s}.pairs.gz"]
    outs += [f"porec/{s}/pairs/{s}.pairs.stats.html"]

    # Multi-resolution mcool
    outs += [f"porec/{s}/cooler/{s}.mcool"]

    # Per-resolution balanced coolers
    outs += expand(
        "porec/{s}/cooler/{s}_{res}_balanced.cool",
        s=[s], res=porec_res
    )

    # QC plots
    outs += expand(
        "porec/{s}/qc/{s}_{res}_diagnostic.png",
        s=[s], res=porec_res
    )
    outs += expand(
        "porec/{s}/qc/plot_vs_counts_{res}.png",
        s=[s], res=porec_res
    )

    # TAD domains
    outs += expand(
        "porec/{s}/tad/{s}_{res}_domains.bed",
        s=[s], res=tad_res
    )

    # Optional .hic file
    if _porec_juicer_enabled():
        outs += [f"porec/{s}/hic/{s}.hic"]

    # Optional target-region plots
    for region in plot_regions:
        region_safe = region.replace(":", "_")
        outs.append(f"porec/{s}/plots/{s}_{plot_res}_{region_safe}.png")

    return outs


rule porec_collect:
    """Sentinel rule: touch a .done file once all PoreC outputs are ready."""
    input:
        lambda wc: _porec_outputs(wc.sample)
    output:
        "porec/{sample}/porec.done"
    shell:
        "touch {output}"
