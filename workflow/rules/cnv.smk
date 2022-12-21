rule create_targets:
    input:
        bed = config['cnvkit']['targets']
    output:
        bed = "cnvkit/ref/targets.bed"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/cnvkit_targets.log"
    params:
        binsize = config['cnvkit']['binsize']
    shell:
        """
        cnvkit.py target \
            {input.bed} \
            --avg-size {params.binsize} \
            --split \
            --output {output.bed} \
            >{log} 2>&1
        """

rule create_antitargets:
    input:
        bed = config['cnvkit']['antitargets']
    output:
        bed = "cnvkit/ref/antitargets.bed"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/cnvkit_antitargets.log"
    params:
        binsize = config['cnvkit']['binsize'],
        minbinsize = config['cnvkit']['minbinsize']           
    shell:
        """
        cnvkit.py antitarget \
            {input.bed} \
            --avg-size {params.binsize} \
            --min-size {params.minbinsize} \
            --output {output.bed} \
            >{log} 2>&1
        """

# Calculate coverages for target bins
rule target_coverage:
    input:
        unpack(get_cnvkit_bam),
        bed = "cnvkit/ref/targets.bed"
    output:
        cnn = "cnvkit/{sample}/{sample}.targetcoverage.cnn"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/{sample}_cnvkit_target_coverage.log"
    threads:
        8
    shell:
        """
        cnvkit.py coverage \
            {input.bam} \
            {input.bed} \
            --processes {threads} \
            --output {output.cnn} \
            >{log} 2>&1
        """
 
# Calculate coverage for anti-target bins
rule antitarget_coverage:
    input:
        unpack(get_cnvkit_bam),
        bed = "cnvkit/ref/antitargets.bed"
    output:
        cnn = "cnvkit/{sample}/{sample}.antitargetcoverage.cnn"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/{sample}_cnvkit_antitarget_coverage.log"
    threads:
        8
    shell:
        """
        cnvkit.py coverage \
            {input.bam} \
            {input.bed} \
            --processes {threads} \
            --output {output.cnn} \
            >{log} 2>&1
        """

# Compile a copy-number reference from the given files or directory (containing normal samples)
rule cnvkit_reference:
    input:
        expand(
            "cnvkit/{sample}/{sample}.{type}.cnn",
            sample = config['cnvkit']['reference_samples'],
            type = ['targetcoverage', 'antitargetcoverage']
        )
    output:
        cnn = "cnvkit/ref/ref.cnn"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/cnvkit_ref_coverage.log"
    params:
        ref = config['ref']['genome']
    threads:
        8
    shell:
        """
        cnvkit.py reference \
            {input} \
            --fasta {params.ref} \
            --output {output.cnn} \
            >{log} 2>&1
        """

# Combine uncorrected target and antitarget coverage tables and correct for biases 
# according to reference. Output: copy number ratio table (.cnr)
rule fix_coverage:
    input:
        cnn_targets = rules.target_coverage.output.cnn,
        cnn_antitargets = rules.antitarget_coverage.output.cnn,
        ref = config['cnvkit']['reference'] if config['cnvkit']['reference'] else "cnvkit/ref/ref.cnn"
    output:
        cnr = "cnvkit/{sample}/{sample}.cnr"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/{sample}_cnvkit_fix.log"
    params:
        ref = config['cnvkit']['reference']
    shell:
        """
        cnvkit.py fix \
            --output {output.cnr} \
            {input.cnn_targets} \
            {input.cnn_antitargets} \
            {input.ref} \
            >{log} 2>&1
        """
# Infer discrete copy number segments from given coverage table
rule segment:
    input:
        cnr = rules.fix_coverage.output.cnr
    output:
        cns = "cnvkit/{sample}/{sample}.cns"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/{sample}_cnvkit_segment.log"
    params:
        method = config['cnvkit']['segment_method']
    threads:
        6
    shell:
        """
        cnvkit.py segment \
            {input.cnr} \
            --method {params.method} \
            --processes {threads} \
            --output {output.cns} \
            >{log} 2>&1
        """
# Given segmented log2 ratio estimates (.cns), derive each segment’s absolute integer copy number 
rule cnvkit_call:
    input:
        cns = rules.segment.output.cns
    output:
        call = "cnvkit/{sample}/{sample}.call.txt"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/{sample}_cnvkit_call.log"
    params:
        center = config['cnvkit']['call']['center'],
        method = config['cnvkit']['call']['method'],
        ploidy = config['cnvkit']['call']['ploidy']
    shell:
        """
        cnvkit.py call \
            {input.cns} \
            --center {params.center} \
            --method {params.method} \
            --ploidy {params.ploidy} \
            --output {output.call} \
            >{log} 2>&1
        """
# TODO: Annotate B-allele frequencies in existing SNP VCF

# Plot scatterplot of all segments
rule cnvkit_scatter:
    input:
        cns = rules.segment.output.cns,
        cnr = rules.fix_coverage.output.cnr
    output:
        "cnvkit/{sample}/{sample}.scatter.pdf"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/{sample}_cnvkit_scatter.log"
    params:
        col = config['cnvkit']['scatter']['color'],
        ymax = config['cnvkit']['scatter']['ymax'],
        ymin = config['cnvkit']['scatter']['ymin'],
    shell:
        """
        cnvkit.py scatter \
            --segment {input.cns} \
            --segment-color {params.col} \
            --y-max {params.ymax} \
            --y-min {params.ymin} \
            --output {output} \
            {input.cnr} \
            >{log} 2>&1
        """

rule convert_pdf:
    input:
        "cnvkit/{sample}/{sample}.{plot}.pdf"
    output:
        "cnvkit/{sample}/{sample}.{plot}.png"
    conda:
        "../env/cnvkit.yml"
    shell:
        """
        pdftoppm -png {input} cnvkit/{wildcards.sample}/$(basename {input} .pdf)
        mv cnvkit/{wildcards.sample}/{wildcards.sample}.{wildcards.plot}-1.png {output}
        """
