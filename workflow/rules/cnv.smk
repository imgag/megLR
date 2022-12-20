

# Calculate coverages for target bins
rule target_coverage:
    input:
        bam = "Sample_{sample}/{sample}.bam"
    output:
        cnn = "cnvkit/{sample}/{sample}.targetcoverage.cnn"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/{sample}_cnvkit_target_coverage.log"
    params:
        bed = config['cnvkit']['targets']
    shell:
        """
        cnvkit.py coverage \
            {input.bam} \
            {params.bed} \
            --output {output.cnn} \
            >{log} 2>&1
        """
# Calculate coverage for anti-target bins
rule antitarget_coverage:
    input:
        bam = "Sample_{sample}/{sample}.bam"
    output:
        cnn = "cnvkit/{sample}/{sample}.antitargetcoverage.cnn"
    conda:
        "../env/cnvkit.yml"
    log:
        "logs/{sample}_cnvkit_antitarget_coverage.log"
    params:
        bed = config['cnvkit']['antitargets']
    shell:
        """
        cnvkit.py coverage \
            {input.bam} \
            {params.bed} \
            --output {output.cnn} \
            >{log} 2>&1
        """
# Combine uncorrected target and antitarget coverage tables and correct for biases 
# according to reference. Output: copy number ratio table (.cnr)
rule fix_coverage:
    input:
        cnn_targets = rules.target_coverage.output.cnn,
        cnn_antitargets = rules.antitarget_coverage.output.cnn
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
            {input.cnn_targets} \
            {input.cnn_antitargets} \
            {params.ref} \
            --output {output.cnr} \
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
        call = "cnvkit/{sample}/{sample}.call.cns",
        vcf = "cnvkit/{sample}/{sample}.cnv.vcf"
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
            --processes {threads} \
            --output {output.call} \
            --vcf {output.vcf} \
            >{log} 2>&1
        """

# Plot scatterplot of all segments
rule cnvkit_scatter:
    input:
        cns = rules.segment.output.cns
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
        pdftoppm -png {input} $(basename {input} .pdf)
        """
