#_______ SV CALLING SNIFFLES _________________________________________________________#

rule sv_sniffles: 
  input:
    use_bam
  output:
    vcf="variant_calling/{sample}/{sample}_sniffles.vcf",
    snf="variant_calling/{sample}/{sample}_sniffles.snf"
  conda:
    "../env/sniffles.yml"
  threads:
    12
  log: 
    "logs/{sample}_sniffles.log"
  shell:
    """
    sniffles \
      --input {input} \
      --vcf {output.vcf} \
      --snf {output.snf} \
      --threads {threads} > {log} 2>&1
    """ 

rule process_sniffles_sv:
  input:
    rules.sv_sniffles.output.vcf
  output:
    "Sample_{sample}/{sample}.sv_sniffles.vcf.gz"
  log:
    "logs/{sample}_process_sv_vcf_sniffles.log"
  threads:
    1
  conda:
    "../env/ngs-bits.yml"
  params:
    coregenome="/mnt/projects/external/promethion/reference_genomes/WGS_grch38.bed"
  shell:
    """
    VcfFilter -in {input} -out variant_calling/{wildcards.sample}/{wildcards.sample}.sniffles.filtered.vcf -reg {params.coregenome} > {log} 2>&1
    VcfSort -in variant_calling/{wildcards.sample}/{wildcards.sample}.sniffles.filtered.vcf -out variant_calling/{wildcards.sample}/{wildcards.sample}.sniffles.sorted.vcf >> {log} 2>&1
    bgzip -c variant_calling/{wildcards.sample}/{wildcards.sample}.sniffles.sorted.vcf  > {output} 2>> {log}
    tabix {output}  2>> {log}
    """

#_______ SV CALLING CUTESV _________________________________________________________#

rule sv_cutesv:
  input:
    use_bam
  output:
    vcf="variant_calling/{sample}/{sample}_cutesv.vcf"
  conda:
    "../env/cutesv.yml"
  threads:
    12
  log: 
    "logs/{sample}_cutesv.log"
  params:
    ref = config['ref']['genome'],
    genotyping = "--genotype" if config['sv_cutesv']['enable_genotyping'] else ""
  shell:
    """
    cuteSV \
      --threads {threads} \
      --sample {wildcards.sample} \
      {input} {params.ref} {params.genotyping} {output.vcf} $(dirname {output.vcf}) \
      > {log} 2>&1
    """ 

rule process_cutesv_vcf:
  input:
    rules.sv_cutesv.output.vcf
  output:
    "Sample_{sample}/{sample}.sv_cutesv.vcf.gz"
  log:
    "logs/{sample}_process_sv_vcf_cutesv.log"
  threads:
    1
  conda:
    "../env/ngs-bits.yml"
  params:
    coregenome="/mnt/projects/external/promethion/reference_genomes/WGS_grch38.bed"
  shell:
    """
    VcfFilter -in {input} -out variant_calling/{wildcards.sample}/{wildcards.sample}.cutesv.filtered.vcf -reg {params.coregenome}
    VcfSort -in variant_calling/{wildcards.sample}/{wildcards.sample}.cutesv.filtered.vcf -out variant_calling/{wildcards.sample}/{wildcards.sample}.cutesv.sorted.vcf
    bgzip -c variant_calling/{wildcards.sample}/{wildcards.sample}.cutesv.sorted.vcf  > {output}
    tabix {output}
    """

#_______ SV HAPDUP DIPDIFF _________________________________________________________#

rule map_to_assembly:
  input:
    fq = "Sample_{sample}/{sample}.fastq.gz",
    fa = lambda wildcards: f"Sample_{wildcards.sample}/{wildcards.sample}.asm.{config['assembler']}.fasta"
  output:
    bam = "variant_calling/{sample}_hapdup/assembly_aln.bam"
  log:
    "logs/{sample}_map_assembly.log"
  conda:
    "../env/minimap2.yml"
  threads:
    config['max_threads']
  shell:
     """
     minimap2 -ax map-ont --eqx -t {threads} {input.fa} {input.fq} | samtools sort -@ 4 -m 4G > {output.bam}
     samtools index -@ 4 {output.bam}
     """

rule hapdup:
  input:
    bam = "variant_calling/{sample}_hapdup/assembly_aln.bam",
    fa = lambda wildcards: f"Sample_{wildcards.sample}/{wildcards.sample}.asm.{config['assembler']}.fasta"
  output:
    f1 = "variant_calling/{sample}_hapdup/haplotype_1.fasta",
    f2 = "variant_calling/{sample}_hapdup/haplotype_2.fasta"
  log:
    "logs/{sample}_hapdup.log"
  threads:
   config['max_threads']
  shell:
    """
    docker run \
    -v "$(realpath $(dirname {input.bam}))":"/mnt/hapdup" \
    -v "$(realpath $(dirname {input.fa}))":"/mnt/assembly" \
    -u `id -u`:`id -g` mkolmogo/hapdup:0.6 \
        hapdup \
    --assembly /mnt/assembly/$(basename {input.fa})\
    --bam /mnt/hapdup/$(basename {input.bam}) \
    --out-dir /mnt/hapdup \
    -t {threads} \
    --rtype ont \
    >{log} 2>&1
   """

rule dipdiff:
  input:
    ref = config['ref']['genome'],
    f1 = "variant_calling/{sample}_hapdup/haplotype_1.fasta",
    f2 = "variant_calling/{sample}_hapdup/haplotype_2.fasta"
  output:
    vcf = "variant_calling/{sample}_dipdiff/variants.vcf"
  log:
    "logs/{sample}_dipdiff.log"
  threads:
    config['max_threads']
  shell: 
    """
    docker run \
        -v "$(realpath $(dirname {input.ref}))":"/mnt/ref" \
        -v "$(realpath $(dirname {input.f1}))":"/mnt/hapdup" \
        -v "$(realpath $(dirname {output.vcf}))":"/mnt/output" \
        -u `id -u`:`id -g` mkolmogo/dipdiff:0.3 \
        dipdiff.py \
        --reference "/mnt/ref/$(basename {input.ref})" \
        --pat "/mnt/hapdup/haplotype_1.fasta" \
        --mat "/mnt/hapdup/haplotype_2.fasta" \
        --out-dir /mnt/output \
        -t {threads} \
        >{log} 2>&1
    """

rule process_dipdiff_vcf:
  input:
    vcf = "variant_calling/{sample}_dipdiff/variants.vcf"
  output:
    filtered = "variant_calling/{sample}_dipdiff/{sample}.filtered.vcf",
    sorted = "variant_calling/{sample}_dipdiff/{sample}.sorted.vcf",
    final = "Sample_{sample}/{sample}.sv_dipdiff.vcf.gz"
  threads:
    1
  conda:
    "../env/ngs-bits.yml"
  params:
    target_region = f"-reg {config['ref']['target_region']}" if config['ref']['target_region'] else ""
  shell:
    """
    VcfFilter -in {input} -out {output.filtered} {params.target_region}
    VcfSort -in {output.filtered} -out {output.sorted}
    bgzip -c {output.sorted}  > {output.final}
    tabix {output.final}
    """