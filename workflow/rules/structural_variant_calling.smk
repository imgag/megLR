#_______ SV CALLING SNIFFLES _________________________________________________________#

rule sv_sniffles: 
  input:
    rules.map_genome_all.output.bam
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

rule process_sv_vcf:
  input:
    rules.sv_sniffles.output.vcf
  output:
    "Sample_{sample}/{sample}.sv_sniffles.vcf.gz"
  threads:
    1
  conda:
    "../env/ngs-bits.yml"
  params:
    coregenome="/mnt/projects/external/promethion/reference_genomes/WGS_grch38.bed"
  shell:
    """
    VcfFilter -in {input} -out variant_calling/{wildcards.sample}/{wildcards.sample}.filtered.vcf -reg {params.coregenome}
    VcfSort -in variant_calling/{wildcards.sample}/{wildcards.sample}.filtered.vcf -out variant_calling/{wildcards.sample}/{wildcards.sample}.sorted.vcf
    bgzip -c variant_calling/{wildcards.sample}/{wildcards.sample}.sorted.vcf  > {output}
    tabix {output}
    """

#_______ SV CALLING CUTESV _________________________________________________________#

rule sv_cutesv:
  input:
    rules.map_genome_all.output.bam
  output:
    vcf="variant_calling/{sample}/{sample}_cutesv.vcf"
  conda:
    "../env/cutesv.yml"
  threads:
    12
  log: 
    "logs/{sample}_cutesv.log"
  params:
    ref = config['ref']['genome']
  shell:
    """
    cuteSV \
      --threads {threads} \
      --sample {wildcards.sample} \
      {input} {params.ref} {output.vcf} $(dirname {output.vcf}) \
      > {log} 2>&1
    """ 

rule process_cutesv_vcf:
  input:
    rules.sv_cutesv.output.vcf
  output:
    "Sample_{sample}/{sample}.sv_cutesv.vcf.gz"
  threads:
    1
  conda:
    "../env/ngs-bits.yml"
  params:
    coregenome="/mnt/projects/external/promethion/reference_genomes/WGS_grch38.bed"
  shell:
    """
    VcfFilter -in {input} -out variant_calling/{wildcards.sample}/{wildcards.sample}.filtered.vcf -reg {params.coregenome}
    VcfSort -in variant_calling/{wildcards.sample}/{wildcards.sample}.filtered.vcf -out variant_calling/{wildcards.sample}/{wildcards.sample}.sorted.vcf
    bgzip -c variant_calling/{wildcards.sample}/{wildcards.sample}.sorted.vcf  > {output}
    tabix {output}
    """