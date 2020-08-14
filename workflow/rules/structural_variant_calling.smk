#_______ SV CALLING SNIFFLES _________________________________________________________#

rule sv_sniffles: 
  input:
    rules.map_genome_all.output.bam
  output:
    vcf="variant_calling/{sample}/{sample}_sniffles.vcf"
  conda:
    "../env/sniffles.yml"
  threads:
    12
  log: 
    "logs/{sample}_sniffles.log"
  shell:
    """
    sniffles -m {input} -v {output.vcf} -t {threads} > {log} 2>&1
    """ 

rule process_sv_vcf:
  input:
    rules.sv_sniffles.output.vcf
  output:
    "{sample}/{sample}.sv_sniffles.vcf.gz"
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