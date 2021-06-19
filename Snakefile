#--------------
configfile: "config.yaml"
#--------------
INPUT_DIR=config["fq_dir"].rstrip("/")
OUTPUT_DIR=config["results_dir"].rstrip("/")
VOTUS_FA=config["votus_fa"].rstrip("/")
RGI_TAB=config["rgi_tab"].rstrip("/")
THREADS=config["threads"]

SAMPLE=glob_wildcards(INPUT_DIR + "/{dir}/{file}.fastq.gz").dir
#---------------

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"
rule all:
    input:
        OUTPUT_DIR + "/count_matrix.tsv",

rule extract_AMG:
    input:
        votus=VOTUS_FA,
        rgi=RGI_TAB,
    output:
        orf=OUTPUT_DIR + "/AMG/AMG.fasta",
        kept=OUTPUT_DIR +  "/AMG/kept.tsv",
    params:
        OUTPUT_DIR + "/AMG"
    conda:
        "envs/biopython.yaml"
    log:
        OUTPUT_DIR + "/logs/extract_AMG.log"
    shell:
        """
        scripts/extract_RGI_fasta.py -b {input.rgi} -i {input.votus} -o {params} >> {log} 2>&1
        """

# Choose representative AMG ORF
# Clustal Omega discarded (only do alignment), use cd-hit instead
# 95% similarity in 90% covered alignment, > 100 bp
rule cd_hit:
    input:
        OUTPUT_DIR + "/AMG/AMG.fasta"
    output:
        OUTPUT_DIR + "/AMG/AMG.cdhit.fasta"
    conda:
        "envs/cdhit.yaml"
    log:
        OUTPUT_DIR + "/logs/cdhit.log"
    threads: THREADS
    shell:
        "cd-hit-est -i {input} -o {output} -c 0.95 -n 8 -l 100 -aS 0.9 -d 0 -B 0 -T 2 -M 10000"

rule bwt2_build:
    input:
        reference=OUTPUT_DIR + "/AMG/AMG.cdhit.fasta"
    output:
        multiext(
            OUTPUT_DIR + "/AMG/AMG.cdhit",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        OUTPUT_DIR + "/logs/bowtie2_build/build.log"
    conda:
        "envs/bwt2sam.yaml" # libtbb.so.2 incompetence due to conda-forge updates
    params:
        extra=""  # optional parameters
    threads: THREADS
    wrapper:
        "v0.75.0/bio/bowtie2/build"

rule cleanup_fq_names:
    input: 
        sample=INPUT_DIR + "/{sample}",
    output:
        r1 = INPUT_DIR + "/{sample}/{sample}_R1.fastq.gz",
        r2 = INPUT_DIR + "/{sample}/{sample}_R2.fastq.gz"
    shell:
        """
        mv {input.sample}/*R1*.fastq.gz {output.r1}
        mv {input.sample}/*R2*.fastq.gz {output.r2}
        """

rule bwt2_map:
    input:
        index=rules.bwt2_build.output,
        sample=[INPUT_DIR + "/{sample}/{sample}_R1.fastq.gz", INPUT_DIR + "/{sample}/{sample}_R2.fastq.gz"]
    output:
        temp(OUTPUT_DIR + "/mapped/{sample}.bam")
    log:
        OUTPUT_DIR + "/logs/bowtie2/{sample}.log"
    conda:
        "envs/bwt2sam.yaml" # libtbb.so.2 incompetence due to conda-forge updates
    params:
        index=OUTPUT_DIR + "/AMG/AMG.cdhit",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: THREADS  # Use at least two threads
    wrapper:
        "v0.75.0/bio/bowtie2/align"

rule samtools_sort:
    input:
        OUTPUT_DIR + "/mapped/{sample}.bam"
    output:
        temp(OUTPUT_DIR + "/mapped/{sample}.sorted.bam")
    conda:
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/logs/samtools/sorted/{sample}.log"
    shell:
        "(samtools view -h {input} | scripts/read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o {output} -) 2>> {log}"

rule samtools_index:        
    input:
        OUTPUT_DIR + "/mapped/{sample}.sorted.bam"
    output:
        temp(OUTPUT_DIR + "/mapped/{sample}.sorted.bam.bai")
    log:
        OUTPUT_DIR + "/logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        THREADS     # This value - 1 will be sent to -@
    wrapper:
        "v0.75.0/bio/samtools/index"

rule header_sample:
    input:
        expand(OUTPUT_DIR + "/mapped/{sample}.sorted.bam.bai", sample=SAMPLE)
    output:
        temp(OUTPUT_DIR + "/header_sample")
    params:
        sample=SAMPLE
    log:
        OUTPUT_DIR + "/logs/header_sample.log"
    run:
        """
        with open(OUTPUT_DIR + "/header_sample", 'w') as f:
             for i in params.sample:
                 f.write("%s\t" % i)
        """

rule gene_names:
    input:
        expand(OUTPUT_DIR + "/mapped/{sample}.sorted.bam.bai", sample=SAMPLE)
    output:
        temp(OUTPUT_DIR + "/gene_names")
    conda:
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/logs/gene_names.log"
    shell:
        """
        samtools idxstats $(ls {input} | head -n 1) | grep -v "*" | cut -f1 > gene_names
        """
       
rule gene_count:
    input:
        bam=OUTPUT_DIR + "/mapped/{sample}.sorted.bam",
        bai=OUTPUT_DIR + "/mapped/{sample}.sorted.bam.bai"
    output:
        temp(OUTPUT_DIR + "/counts/{sample}.count")
    conda:
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/logs/counts/{sample}.log"
    shell:
        """
        samtools idxstats {input.bam} | grep -v "*" | cut -f3 > {output}
        """

rule count_matrix:
    input:
        gene_count=expand(OUTPUT_DIR + "/counts/{sample}.count", sample=SAMPLE),
        header_sample=OUTPUT_DIR + "/header_sample",
        gene_name=OUTPUT_DIR + "/gene_names",
    output:
        OUTPUT_DIR + "/count_matrix.tsv"
    log:
        OUTPUT_DIR + "/logs/count_matrix.log"
    shell:
        """
        paste {input.gene_name} {input.gene_count} | cat {input.header_sample} - > {output}
        """