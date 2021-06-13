#--------------
configfile: "config.yaml"
#--------------
INPUT_DIR = config["fq_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
VOTUS_FA = config["votus_fa"].rstrip("/")
RGI_TAB = config["rgi_tab"].rstrip("/")

BATCH, LIB_ID = glob_wildcards(INPUT_DIR + "/{batch, .*\d+}_{lib_id, .*\d+}")
#---------------

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

rule all:
    input: 
        #expand(OUTPUT_DIR + "/depth/{batch}/{lib_id}.depth", batch = BATCH, lib_id = LIB_ID),
        expand(OUTPUT_DIR + "/depth_analysis/{batch}", batch = BATCH)

rule extract_AMG:
    input:
        votus=VOTUS_FA,
        rgi=RGI_TAB,
    output:
        orf=OUTPUT_DIR + "/AMG/AMG.fasta",
        kept=OUTPUT_DIR +  "/AMG/kept.tsv",
    params:
        prefix = OUTPUT_DIR + "/AMG"
    script: 
        "scripts/extract_RGI_fasta.py"
    conda: 
        "envs/biopython.yaml"
    log:
        OUTPUT_DIR + "/logs/extract_AMG.log"
    shell:
        "python {script} -b {input.rgi} -i {input.votus} -o {prefix} >> {log} 2>&1"

# Choose representative AMG ORF (Clustal Omega)
rule clustalo:
    input:
        OUTPUT_DIR + "/AMG/AMG.fasta"
    output:
        OUTPUT_DIR + "/AMG/AMG.msa.fasta"
    params:
        extra=""
    log:
        OUTPUT_DIR + "/logs/clustalo.log"
    threads: 8
    wrapper:
        "v0.75.0/bio/clustalo"

rule bwt2_build:
    input:
        reference=OUTPUT_DIR + "/AMG/AMG.msa.fasta"
    output:
        multiext(
            OUTPUT_DIR + "/AMG/AMG.msa",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        OUTPUT_DIR + "/logs/bowtie2_build/build.log"
    params:
        extra=""  # optional parameters
    threads: 4
    wrapper:
        "v0.75.0/bio/bowtie2/build"

rule cleanup_fq_names:
    input: 
        index=rules.bwt2_build.output,
        sample=INPUT_DIR + "/{sample}"
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
        sample=[INPUT_DIR + "/{sample}/{sample}_R1.fastq.gz", INPUT_DIR + "/{sample}/{sample}_R2.fastq.gz"]
    output:
        temp(OUTPUT_DIR + "/mapped/{sample}.bam")
    log:
        OUTPUT_DIR + "/logs/bowtie2/{sample}.log"
    params:
        index=OUTPUT_DIR + "/AMG/AMG.msa",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "v0.75.0/bio/bowtie2/align"

rule samtools_sort:
    input:
        OUTPUT_DIR + "/mapped/{sample}.bam"
    output:
        temp(OUTPUT_DIR + "/mapped/{sample}.sort.bam")
    script:
        "scripts/read_count_bam.pl"
    conda:
        "envs/samtools.yaml"
    log:
        OUTPUT_DIR + "/logs/samtools/sorted/{sample}.log"
    shell:
        "(samtools view -h {input} | {script} | samtools view -Su -q30 - | samtools sort -O BAM -o {output} -) 2>> {log}"

rule samtools_index:        