# RGIcat

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.4.1-brightgreen.svg)](https://snakemake.readthedocs.io/)
[![Build Status](https://travis-ci.org/snakemake-workflows/RGIcat.svg?branch=master)](https://travis-ci.org/snakemake-workflows/RGIcat)

A snakemake workflow to build resistence gene catelogue using RGI output from CARD database, and align raw reads against it (dereplicated by cd-hit) to calculate the abundance matrix.

## Rule graph

![rule graph](/rules.png)

## INSTALLATION

This tiny pipeline is tested using snakemake v6.4.1, which shall be compatible with old versions. V6 snakemakw is highly recommended, in that `mamba` substitutes `conda` to speed up the dependency solving.

[Install snakemake and mamba](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

```bash
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

## USAGE

The pipeline follows the Snakemake usage, more details can be found on the [official tutorials](https://snakemake.readthedocs.io).

**To start a new run**, please revise the config file accordingly. The input R1 and R2 fastq files shall be placed in the directories separated by the sample names. Below is a demo for the `config.yaml`.

```
# This file should contain everything to configure the workflow on a global scale.
fq_dir: "/mnt/md0/RGIcat_test/fqs"
results_dir: "/mnt/md0/RGIcat_test/results"
votus_fa: "/mnt/md0/RGIcat_test/vOTUs.fasta"
rgi_tab: "/mnt/md0/RGIcat_test/Predicted_AMG.txt"
threads: 4
```

```
fq_dir: "/path/to/fqs"
results_dir: "/path/to/results"
votus_fa: "/path/to/contig.fasta"
rgi_tab: "/path/to/RGI.tab"
threads: the number of threads
```

**To run the workflow**

```
snakemake --core 6 --use-conda
```