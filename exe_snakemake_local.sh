#!/bin/bash
mkdir -p ~/tmpdir
export TMPDIR=~/tmpdir
snakemake --rerun-incomplete \
    --configfile $1 \
    --use-conda \
    --conda-prefix $HOME/envs_snakemake \
    --default-resources 'mem_mb=2048' \
    --resources 'mem_mb=30720' \
    --cores 10
