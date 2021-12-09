#!/bin/bash
mkdir -p ~/tmpdir
export TMPDIR=~/tmpdir
snakemake --rerun-incomplete \
    --configfile $1 \
    --use-conda \
    --default-resources 'mem_mb=2048' \
    --resources 'mem_mb=8192' \
    --cores 4
