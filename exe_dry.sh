#!/usr/bin/bash
# use rule all 
snakemake \
    --use-conda \
    --configfile $1 \
    --dryrun \
    --printshellcmds \
    --cores 12 
