#!/usr/bin/bash
# use rule all 
snakemake --rerun-incomplete \
    --use-conda \
    --configfile $1 \
    --cluster-config "config/cluster.json" \
    --cluster "sbatch \
        --cpus-per-task={threads} \
        --mem={resources.mem_mb}M \
        --time={cluster.time} \
        --qos={cluster.queue} \
        --job-name={cluster.name} \
        -o {cluster.out}" \
    --default-resources 'mem_mb=3072' \
    --cores 24 \
    --local-cores 1
