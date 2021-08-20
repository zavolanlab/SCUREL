#!/usr/bin/bash
snakemake \
    --configfile $1 \
    --dag -np | dot -Tsvg > images/dag.svg
EXIT_STATUS_1=$?

snakemake \
    --configfile $1 \
    --dag -np | dot -Tpng > images/dag.png
EXIT_STATUS_2=$?

EXIT EXIT_STATUS_1 || EXIT_STATUS_2
