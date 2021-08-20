#!/usr/bin/bash
set -o pipefail

snakemake \
    --configfile $1 \
    --rulegraph -np | dot -Tsvg > images/rule_graph.svg

snakemake \
    --configfile $1 \
    --rulegraph -np | dot -Tpng > images/rule_graph.png
