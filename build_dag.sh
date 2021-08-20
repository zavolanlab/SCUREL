snakemake \
    --configfile $1 \
    --dag -np | dot -Tsvg > images/dag.svg

snakemake \
    --configfile $1 \
    --dag -np | dot -Tpng > images/dag.png
