#!/usr/bin/env python


# with curl -sS URL | gunzip > genome.gff

ORGANISM = 'human'

# mouse: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/all_assembly_versions/GCF_000001635.27_GRCm38.p6/
# human: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
if ORGANISM == 'mouse':
    gtffile = 'genome.GRCm38.p6.gff'
    outfile = 'genome.GRCm38.p6.adj.gff'
    COL_CHR = 9
elif ORGANISM == 'human':
    gtffile = 'genome.GRCh38.p13.gff' 
    outfile = 'genome.GRCh38.p13.adj.gff'
    COL_CHR = 9
# download assembly report from same directory as genome annotation
reportfile = 'assembly_report.txt'

# load report
# create dictionary of RefSeq-Accn to Sequence-Name
# use column [9] for mm10 (using chr), and column [9] for hg38
# Only include translated 'chr' chromosome names
# and remove other chromosome names (e.g. unlocalized-scaffold, fix-patch, novel-patch)
# by only including chromosomes without "_"

seqs = {}
with open(reportfile, 'r') as report:
    for line in report.read().splitlines():
        if line[0] != "#":
            entries = line.split("\t")
            if 'chr' in entries[COL_CHR] and '_' not in entries[COL_CHR]:
                seqs[entries[6]] = entries[COL_CHR]

if ORGANISM == 'human':
    # remove chr from name to be consistent with cellranger
    # and rename chrM to MT
    for k, e in seqs.items():
        if e == 'chrM':
            seqs[k] = 'MT'
        else:
            seqs[k] = e.replace('chr', '')


# adjust gtffile 
# 1. remove Gnomon lines
# 2. adjust sequence name from NC_ > chr
with open(gtffile, 'r') as infile, open(outfile, 'w') as o_file:
    for line in infile:
        if line[0] == "#":
            o_file.write(line)
        else:
            entries = line.split("\t")
            # remove Gnomon lines
            # only include lines in chromosome names
            if entries[1] != "Gnomon" and entries[0] in seqs.keys():
                # adjust Sequence-name
                entries[0] = seqs[entries[0]]
                line = '\t'.join(entries)
                o_file.write(line)
        

# sort new gtf (with igvtools)
# index new gtf