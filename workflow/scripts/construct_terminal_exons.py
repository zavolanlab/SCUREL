#!/usr/bin/env python
"""
Script for obtaining terminal exons from a genome annotation file (GFF)
"""

import logging
import argparse
import os

import numpy as np
import pandas as pd

from intervaltree import Interval, IntervalTree

__author__ = "Dominik Burri"
__date__ = "2021-03-30"
__license__ = "Apache"
__version__ = "0.0.1"
__maintainer__ = "Dominik Burri"
__email__ = "dominik.burri@unibas.ch"
__status__ = "Development"

def downstream_exon(current, challenge):
    '''
    Return for two exons the one more downstream (3p direction)
    current and challenge are each the list of entries of the gff line
    '''
    # set current exon as downstream
    downstream = current
    # sanity check: require same chromosome and strand
    if current[0] != challenge[0] or current[6] != challenge[6]:
        raise Exception("Problem with exon definitions.\ncurrent=%s\nchallenge=%s" %(str(current), str(challenge)))
    # plus strand: higher position relevant
    # minus strand: lower position relevant
    if current[6] == '+':
        if current[4] < challenge[4]:
            downstream = challenge
    elif current[6] == '-':
        if current[3] > challenge[3]:
            downstream = challenge
    else:
        # TODO: raise error, no valid strand provided
        pass
    return(downstream)

def non_overlapping_exons(df_exon):
    '''
    Return df of non-overlapping exons
    df_exon: df of exons in gff format, all on same strand and chromosome
    '''
    df_exon['length'] = df_exon['end'] - df_exon['start']
    # create intervaltree of exons
    tree = IntervalTree()
    for i, row in df_exon.iterrows():
        tree[row.start:row.end] = (row.name, row.length)
    # merge overlapping regions and keep the longest exon
    tree.merge_overlaps(data_reducer=lambda current, nex: current if current[1] > nex[1] else nex)
    res_df = df_exon.loc[[i.data[0] for i in tree]]
    return(res_df)

def get_terminal_exons(genome_file, type_id, gene_type, transcript_source):
    '''
    Reads the genome file (gff v3) and extracts terminal exons
    returns pd.DataFrame with terminal exons
    '''
    gene_ids = [] # created by genes
    genes = {} # created by transcripts
    transcripts = {} # created by exons
    # construct gene, transcript and exon relationship as nested dictionary
    with open(genome_file, 'r') as infile:
        for line in infile:
            if line[0] != "#":
                entries = line.split("\t")
                attributes = {attr.split("=")[0]:attr.split("=")[1].rstrip() for attr in entries[8].split(";")}
                if entries[2] == 'gene':
                    # Gather info about selected gene biotype
                    if attributes[type_id] in gene_type:
                            gene_ids.append(attributes['ID'])
                elif entries[2] in transcript_source:
                    # add gene to dictionary
                    if attributes['Parent'] not in genes:
                        genes[attributes['Parent']] = [attributes['ID']]
                    else:
                        # append transcript ID to gene
                        genes[attributes['Parent']].append(attributes['ID'])
                # consider type exon (column 3)
                elif entries[2] == 'exon':
                    # add parent transcript to dict if not available
                    if attributes['Parent'] not in transcripts:
                        # add transcript if not available
                         transcripts[attributes['Parent']] = entries
                    else:
                        # add exon to transcript if exon more downstream than current exon
                        transcripts[attributes['Parent']] = downstream_exon(transcripts[attributes['Parent']], entries)
    # Create dataframe of terminal exons
    df = pd.DataFrame()
    for g in gene_ids:
        if g in genes:
            for t in genes[g]:
                dfe = pd.DataFrame.from_dict({t:[g]+transcripts[t]}, orient = 'index')
                df = pd.concat([df, dfe])
    # define column names: 0 gene id, 1:9 from gff specification
    df.columns = ['gene', 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    df = df.astype({'start': int, 'end': int})
    return(df)

#### main program
def main():
    parser = argparse.ArgumentParser(description="Obtain terminal exons (TE) for genome.")
    parser.add_argument('--type-id', dest='type_id', 
      required=False, help='Specify for which attribute to filter. Default (refseq): gene_biotype',
      default='gene_biotype')
    parser.add_argument('--gene-type', dest='gene_type', 
      required=False, help='Specify for which attribute value to filter. Default (refseq): protein_coding lncRNA',
      default='protein_coding lncRNA', nargs='*')
    parser.add_argument('--source', dest='transcript_source', 
      required=False, help='Specify for which transcript type (column 3) to filter. Default: mRNA lnc_RNA',
      default='mRNA lnc_RNA', nargs='*')
    parser.add_argument('--log-file', dest="log_file", 
      required=False, help='Logging file. Default: construct_terminal_exons.log.', default='construct_terminal_exons.log')
    parser.add_argument('GFF_IN',
      help='Input genome annotation file in gff v3.')
    parser.add_argument('BED_OUT',
      help='Output file in BED format.')
    args = parser.parse_args()

    # logging settings: from snakemake
    logging.basicConfig(filename=args.log_file, level=logging.INFO, 
                    format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')
    
    # read genome annotation file (GFF v3)
    # No sorting is assumed
    genome_file = args.GFF_IN
    logging.info("Will use GFF %s" % genome_file)
    out_file = args.BED_OUT
    logging.info("Will write output to %s" % out_file)

    # filter entries by gene products (keywords might change)
    type_id = args.type_id
    gene_type = args.gene_type
    transcript_source = args.transcript_source # column 3
    logging.info("Parameters: type_id: %s, gene_type: %s, source: %s" %(type_id, gene_type, transcript_source))

    # Obtain and filter terminal exons
    df = get_terminal_exons(genome_file, type_id, gene_type, transcript_source)
    logging.info("Obtained %d terminal exons" % df.shape[0])
    # Obtain non-overlapping TEs per strand and gene
    chroms = list(set(df.seqid))
    chroms.sort()
    res = pd.DataFrame()
    for chrom in chroms:
        df_chr_plus = df.loc[(df.seqid == chrom) & (df.strand == "+")].iloc[:,1:]
        res_plus = non_overlapping_exons(df_chr_plus)
        res = res.append(res_plus)
        df_chr_minus = df.loc[(df.seqid == chrom) & (df.strand == "-")].iloc[:,1:]
        res_minus = non_overlapping_exons(df_chr_minus)
        res = res.append(res_minus)
    logging.info("Retained %d terminal exons after removing overlaps" % res.shape[0])
    # sort dataframe by chromosome and start position
    res.sort_values(by=['seqid', 'start'], inplace = True)
    # add transcript id as column (might depend on genome annotation)
    res['id'] = [x.split('-')[1].split('.')[0] for x in res.index]
    # save Terminal exons to BED
    res.to_csv(out_file, sep = '\t', header = False, index = False,
        columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    logging.info("Wrote results to: " + out_file)
    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()