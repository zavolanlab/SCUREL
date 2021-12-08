#!/usr/bin/env python
"""
Script for obtaining cumulative coverage per annotated TE
The coverage tracks need to include zero values to capture all base pairs
"""

import logging
import argparse
import os
import multiprocessing

import numpy as np
import pandas as pd
from itertools import repeat

import intervaltree

__author__ = "Dominik Burri"
__date__ = "2019-12-18"
__license__ = "Apache"
__version__ = "0.0.1"
__maintainer__ = "Dominik Burri"
__email__ = "dominik.burri@unibas.ch"
__status__ = "Development"

def get_coverage(te, sample):
    'Obtain 3p end coverage for sample in correct strand'
    tmp = sample[sample.chrom == te.chr]
    tmp = tmp[tmp.start >= te.TE_start]
    tmp = tmp[tmp.end <= te.TE_end]
    return(tmp)

def cumulative_coverage(TE, ivt_samples):
    ''' compute cumulative coverage in TE regions '''
    # manually intersect
    # return empty dataframe if no coverage in TE chromosome 
    if TE.chr not in ivt_samples:
        return pd.DataFrame()
    # find the regions completely in the frame
    df = ivt_samples[TE.chr].overlap(TE.TE_start, TE.TE_end)
    # check if any coverage is found
    if len(df) == 0: return(df)
    # next, find the regions just adjacent to the TE region and include it into df
    df_sorted = sorted(df)
    # create regions encompasing the TE region
    if (df_sorted[0].begin > TE.TE_start):
        df.add(intervaltree.Interval(TE.TE_start, df_sorted[0].begin, 0))
    if (df_sorted[-1].end < TE.TE_end):
        df.add(intervaltree.Interval(df_sorted[-1].end, TE.TE_end, 0))
    # cast to pandas df
    df = pd.DataFrame(sorted(df))
    df.columns = ["start", "end", "coverage"]
    df.insert(0, "chrom", TE.chr)
    # compute cumulative sum
    if TE.strand == '+':
        df.insert(4, "cumsum", np.flip(np.cumsum(np.flip(df.coverage))))
    elif TE.strand == '-':
        df.insert(4, "cumsum", np.cumsum(df.coverage))
    else:
        df.insert(4, "cumsum", -1)
    return df

def compute_cumsum_transcripts(ivt_samples, strand, TEs):
    ''' compute and return df of cumsum per transcript '''
    res = pd.DataFrame({'chrom': pd.Series([], dtype = 'str'),
        'start': pd.Series([], dtype = np.int64),
        'end': pd.Series([], dtype = np.int64),
        'coverage': pd.Series([], dtype = np.int64),
        'cumsum': pd.Series([], dtype = np.int64)})
    for i in TEs.index:
        TE = TEs.loc[i]
        # continue, if TE not on the same strand as the df coverage as then cumulative coverage not of interest
        if TE.strand != strand:
            logging.debug("Transcript %s on opposite strand %s." %(TE.TE_id, TE.strand))
            continue
        df = cumulative_coverage(TE, ivt_samples)
        # write to file directly
        if len(df) == 0:
            logging.debug("No coverage for transcript %s." % TE.TE_id)
        else:
            res = res.append(df, ignore_index = True)
            logging.debug("Added cumsum of transcript %s with %d lines." %(TE.TE_id, len(df.index)))
    return res

def interval_tree(df):
    'Function to create interval tree for given df'
    tree = intervaltree.IntervalTree()
    for index, row in df.iterrows():
        tree[row.start:row.end] = row.coverage
    return(tree)

def interval_tree_from_bedgraph(df_sample, n_cores):
    ''' Create for each chromosome interval tree containing coverages '''
    trees = {}
    chroms = list(set(df_sample.chrom))
    # slice df_sample into chroms
    dfs = []
    for chrom in chroms:
        df = df_sample[df_sample.chrom == chrom]
        dfs.append(df)
    try:
        pool = multiprocessing.Pool(n_cores)
        logging.info("Started multiprocessing pool with %s cores." % n_cores)
        res = pool.imap(interval_tree, dfs, chunksize = 6)
        for chrom, entry in zip(chroms, res):
            trees[chrom] = entry
    except Exception as e:
        logging.error("Error occurred during multiprocessing: " + str(e))
        pool.terminate()
        sys.exit(1)
    finally:
        pool.close()
        pool.join() 
        logging.info("multiprocessing.Pool closed.")
    return(trees)

#### main program
def main():
    parser = argparse.ArgumentParser(description="Compute reverse cumulative distribution for bedgraph tracks in intervals of BED file")
    parser.add_argument('-bed', '--bed-file', dest='bed_file', 
      required=True, help='Input BED file')
    parser.add_argument('-i', '--coverage', dest='coverage', 
      required=True, help='Specify input coverage bedgraph.')
    parser.add_argument('-o', '--output', dest='output', 
      required=True, help='Specify output cumulative coverage bedgraph.')
    parser.add_argument('-s', '--strand', dest='strand',
      required=True, type=str, help='Specify genome orientation of coverage tracks, either plus or minus.')
    parser.add_argument('-c', '--cores', dest='NUM_CORES',
      required=False, default=1, type=int,
      help='Specify number of cores for parallelisation. Default: 1.')
    parser.add_argument('-log', '--log-file', dest="log_file", 
      required=False, help='Logging file', default='cumulative_coverage.log')
    args = parser.parse_args()

    # logging settings
    logging.basicConfig(filename=args.log_file, level=logging.INFO, 
                    format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')

    if not os.path.isfile(args.bed_file):
        logging.error("BED file does not exist")
        exit(-1)

    if not os.path.isfile(args.coverage):
        logging.error("bedgraph file does not exist")
        exit(-1)

    if args.strand == "plus":
        strand = "+"
    elif args.strand == "minus":
        strand = "-"

    # load Terminal Exon annotations
    TEs = pd.read_csv(args.bed_file,
        sep="\t", header=None, 
        names=["chr", "TE_start", "TE_end", "TE_id", "score", "strand"],
        dtype={"chr": str, "TE_start": np.int64, "TE_end": np.int64, "TE_id": str, "score": str, "strand": str})
    logging.info("Loaded: " + args.bed_file)
    # include only TEs with proper strand
    TEs_strand = TEs[TEs.strand == strand]
    logging.info("Will use %d TEs on %s strand." %(len(TEs_strand), strand))

    # load sample coverages
    samples = pd.read_csv(args.coverage,
                sep="\t", header=None, 
                names = ["chrom", "start", "end", "coverage"],
                dtype={0: str, 1: np.int64, 2: np.int64, 3: np.int64})
    logging.info("Loaded: " + args.coverage)

    # filter sample coverage by TE region
    logging.info("Start filter sample coverage by BED regions.")
    try:
        pool = multiprocessing.Pool(args.NUM_CORES)
        logging.info("Started multiprocessing pool with %s cores." % args.NUM_CORES)
        res = pool.starmap(get_coverage,
            zip([TEs.iloc[i] for i in TEs_strand.index], repeat(samples)), 
            chunksize = 50)
        samples = pd.concat(res, ignore_index = True)
    except Exception as e:
        logging.error("Error occurred during multiprocessing: " + str(e))
        pool.terminate()
        sys.exit(1)
    finally:
        pool.close()
        pool.join() 
        logging.info("multiprocessing.Pool closed.")
    logging.info("Finished filter sample coverage by BED regions.")

    # use interval tree
    logging.info("Start computing intervaltree for: " + args.coverage)
    samples = interval_tree_from_bedgraph(samples, args.NUM_CORES)
    logging.info("Finished computing intervaltree for: " + args.coverage)

    # Compute cumsum in TE regions
    logging.info("Start with computation of cumulative sum.")
    res = compute_cumsum_transcripts(samples, strand, TEs_strand)
    logging.info("Finished with computation of cumulative sum.")
    res.to_csv(args.output, sep = "\t", index = False, 
        header = False, 
        columns = ["chrom", "start", "end", "cumsum"])
    logging.info("Saved: " + args.output)
    logging.info("Finished script.")

if __name__ == "__main__":
    main()