"""
Script for randomising sample labels from 3' end coverage in BED intervals
Multithreading: Split BED intervals into 100 chunks for one process. 
For files between 500 - 1000 MB, each core needs approx. 5GB RAM.
"""

import logging
import argparse
import os
import random
import sys
from itertools import repeat
import multiprocessing

import numpy as np
import pandas as pd

__author__ = "Dominik Burri"
__date__ = "2020-11-25"
__license__ = "Apache"
__version__ = "0.0.1"
__maintainer__ = "Dominik Burri"
__email__ = "dominik.burri@unibas.ch"
__status__ = "Development"

def load_sample(file):
    'Load sample with cumulative coverage'
    sample = pd.read_csv(file,
                sep="\t", header=None, 
                names = ["chrom", "start", "end", "coverage"],
                dtype={0: str, 1: np.int64, 2: np.int64, 3: np.int64})
    logging.info("Loaded: " + file)
    return sample

def write_sample(df, file):
    'Write df to file'
    df.to_csv(file, sep = "\t", index = False, header = False, columns = ["chrom", "start", "end", "coverage"])
    logging.info("Wrote: " + file)

def get_coverage(te, sample):
    'Obtain 3p end coverage for sample in correct strand'
    tmp = sample[sample.chrom == te.chr]
    tmp = tmp[tmp.start >= te.TE_start]
    tmp = tmp[tmp.end <= te.TE_end]
    return(tmp)

def unique(list1): 
    x = np.array(list1) 
    return(np.unique(x).tolist()) 

def fill_coverage(rnd):
    'Create new pandas Dataframe from list of end positions'
    rows = unique(rnd)
    res = pd.DataFrame({'start': [x - 1 for x in rows], 
        'end': rows, 
        'coverage': [rnd.count(x) for x in rows]})
    return res

def extend_coverage(TE, df):
    'Extend dataframe to include zero coverage and TE borders'
    df.reset_index(drop = True, inplace = True)
    # include zero coverage
    for i in range(len(df.index)-1):
        if df.loc[i+1].start - df.loc[i].end > 0:
            df = df.append({'chrom': TE.chr, 'start': df.loc[i].end, 'end': df.loc[i+1].start, 'coverage': 0}, ignore_index = True)
    df.sort_values(by = 'start', inplace = True)
    # extend df region to include TE start and end positions
    if df.loc[0].start > TE.TE_start:
        tmp = pd.DataFrame({'chrom': TE.chr, 'start': TE.TE_start, 'end': df.loc[0].start, 'coverage': 0}, 
            index = [0])
        df = pd.concat([tmp, df], ignore_index = True)
    if df.loc[df.index[-1]].end < TE.TE_end:
        tmp = pd.DataFrame({'chrom': TE.chr, 'start': df.loc[df.index[-1]].end, 'end': TE.TE_end, 'coverage': 0}, 
            index = [df.index[-1]+1])
        df = pd.concat([df, tmp])
    return df

def get_randomised_reads(cov1, cov2, n_reads):
    'Randomise labels in coverage dataframes cov1 and cov2, with n_reads reads for sample1.'
    # obtain chromosome
    chrom = set(cov1.chrom).pop()
    # create vector of bp positions, one for each read
    vec = np.append(np.repeat(np.array(cov1.end), cov1.coverage), 
        np.repeat(np.array(cov2.end), cov2.coverage)).tolist()
    # sample reads from vector
    rnd = random.sample(vec, n_reads)
    # sort base positions
    rnd.sort()
    # prepare new df
    tmp1 = fill_coverage(rnd)
    tmp1.insert(loc = 0, column = 'chrom', value = chrom)
    # get coverage for sample2 by removing elements from vec
    for entry in rnd:
        vec.remove(entry)
    vec.sort()
    tmp2 = fill_coverage(vec)
    tmp2.insert(loc = 0, column = 'chrom', value = chrom)
    return([tmp1, tmp2])

def compute_random_labels(zip_tuple):
    global sample1
    global sample2
    te = zip_tuple[1]
    cov1 = get_coverage(te, sample1)
    cov2 = get_coverage(te, sample2)
    # obtain nr. reads per sample
    n1 = cov1.coverage.sum()
    n2 = cov2.coverage.sum()
    if n1 == 0 and n2 == 0:
        return [cov1, cov2]
    elif n1 == 0:
        return [cov1, extend_coverage(te, cov2)]
    elif n2 == 0:
        return [extend_coverage(te, cov1), cov2]
    if n1+n2 > sys.maxsize:
        logging.warn("Sum of coverage would exceed maximum size of python list! \
            Will reduce coverage by 2 until the sum is below the maximum. \
            TE affected: %s" % te.TE_id)
        while(n1+n2 > sys.maxsize):
            # Integer ceiling division
            cov1.coverage = -(cov1.coverage // -2)
            cov2.coverage = -(cov2.coverage // -2)
            n1 = cov1.coverage.sum()
            n2 = cov2.coverage.sum()
    # Main logic: randomise reads in the two samples
    tmp1, tmp2 = get_randomised_reads(cov1, cov2, n1)
    # add zero coverage to capture whole TE
    tmp1 = extend_coverage(te, tmp1)
    tmp2 = extend_coverage(te, tmp2)
    return [tmp1, tmp2]

def pool_starmap(TEs, num_cores):
    out1 = pd.DataFrame({'chrom': pd.Series([], dtype = 'str'),
        'start': pd.Series([], dtype = np.int64),
        'end': pd.Series([], dtype = np.int64),
        'coverage': pd.Series([], dtype = np.int64)})
    out2 = pd.DataFrame({'chrom': pd.Series([], dtype = 'str'),
        'start': pd.Series([], dtype = np.int64),
        'end': pd.Series([], dtype = np.int64),
        'coverage': pd.Series([], dtype = np.int64)})
    try:
        pool = multiprocessing.Pool(num_cores)
        logging.info("Started multiprocessing.Pool.imap with %s cores." % num_cores)
        out = pool.imap(compute_random_labels, 
            TEs.iterrows(),
            chunksize = 100)
        for entry in out:
            out1 = out1.append(entry[0], ignore_index = True)
            out2 = out2.append(entry[1], ignore_index = True)
    except Exception as e:
        logging.error("Error occurred during multiprocessing: " + str(e))
        pool.terminate()
        sys.exit(1)
    finally:
        pool.close()
        pool.join() 
        logging.info("multiprocessing.Pool closed.")
    return((out1, out2))

#### main program
def main():
    parser = argparse.ArgumentParser(description="Randomise 3p end positions in BED intervals (stranded) for two samples (stranded).")
    parser.add_argument('-bed', '--bed-file', dest='bed_file', 
      required=True, help='Input BED file')
    parser.add_argument('-s1p', '--one-plus', dest='sample1_plus', 
      required=True, help='Sample 1 3p end coverage bedgraph plus strand')
    parser.add_argument('-s1m', '--one-minus', dest='sample1_minus', 
      required=True, help='Sample 1 3p end coverage bedgraph minus strand')
    parser.add_argument('-s2p', '--two-plus', dest='sample2_plus', 
      required=True, help='Sample 2 3p end coverage bedgraph plus strand')
    parser.add_argument('-s2m', '--two-minus', dest='sample2_minus', 
      required=True, help='Sample 2 3p end coverage bedgraph minus strand')
    parser.add_argument('-o1p', dest='out1_plus', 
      required=False, default="randomised1.plus.bedgraph", help='Output file for sample1 plus strand')
    parser.add_argument('-o1m', dest='out1_minus', 
      required=False, default="randomised1.minus.bedgraph", help='Output file for sample1 minus strand')
    parser.add_argument('-o2p', dest='out2_plus', 
      required=False, default="randomised2.plus.bedgraph", help='Output file for sample2 plus strand')
    parser.add_argument('-o2m', dest='out2_minus', 
      required=False, default="randomised2.minus.bedgraph", help='Output file for sample2 minus strand')
    parser.add_argument('-c', '--cores', dest='NUM_CORES',
      required=False, default=1, type=int,
      help='Specify number of cores for parallelisation. Default: 1.')
    parser.add_argument('-log', '--log-file', dest="log_file", 
      required=False, help='Logging file. Default: label_randomisation.log.', default='label_randomisation.log')
    args = parser.parse_args()
    # logging settings: from snakemake
    logging.basicConfig(filename=args.log_file, level=logging.DEBUG, 
                    format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')
    TEs = pd.read_csv(args.bed_file,
            sep="\t", header=None, 
            names=["chr", "TE_start", "TE_end", "TE_id", "score", "strand"],
            dtype={"chr": str, "TE_start": np.int64, "TE_end": np.int64, "TE_id": str, "score": str, "strand": str})
    
    global sample1 
    global sample2
    
    sample1_plus = load_sample(args.sample1_plus)
    sample2_plus = load_sample(args.sample2_plus)
    sample1 = sample1_plus
    sample2 = sample2_plus
    # obtain coverage in TE region for plus strand
    TEs_plus = TEs[TEs.strand == '+']
    out1_plus, out2_plus = pool_starmap(TEs_plus, args.NUM_CORES)
    logging.info("Finished processing %d TEs on plus strand" % (TEs_plus.shape[0]))
    # write to file
    write_sample(out1_plus, args.out1_plus)
    write_sample(out2_plus, args.out2_plus)
    del out1_plus
    del out2_plus
    # for minus strand
    sample1_minus = load_sample(args.sample1_minus)
    sample2_minus = load_sample(args.sample2_minus)
    sample1 = sample1_minus
    sample2 = sample2_minus
    TEs_minus = TEs[TEs.strand == '-']
    out1_minus, out2_minus = pool_starmap(TEs_minus, args.NUM_CORES)
    logging.info("Finished processing %d TEs on minus strand" % (TEs_minus.shape[0]))
    write_sample(out1_minus, args.out1_minus)
    write_sample(out2_minus, args.out2_minus)
    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()