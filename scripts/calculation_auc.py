#!/usr/bin/env python
"""
Script for obtaining AUC per annotated TE
Between sample1 and sample2
"""

import logging
import argparse
import os

import numpy as np
import pandas as pd
import pickle

__author__ = "Dominik Burri"
__date__ = "2020-11-19"
__license__ = "Apache"
__version__ = "0.0.1"
__maintainer__ = "Dominik Burri"
__email__ = "dominik.burri@unibas.ch"
__status__ = "Development"

def compute_area(x,y):
    "Compute Area under the curve by integration"
    return np.trapz(y,x)

def compute_cpm(counts):
    "Compute counts per million for a vector of counts"
    return counts / np.sum(counts) * 1e6

def load_sample(file):
    'Load sample with cumulative coverage'
    sample = pd.read_csv(file,
                sep="\t", header=None, 
                names = ["chrom", "start", "end", "coverage"],
                dtype={0: str})
    logging.info("Loaded: " + file)
    return sample

def save_sample(filePath, pyobject):
    'Save python object to file'
    fileObj = open(filePath, 'wb')
    pickle.dump(pyobject, fileObj)
    fileObj.close()

class TE:
    "Object holding one terminal exon"
    def __init__(self, TE):
        self.chr = TE.chr
        self.start = TE.TE_start
        self.end = TE.TE_end
        self.id = TE.TE_id
        self.score = TE.score
        self.strand = TE.strand
        self.df = None
        self.max_coverage = None
        self.cpm = None
        self.norm_coverage = None
        self.coverage_per_bp = None
    
    def get_TE(self):
        'Return TE object'
        return pd.Series({"chr": self.chr, "TE_start": self.start, "TE_end": self.end, 
          "TE_id": self.id, "score": self.score, "strand": self.strand, 
          "max_coverage": self.max_coverage, "cpm": self.cpm})
    
    def set_max_coverage(self):
        return np.max(self.df.coverage)
    
    def set_sample(self, sample):
        "Get the cumulative coverage in the TE region. Ensure same strandedness"
        tmp = sample[sample.chrom == self.chr]
        tmp = tmp[tmp.start >= self.start]
        tmp = tmp[tmp.end <= self.end]
        self.df = tmp 
        if len(tmp) != 0:
            self.max_coverage = self.set_max_coverage()
        else:
            self.max_coverage = 0
    
    def set_norm_coverage(self):
        'Compute normalised coverage. Only valid for non-zero coverage. Else keeps None'
        if self.max_coverage != 0:
            self.norm_coverage = self.df.coverage / self.max_coverage
    
    def compute_coverage_per_bp(self):
        "Construct vector of length TE.TE_end - TE.TE_start"
        if self.norm_coverage is not None:
            out = np.empty(0, dtype=int)
            for i in self.df.index:
                arr = np.repeat(self.norm_coverage.loc[i], self.df.loc[i].end - self.df.loc[i].start)
                out = np.append(out, arr)
            self.coverage_per_bp = out
    
    def compute_auc(self, other_TE):
        "Compute Area under the curve by integration"
        # ensure the other TE is the same
        if self.id != other_TE.id:
            raise Exception("Not the same terminal exon!")
        if self.coverage_per_bp is None or other_TE.coverage_per_bp is None:
            return float('nan')
        if self.strand == "+":
            res = -1 * compute_area(self.coverage_per_bp, other_TE.coverage_per_bp)
        else:
            res = compute_area(self.coverage_per_bp, other_TE.coverage_per_bp)
        return res
    
    def get_TE_length(self):
        "Return length of TE in base pairs"
        return(self.end - self.start)
    def get_quantile_TE_length(self, quantile):
        "Compute and return the TE length for a given quantile in range [0-1]"
        # Will return length of 0 for quantile 1, and almost max TE length for 0
        quantile_TE_length = np.NaN
        # quantile TE length can only be computed if coverage > 1 (i.e. at least 2 rows in df)
        if self.max_coverage > 1 and self.norm_coverage is not None:
            if self.strand == '+':
                pos = self.df[self.norm_coverage >= quantile].iloc[-1].start
                quantile_TE_length = pos - self.start
            elif self.strand == '-':
                pos = self.df[self.norm_coverage >= quantile].iloc[0].end
                quantile_TE_length = self.end - pos
        return(quantile_TE_length)
    
    def info(self):
        return str("%s, strand: %s, max coverage: %d, cpm: %f" 
          % (self.id, self.strand, self.max_coverage, self.cpm))

class Sample_TE_coverage:
    """
    Object holding the cumulative coverage for one sample and n TEs
    The object is strand-specific.
    """
    def __init__(self, sample_df_plus, sample_df_minus):
        self.df_plus = sample_df_plus
        self.df_minus = sample_df_minus
        self.TEs = []
        self.norm_coverage = None
        self.out = None
    
    def add_TE(self, new_TE):
        "Add a terminal exon by creating an instance of class TE"
        
        # only add TE if ID not in list
        if new_TE.TE_id in self.get_ids():
            return
        new = TE(new_TE)
        if new.strand == "+":
            new.set_sample(self.df_plus)
        else:
            new.set_sample(self.df_minus)
        self.TEs.append(new)
    
    def get_max_coverage(self, TE_id):
        "Obtain the max coverage for TE_id"
        for te in self.TEs:
            if te.id == TE_id:
                return te.max_coverage
    
    def get_TE_object(self, TE_id):
        # use this function to plot the AUC plot for a given TE between samples
        'Obtain the reference for the TE_id'
        for te in self.TEs:
            if te.id == TE_id:
                return te
    
    def get_ids(self):
        return [te.id for te in self.TEs]
    
    def get_all_max_coverage(self):
        "Return df with max_coverage for all TEs"
        all_max_coverage = [te.max_coverage for te in self.TEs]
        all_ids = self.get_ids()
        return pd.DataFrame(data={'id': all_ids, 'max_coverage': all_max_coverage})
    
    def get_all_cpm(self):
        "Return df with cpm for all TEs"
        all_cpm = [te.cpm for te in self.TEs]
        all_ids = self.get_ids()
        return pd.DataFrame(data={'id': all_ids, 'cpm': all_cpm})
    
    def get_all_TEs(self):
        'Obtain all TEs including max_coverage and CPM'
        df = pd.DataFrame([te.get_TE() for te in self.TEs])
        return df
    def get_all_TE_lengths(self):
        'Obtain all TE lengths'
        all_te_lengths = [te.get_TE_length() for te in self.TEs]
        all_ids = self.get_ids()
        return(pd.DataFrame(data={'id': all_ids, 'max_coverage': all_te_lengths}))
    
    def get_quantiles(self, quantiles):
        "Obtain pd.dataframe of all TEs with additional columns of quantile TE lengths and TE length"
        df = pd.DataFrame(data={'id': self.get_ids()})
        qnames = ['length_q' + str(x) for x in quantiles]
        for i in range(len(quantiles)):
            df[qnames[i]] = [te.get_quantile_TE_length(quantiles[i]) for te in self.TEs]
        return(df)
    
    def set_cpm(self):
        'Compute counts per million for all available TEs'
        df = self.get_all_max_coverage()
        df['cpm'] = compute_cpm(df.max_coverage)
        for i in range(len(df)):
            self.get_TE_object(df.loc[i].id).cpm = df.loc[i].cpm
    
    def remove_TEs(self, TE_ids):
        "Remove TEs with given ID from TEs list"
        for TE_id in TE_ids:
            if TE_id in self.get_ids():
                entry = [te for te in self.TEs if te.id == TE_id][0]
                self.TEs.remove(entry)
    
    def set_all_norm_coverage(self):
        "Compute coverage per bp for all TEs"
        [te.set_norm_coverage() for te in self.TEs]
    
    def set_all_coverage_per_bp(self):
        "Compute coverage per bp for all TEs"
        [te.compute_coverage_per_bp() for te in self.TEs]
    
    def compute_auc_to(self, other_sample):
        "Compute AUC for all TEs and return dataframe"
        df1 = self.get_all_TEs()
        df1.rename(columns={'max_coverage':'coverage_sample1'}, inplace=True)
        df1.rename(columns={'cpm':'cpm_sample1'}, inplace=True)
        df2 = other_sample.get_all_TEs()
        df2.rename(columns={'max_coverage':'coverage_sample2'}, inplace=True)
        df2.rename(columns={'cpm':'cpm_sample2'}, inplace=True)
        df = df1.merge(df2, on = ['chr','TE_start','TE_end','TE_id','score','strand'])
        df['auc'] = np.full(len(df), float("nan"))
        for i in df.index:
            te1 = self.get_TE_object(df.loc[i,'TE_id'])
            te2 = other_sample.get_TE_object(df.loc[i,'TE_id'])
            df.loc[i, 'auc'] = te1.compute_auc(te2)
        return df
    
    def info(self):
        "Return some info about the object"
        return str("Nr. TEs: %s" % (len(self.TEs)))

#### main program
def main():
    parser = argparse.ArgumentParser(description="Compute the AUC of the TE coverage between two samples.")
    parser.add_argument('-bed', '--bed-file', dest='bed_file', 
      required=True, help='Input BED file')
    parser.add_argument('--one-plus', dest='sample1_plus', 
      required=True, help='Sample 1 cumulative coverage bedgraph plus strand')
    parser.add_argument('--one-minus', dest='sample1_minus', 
      required=True, help='Sample 1 cumulative coverage bedgraph minus strand')
    parser.add_argument('--two-plus', dest='sample2_plus', 
      required=True, help='Sample 2 cumulative coverage bedgraph plus strand')
    parser.add_argument('--two-minus', dest='sample2_minus', 
      required=True, help='Sample 2 cumulative coverage bedgraph minus strand')
    parser.add_argument('--out-file', dest='out_file', 
      required=False, help='Desination for main output file. Will create dirs as necessary. Default: calculation_auc_out/auc_sample1_sample2.tsv.', 
      default=os.path.join("calculation_auc_out", "auc_sample1_sample2.tsv"))
    parser.add_argument('-cpm', '--cpm-threshold', dest='cpm_threshold', 
      required=False, help='Exclude TEs with < cpm counts per million (Integer; calculated per sample). Default: 0, perform no filtering',
      default = 0, type = int)
    parser.add_argument('-quants', '--quantiles', dest='quantiles', 
      required=False, help='Obtain TE length at provided quantiles and save to TSV. Default: no quantiles used and no TSV created. Example: 0.4 0.5',
      nargs = "+", type = float)
    parser.add_argument('-log', '--log-file', dest="log_file", 
      required=False, help='Logging file. Default: calculation_auc.log.', default='calculation_auc.log')
    args = parser.parse_args()

    # create out dir if not exists
    if not os.path.exists(os.path.dirname(args.out_file)):
        os.makedirs(os.path.dirname(args.out_file))

    # logging settings: from snakemake
    logging.basicConfig(filename=args.log_file, level=logging.INFO, 
                    format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')

    # load Terminal Exon annotations
    TEs = pd.read_csv(args.bed_file,
        sep="\t", header=None, 
        names=["chr", "TE_start", "TE_end", "TE_id", "score", "strand"],
        dtype={"chr": str, "TE_start": np.int64, "TE_end": np.int64, "TE_id": str, "score": str, "strand": str})
    logging.info("Loaded: " + args.bed_file)
    
    # load sample files
    cpm_threshold = int(args.cpm_threshold)
    sample1_plus = load_sample(args.sample1_plus)
    sample1_minus = load_sample(args.sample1_minus)
    sample2_plus = load_sample(args.sample2_plus)
    sample2_minus = load_sample(args.sample2_minus)
    # create object holding the sample
    sample1 = Sample_TE_coverage(sample1_plus, sample1_minus)
    logging.info("Created sample1 object from: %s and %s" %(args.sample1_plus, args.sample1_minus))
    sample2 = Sample_TE_coverage(sample2_plus, sample2_minus)
    logging.info("Created sample2 object from: %s and %s" %(args.sample2_plus, args.sample2_minus))
    
    # add TEs to sample object
    for ind in TEs.index:
        te = TEs.loc[ind]
        sample1.add_TE(te)
        sample2.add_TE(te)
    logging.info("Added %d Terminal Exons to samples object." % len(TEs))

    # set CPM
    sample1.set_cpm()
    cpm1 = sample1.get_all_cpm()
    sample2.set_cpm()
    cpm2 = sample2.get_all_cpm()
    logging.info("Set CPM for sample1 and sample2.")

    # remove TEs with < X CPM
    if cpm_threshold == 0:
        logging.info("CPM threshold: %d, no filtering of TEs performed" % cpm_threshold)
    elif cpm_threshold > 0:
        ids1 = cpm1[cpm1.cpm < cpm_threshold].id.tolist()
        ids2 = cpm2[cpm2.cpm < cpm_threshold].id.tolist()
        # take intersection: only remove TE if CPM < threshold in both samples 
        ids = list(set(ids1).intersection(ids2))
        logging.info("Filter %d TEs with < %d CPM in both samples." %(len(ids), cpm_threshold))
        sample1.remove_TEs(ids)
        sample2.remove_TEs(ids)
        logging.info("Finished filtering TEs")
    else:
        logging.warn("No valid CPM threshold provided: %d. Will not perform filtering and continue." % cpm_threshold)

    # Compute and set normalised coverage
    sample1.set_all_norm_coverage()
    sample2.set_all_norm_coverage()
    logging.info("Computed normalised coverage per BED region.")

    # compute coverage per bp, based on norm coverage
    sample1.set_all_coverage_per_bp()
    sample2.set_all_coverage_per_bp()
    logging.info("Finished computing TE coverage per bp.")
    # save samples to file, in format sample_name_one.obj
    save_sample(os.path.join(os.path.dirname(args.out_file), 
        "_".join(os.path.basename(args.out_file).split(".")[0].split("_")[1:]) + '_sample1.obj'), sample1)
    save_sample(os.path.join(os.path.dirname(args.out_file), 
        "_".join(os.path.basename(args.out_file).split(".")[0].split("_")[1:]) + '_sample2.obj'), sample2)

    # compute AUC for each TE
    res = sample1.compute_auc_to(sample2)
    logging.info("Finished computing AUC between sample1 and sample2.")

    # get quantile TE lengths and write to file
    quantiles = list(args.quantiles)
    if len(quantiles) != 0:
        df1 = sample1.get_quantiles(quantiles)
        df1.columns = ['TE_id'] + ['sample1_length_q' + str(x) for x in quantiles]
        res = pd.merge(res, df1, how = 'left', on = 'TE_id')
        df2 = sample2.get_quantiles(quantiles)
        df2.columns = ['TE_id'] + ['sample2_length_q' + str(x) for x in quantiles]
        res = pd.merge(res, df2, how = 'left', on = 'TE_id')
        logging.info("Finished computing TE lengths at quantiles: %s" % quantiles)
    
    res.to_csv(args.out_file, sep = "\t", index = False)
    logging.info("Wrote results to: " + args.out_file)
    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()