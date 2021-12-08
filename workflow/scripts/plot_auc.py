#!/usr/bin/env python
"""
Script for obtaining AUC per annotated TE
Between sample1 and sample2
"""

import argparse
import os

import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

# load class objects and functions
from calculation_auc import TE, Sample_TE_coverage, compute_area, compute_cpm, load_sample, save_sample

__author__ = "Dominik Burri"
__date__ = "2021-02-03"
__license__ = "Apache"
__version__ = "0.0.1"
__maintainer__ = "Dominik Burri"
__email__ = "dominik.burri@unibas.ch"
__status__ = "Development"

def load_pyobject(filePath):
    'Load python object from file'
    fileObj = open(filePath, 'rb')
    pyobject = pickle.load(fileObj)
    fileObj.close()
    return(pyobject)

# TODO: load object classses and helper functions

#### main program
def main():
    parser = argparse.ArgumentParser(description="Plot AUC graphs for selected terminal exon (TE).")
    parser.add_argument('--out-dir', dest='out_dir', 
      required=False, 
      help='Specify output directory to save files. Default: plot_auc',
      default='plot_auc')
    parser.add_argument('--sample1', dest='sample1', 
      required=True, 
      help='Specify sample 1, generated by calculation_auc.py.')
    parser.add_argument('--sample2', dest='sample2', 
      required=True, 
      help='Specify sample 2, generated by calculation_auc.py.')
    parser.add_argument('--te-id', dest='te_id', 
      required=True, 
      help='Specify TE id to plot. Usually a refseq identifier, e.g. NM_176917.')
    args = parser.parse_args()
    # Parse arguments and load files
    out_dir = args.out_dir
    # create out_dir if not exists
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    # load python object
    sample1 = load_pyobject(args.sample1)
    sample2 = load_pyobject(args.sample2)

    # plot graph that is used to compute AUC
    te_id = args.te_id
    gene_activated = sample1.get_TE_object(te_id)
    gene_naive = sample2.get_TE_object(te_id)
    auc = gene_activated.compute_auc(gene_naive)
    mean_cpm = np.mean([gene_activated.cpm, gene_naive.cpm])
    title_text = "%s, AUC %.3f, mean CPM %.0f" %(te_id, auc, mean_cpm)
    plt.figure(figsize=[3.3,3])
    plt.plot(gene_activated.coverage_per_bp, gene_naive.coverage_per_bp, 'mo-', markersize = 3)
    plt.plot([0,1], [0,1], color = 'grey')
    plt.title(title_text)
    plt.savefig(os.path.join(out_dir, 'plot_auc_%s.png' % te_id))
    plt.close()

    # plot coverage per bp
    plt.figure(figsize=[4.4,4])
    plt.plot(range(1, len(gene_activated.coverage_per_bp)+1), gene_activated.coverage_per_bp, label = 'sample1', color = '#6666FF')
    plt.plot(range(1, len(gene_naive.coverage_per_bp)+1), gene_naive.coverage_per_bp, label = 'sample2', color = '#FF6666')
    plt.xlabel('position (increasing)')
    plt.ylabel('normalised cumulative fraction')
    plt.legend(loc="best")
    plt.savefig(os.path.join(out_dir, 'plot_fraction_%s.png' % te_id))
    plt.close()

if __name__ == "__main__":
    main()