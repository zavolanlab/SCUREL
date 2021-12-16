#!/usr/bin/env python
"""
Script for separating aligned reads by cell type annotation
Uses snakemake for input and output
"""

import logging
import os
import sys

import numpy as np
import pandas as pd

import pysam

__author__ = "Dominik Burri"
__date__ = "2019-11-13"
__license__ = "Apache"
__version__ = "0.0.1"
__maintainer__ = "Dominik Burri"
__email__ = "dominik.burri@unibas.ch"
__status__ = "Development"

# logging settings: from snakemake
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, 
                    format='%(asctime)s %(levelname)s:%(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S ')
logging.info("Started script.")

# read files
try:
        # load cell type annotation
        annotations = pd.read_csv(snakemake.input["annots"], sep=",")
        # load .bam
        bam_reader = pysam.AlignmentFile(snakemake.input["sample"], "rb")
except:
        logging.error("Problems with reading files.")
        sys.exit()

logging.info("Succesfully read files.")

# change barcode to match the one from bam (contains -1 at end)
if "-" not in annotations.cell[0]:
    annotations.cell = annotations.cell + "-1"

sample_name = snakemake.params["sample_name"]
logging.info("Sample name: %s" % sample_name)
# make dictionary of annotations
ct_dict = {}
for i, row in annotations.loc[(annotations["orig.ident"] == sample_name)].iterrows():
    ct_dict[row.cell] = row.cell_type
# create cell type specific bam files
bams_out = {}
cell_types = annotations.cell_type.unique()
cell_types_used = []
# create dictionary of cell type specific bam files
for i in range(len(snakemake.output)):
    ct = snakemake.output[i].split("/")[-1].split(".")[0]
    if ct not in cell_types:
        logging.error("Cell type %s not in the cell types of %s" % (ct, snakemake.input["annots"]))
        sys.exit(-1)
    bams_out[ct] = pysam.AlignmentFile(snakemake.output[i], "wb", 
        template=bam_reader)
    cell_types_used.append(ct)
    logging.info("Created file: %s" % snakemake.output[i])

# grep the barcodes from .bam
counter = 0
lines = 0
for aln in bam_reader:
    lines += 1
    if lines % 1e6 == 0: logging.debug("Nr. lines processed: %d" % lines)
    cb_tag = aln.get_tag("CB")
    # retrieve cell type annotation; should only be one cell barcode, enforce
    cell_type = ct_dict.get(cb_tag) # returns None if key not existant
    if cell_type is not None and cell_type in cell_types_used:
        # write read to appropriate file
        bams_out[cell_type].write(aln)
    else:
        counter += 1

bam_reader.close()
logging.info("Nr. reads not annotated to any cell type: %d" % counter)
for i in range(len(snakemake.output)):
    ct = snakemake.output[i].split("/")[-1].split(".")[0]
    bams_out[ct].close()
    logging.info("Finished writing to %s" % snakemake.output[i])

logging.info("End of file")

