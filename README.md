<!-- 
Title: SCUREL description
Author: Dominik Burri
Date: 2021-06-29 
-->

<img align="right" width="50" height="50" src="images/logo.128px.png">

# SCUREL

[![ci](https://github.com/zavolanlab/SCUREL/workflows/ci/badge.svg?branch=main)](https://github.com/zavolanlab/SCUREL/actions?query=workflow%3Aci)
[![GitHub issues](https://img.shields.io/github/issues/zavolanlab/scurel)](https://github.com/zavolanlab/scurel/issues)
[![GitHub license](https://img.shields.io/github/license/zavolanlab/scurel)](https://github.com/zavolanlab/SCUREL/blob/main/LICENSE)

SCUREL is a method to detect 3'UTR changes from scRNA-seq data. It is built as a Snakemake workflow and uses virtual environments for execution.

# Table of Contents

- [SCUREL](#scurel)
  - [Table of Contents](#table-of-contents)
  - [General information](#general-information)
    - [Input](#input)
    - [Output](#output)
  - [Installation](#installation)
    - [Step 1: Clone this pipeline](#step-1-clone-this-pipeline)
    - [Step 2: Install Miniconda 3](#step-2-install-miniconda-3)
    - [Step 3: Install Snakemake](#step-3-install-snakemake)
  - [Prerequesites](#prerequesites)
    - [Map reads](#map-reads)
    - [Annotate cell barcodes](#annotate-cell-barcodes)
  - [Pipeline execution](#pipeline-execution)
    - [Dry run](#dry-run)
    - [Local](#local)
    - [Slurm](#slurm)
    - [Remove output](#remove-output)
  - [Analysis](#analysis)
  - [Contributing](#contributing)
  - [Contact](#contact)

## General information

SCUREL (Single Cell 3'Untranslated REgion Length analysis) performs differential 3'UTR length analysis on the level of group of cells.
Take a look at the publication for more information regarding the method and its application: http://dx.doi.org/10.1261/rna.078886.121.

Starting from scRNA-seq data generated with the 10X Genomics platform and cell type annotations, the framework detects 3'UTR length change events by computing the Area under the curve (AUC) for a given terminal exon region (TE) between two samples and reports TEs with significant changes, separate for 3'UTR shortening and lengthening.

Two execution modes are possible:
- `cell_type_comparison`: compare cell types specified in a separate table. Cell types with or without `sample_origin` can be compared.
- `cell_state_comparison`: compare same cell type from different *sample_origin*, latter of which has to be specified for each sample in the samples table `config/samples.tsv`. Only pairwise comparisons are possible and only the first two entries are considered. If more than two sample origins available, perform them separately. The entries are order-sensitive, i.e. the analysis will always compare the first versus the second entry. This mode will also create a quasi-bulk sample by gathering all reads from sample origin, irrespective of cell type. This additional analysis is called *merged* in the results.

The framework is written as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and therefore individual start- and end-points are possible. The workflow uses [conda](https://docs.conda.io/en/latest/) for creating virtual environments for individual rules. 

This repository also contains a script for cell type annotation based on marker genes using [Seurat](https://satijalab.org/seurat/).

### Input

* A configuration file (`config/config.yaml`) with fields for file paths, directory locations and parameters. This file and all other files can be copied, adjusted and renamed. It is important to adjust all file paths in this configuration file.
* A samples table (see skeleton `config/samples.tsv`) with tab separated columns **sample**, **name**, **fastqs** and **origin** for site of origin. The file is specified in field **defsamples**.
  * The field **name** is used for the sample names as occurring in the 10X sequencing files. The actual file names must follow cellrangers naming convention (see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input)).
  * **fastqs** denotes the directory where the samples reside. 
  * Mapped reads (BAM) obtained with *cellranger count* can be used. The BAM files need to be stored in directory `cellranger_count/{sample}/outs/`.
  * Alternatively, samples can be provided as raw reads (FASTQ) from 10x Genomics 3' end sc-RNA-seq. A helper rule can execute *cellranger count*, for which the reference transcriptome for *cellranger* is necessary. *cellranger* must be obtained separately and can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).
* Cell type annotations (`cell_type_annotations.csv`), a comma separated table with assignment of each cell to a cell type. For more details on how to obtain this, see section *Prerequesites* below.
* Genome annotation for obtaining terminal exon regions. Can be independent of reference transcriptome from *cellranger*. 
  * The script `scripts/build_refseq.py` is intended to build a suitable genome annotation file from Refseq reads mapped with *cellranger*. In particular, it changes chromosome names to the one used by *cellranger* and subsets by removing other chromosomes such as unlocalized-scaffolds or fixed-patches. Additionally, it filters by predicted (Gnomon) entries.


### Output

All intermediate files and results are saved in the directory specified in field **out_dir** in `config/config.yaml`.
Results are saved in different subdirectories for execution modes:

| mode | subdirectory |
| --- | --- |
|`cell_type_comparison` | `auc_comparisons`|
| `cell_state_comparison` | `auc` |

Each subdirectory contains a directory `analysis_out` with plots and list of TEs with significant changes in 3'UTR length. The change in 3'UTR length is always reported in respect to the first sample (either cell type or first entry in *sample_origin*). For example, one compares T cells in tumor versus control tissue. Then, the AUC analyses that report 3'UTR lengthening mean that this lengthening in observed in tumor tissue in comparison to control tissue. Also, AUC values < 0.5 indicate lengthening in the first sample (tumor tissue). The same applies when performing cell type comparison, the first cell type is compared against the second one.

## Installation

Installation procedure on Linux (tested on CentOS 7).

### Step 1: Clone this pipeline

Clone the pipeline and follow the instructions below. This requires *git* installed and configured (see [github git cheat sheet](https://github.github.com/training-kit/downloads/github-git-cheat-sheet.pdf)).

```bash
git clone https://github.com/zavolanlab/SCUREL.git
cd SCUREL
```

### Step 2: Install Miniconda 3

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source .bashrc
```

Test installation with

```bash
conda list
```

#### Install mamba

[Mamba](https://mamba.readthedocs.io/en/latest/installation.html) is a faster way of installing conda packages.
Install with

```bash
conda install mamba -n base -c conda-forge
```

Installation of mamba is not necessary. If not installed, replace `mamba` commands by `conda`.

### Step 3: Install Snakemake

The easiest way to install snakemake and all dependencies for this pipeline, is to create and activate the conda environment named *scurel* from the supplied `install/scurel.yaml` file via

```bash
mamba env create --name scurel --file install/scurel.yaml
conda activate scurel
```

Alternatively, snakemake can be installed as standalone package into a new environment *scurel* via

```bash
conda create -n scurel -c bioconda snakemake
```

Then the environment must be activated with

```bash
conda activate scurel
```

## Prerequesites

### Map reads

If one starts from FASTQ, the reads need to be mapped. We provide a wrapper for *cellranger* in `workflow/rules/mapping.smk`. The rules map FASTQ files in the given directories set in `config/samples.tsv` and filters reads. 

In the `Snakefile` put the following as finish rule:

```bash
rule finish:
    input: 
      expand(os.path.join("cellranger_count", "{sample}", 
        "outs", "possorted_genome_bam.bam"),
        sample = samples.index)
```

### Annotate cell barcodes

In order to run the pipeline, a file named `cell_type_annotations.csv` with columns **orig.ident** (identical to **sample**), **cell**, **cell_type** is needed. 

One way to create this is by running `sc_quant_and_cell_types.Rmd` step-by-step.
> **Important**: The steps below require manual execution! 
>
> Also, in order to define cell types of interest marker genes for each cell type are needed.

#### Cell type identification with marker genes

The procedure in `sc_quant_and_cell_types.Rmd` implements single cell gene expression quantification and determination of the major cell types based on [Lambrechts et al. 2018](https://doi.org/10.1038/s41591-018-0096-5). The best way to explore it is by using [RStudio](https://www.rstudio.com/). It is adviced to carefully run the script on new datasets and check the parameters and adjust as needed (e.g. the dimensionality of the dataset).
The output will be a file containing for each cell its cell type annotation (`cell_type_annotations.csv`).


## Pipeline execution

The following files should be copied and adjusted:
- `config/config.yaml`
- `config/samples.tsv`

### Dry run

A dry run can be executed to check if everything is in place and the correct directories are used.

Run in the directory of this pipeline and replace CONFIGFLIE with an actual file

```bash
bash exe_dry.sh CONFIGFILE
```

> Note: the bash script `exe_dry.sh` takes one command line argument, namely the configuration file.

Additionally, the directed acyclic graph (DAG) and the rulegraph for CONFIGFILE can be visualised by executing the following lines respectively. 

```bash
bash build_dag.sh CONFIGFILE
bash build_rulegraph.sh CONFIGFILE
```

### Local

When executing the pipeline locally (i.e. on personal machine), *cellranger v3* has to be installed to map 10x raw fastq files to the genome: [link](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation).
Please note that *cellranger* requires a minimum of 64GB RAM (see [system requirements](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/system-requirements)).

Execute *cellranger* as standalone pipeline and generate BAM files of mapped reads. To ensure the subsequent execution of this pipeline, create a directory `cellranger_count` inside the pipeline directory. Otherwise, adjust the relative paths for rule *filter_high_quality* in `workflow/rules/mapping.smk`.

Run in the directory of this pipeline and specify CONFIGFILE

```bash
bash exe_snakemake_local.sh CONFIGFILE
```

### Slurm

The pipeline is currently tailored to execution on a slurm cluster with CentOS 7. Like this, only `config/config.yaml` has to be adjusted.
If necessary, adjust queue and time constraints in `config/cluster.json`. 
Additionally, memory and thread usage can be changed in the individual snakemake rules.

> Note that *cellranger* is loaded manually on the slurm cluster. This has to be specified within the file `workflow/rules/mapping.smk` inside the rule *cellranger_count* as a first statement in the shell command (e.g. `ml CellRanger/5.0.0;`).

Run in the directory of this pipeline and specify CONFIGFILE

```bash
bash exe_snakemake_cluster.sh CONFIGFILE
```

### Remove output

With

```bash
bash remove_logs.sh
```

all log files from the snakemake runs will be removed.

## Analysis

The following scripts are used for the pathway analysis.

* `workflow/scripts/analysis/dataset_and_patient_comparison.Rmd`: Interactive script to compare datasets or patients for the number of TEs with shortening and lengthening. Converts TE ids to gene names.
* `workflow/scripts/convert_transcripts.R`: Command line script to convert Refseq identifiers to gene names.
* `workflow/scripts/analysis/comparison_with_scapa.R`: Interactive script to compare results of SCUREL to scAPA. 
* `workflow/scripts/plotting/plot_enrichment.R`: Command line script to plot pathway enrichment.
* `workflow/scripts/plotting/get_bg_genome.R`: Command line script to get background genome from expressed transcripts.
* `workflow/scripts/plot_auc.py`: Command line script to plot graph to compute AUC.
* `workflow/scripts/plotting/heatmap_enrichment.R`: Command line script to plot pathway enrichment for multiple pathways as heatmap.
* `workflow/scripts/analysis/table_all_APA.R`: Interactive script to gather all significant APA shortening events and create binary table. 

## Contributing

This project lives off your contributions, be it in the form of bug reports, feature requests, discussions, or fixes and other code changes. Please refer to the [contributing](CONTRIBUTING.md) guidelines if you are interested to contribute. Please mind the [code of conduct](CODE_OF_CONDUCT.md) for all interactions with the community.

## Contact

For questions or suggestions regarding the code, please use the [issue tracker](https://github.com/zavolanlab/SCUREL/issues). For any other inquiries, please contact us by email: zavolab-biozentrum@unibas.ch.

2021 [Zavolab, Biozentrum, University of Basel](https://zavolan.biozentrum.unibas.ch)