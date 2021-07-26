<!-- 
Title: sc-polyAsite description 
Author: Dominik Burri
Date: 2021-06-29 
-->

# Overview

SCUREL (Single Cell 3'Untranslated REgion Length analysis) performs differential 3'UTR length analysis on the level of group of cells.
Check the preprint for more information: https://doi.org/10.1101/2021.06.30.450496.

Starting from scRNA-seq data generated with the 10X Genomics platform and cell type annotations, the framework detects 3'UTR length change events by computing the Area under the curve (AUC) for a given terminal exon region (TE) between two samples and reports TEs with significant changes, separate for 3'UTR shortening and lengthening.

Two execution modes are possible:
- `cell_type_comparison`: compare cell types specified in a separate table. Cell types with or without `sample_origin` can be compared.
- `cell_state_comparison`: compare same cell type from different *sample_origin*, both of which have to be specified for each sample in the samples table `config/samples.tsv`. This mode will also create a quasi-bulk sample by gathering all reads from sample origin, irrespective of cell type. This additional analysis is called *merged* in the results.

The framework is written as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline and therefore individual start- and end-points are possible. The pipeline uses [conda](https://docs.conda.io/en/latest/) for creating virtual environments for individual rules. 

This repository also contains a script for cell type annotation based on marker genes using [Seurat](https://satijalab.org/seurat/).

## Input

* A configuration file (`config/config.yaml`) with fields for file paths, directory locations and parameters. This file and all other files can be copied, adjusted and renamed. It is important to adjust all file paths in this configuration file.
* A samples table (see skeleton `config/samples.tsv`) with tab separated columns **sample**, **name**, **fastqs** and **origin** for site of origin. The file is specified in field **defsamples**.
  * The field **name** is used for the sample names as occurring in the 10X sequencing files. The actual file names must follow cellrangers naming convention (see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input)).
  * **fastqs** denotes the directory where the samples reside. 
  * Mapped reads (BAM) obtained with *cellranger count* can be used. The BAM files need to be stored in directory `cellranger_count/{sample}/outs/`.
  * Alternatively, samples can be provided as raw reads (FASTQ) from 10x Genomics 3' end sc-RNA-seq. A helper rule can execute *cellranger count*, for which the reference transcriptome for *cellranger* is necessary. *cellranger* must be obtained separately and can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).
* Cell type annotations (`cell_type_annotations.csv`), a comma separated table with assignment of each cell to a cell type. For more details on how to obtain this, see section *Prerequesites* below.
* Genome annotation for obtaining terminal exon regions. Can be independent of reference transcriptome from *cellranger*. 
  * The script `scripts/build_refseq.py` is intended to build a suitable genome annotation file from Refseq reads mapped with *cellranger*. In particular, it changes chromosome names to the one used by *cellranger* and subsets by removing other chromosomes such as unlocalized-scaffolds or fixed-patches. Additionally, it filters by predicted (Gnomon) entries.


## Output

All intermediate files and results are saved in the directory specified in field **out_dir** in `config/config.yaml`.
Results are saved in different subdirectories for execution modes:

| mode | subdirectory |
| --- | --- |
|`cell_type_comparison` | `auc_comparisons`|
| `cell_state_comparison` | `auc` |

Each subdirectory contains a directory `analysis_out` with plots and list of TEs with significant changes in 3'UTR length.

# Installation

Installation procedure on Linux (tested on CentOS 7).

## Step 1: Clone this pipeline

Clone the pipeline and follow the instructions below. This requires *git* installed and configured (see [github git cheat sheet](https://github.github.com/training-kit/downloads/github-git-cheat-sheet.pdf)).

```bash
git clone https://github.com/zavolanlab/SCUREL.git
cd SCUREL
```

## Step 2: Install Miniconda 3

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source .bashrc
```

Test installation with

```bash
conda list
```

## Step 3: Install snakemake

The easiest way to install snakemake and all dependencies for this pipeline, is to create and activate the conda environment named *scurel* from the supplied `envs/scurel.yaml` file via

```bash
conda env create --name scurel --file envs/scurel.yaml
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

# Prerequesites

## Map reads

If one starts from FASTQ, the reads need to be mapped. We provide a wrapper for *cellranger* in `workflow/mapping.smk`. The rules map FASTQ files in the given directories set in `config/samples.tsv` and filters reads. 

In the `Snakefile` put the following as finish rule:

```bash
rule finish:
    input: 
      expand(os.path.join("cellranger_count", "{sample}", 
        "outs", "possorted_genome_bam.bam"),
        sample = samples.index)
```

## Annotate cell barcodes

In order to run the pipeline, a file named `cell_type_annotations.csv` with columns **orig.ident** (identical to **sample**), **cell**, **cell_type** is needed. 

One way to create this is by running `sc_quant_and_cell_types.Rmd` step-by-step.
> **Important**: The steps below require manual execution! 
>
> Also, in order to define cell types of interest marker genes for each cell type are needed.

### Cell type identification with marker genes

The procedure in `sc_quant_and_cell_types.Rmd` implements single cell gene expression quantification and determination of the major cell types based on [Lambrechts et al. 2018](https://doi.org/10.1038/s41591-018-0096-5). The best way to explore it is by using [RStudio](https://www.rstudio.com/). It is adviced to carefully run the script on new datasets and check the parameters and adjust as needed (e.g. the dimensionality of the dataset).
The output will be a file containing for each cell its cell type annotation (`cell_type_annotations.csv`).


# Pipeline execution

The following files should be copied and adjusted:
- `config/config.yaml`
- `config/samples.tsv`

## Dry run

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

## Local

When executing the pipeline locally (i.e. on personal machine), *cellranger v3* has to be installed to map 10x raw fastq files to the genome: [link](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation).
Please note that *cellranger* requires a minimum of 64GB RAM (see [system requirements](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/system-requirements)).

Execute *cellranger* as standalone pipeline and generate BAM files of mapped reads. To ensure the subsequent execution of this pipeline, create a directory `cellranger_count` inside the pipeline directory. Otherwise, adjust the relative paths for rule *filter_high_quality* in `workflow/mapping.smk`.

Run in the directory of this pipeline and specify CONFIGFILE

```bash
bash exe_snakemake_local.sh CONFIGFILE
```

## Slurm

The pipeline is currently tailored to execution on a slurm cluster with CentOS 7. Like this, only `config/config.yaml` has to be adjusted.
If necessary, adjust queue and time constraints in `config/cluster.json`. 
Additionally, memory and thread usage can be changed in the individual snakemake rules.

> Note that *cellranger* is loaded manually on the slurm cluster. This has to be specified within the file `workflow/mapping.smk` inside the rule *cellranger_count* as a first statement in the shell command (e.g. `ml CellRanger/5.0.0;`).

Run in the directory of this pipeline and specify CONFIGFILE

```bash
bash exe_snakemake_cluster.sh CONFIGFILE
```

## Remove output

With

```bash
bash remove_logs.sh
```

all log files from the snakemake runs will be removed.

# Analysis

The following scripts are used for the pathway analysis.

* `scripts/analysis/patient_comparison.Rmd`: Interactive script to compare datasets or patients for the number of TEs with shortening and lengthening. Converts TE ids to gene names.
* `scripts/convert_transcripts.R`: Command line script to convert Refseq identifiers to gene names.
* `scripts/analysis/comparison_with_scapa.R`: Interactive script to compare results of SCUREL to scAPA. 
* `scripts/plotting/plot_enrichment.R`: Command line script to plot pathway enrichment.
* `scripts/plotting/get_bg_genome.R`: Command line script to get background genome from expressed transcripts.
* `scripts/plot_auc.py`: Command line script to plot graph to compute AUC.
* `scripts/plotting/heatmap_enrichment.R`: Command line script to plot pathway enrichment for multiple pathways as heatmap.
* `scripts/analysis/table_all_APA.R`: Interactive script to gather all significant APA shortening events and create binary table. 
