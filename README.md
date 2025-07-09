# roquin-iclip
Repository associated to the paper "ROQUIN-1 binding to the transcriptome of T cells reveals conserved targets in the control of lymphopoiesis"

## Project structure
The repository is structured into two main folders:
- ``analyses/`` - contains the code for the pipelines developed to analyze the iCLIP data and for downstream tasks. This is divided into subfolders for each analysis step, such as peak calling, sequence conservation, and structure conservation, with their own README files with instructions on how to run them.
- ``figures/`` - contains the code to generate the figures presented in the paper, including the scripts for generating the plots and the processed input data used in them.

---

### System requirements
To run the full code (including genome alignment and peak calling) a Linux multi-core system with at least 32 GB of RAM is required. The code has been developed and tested on a CentOS 7/8 system but should run on any Linux distribution with the required dependencies installed.

## Setup

### Data collection
Download the experimental data associated to this work from GEO (accession no. GSE295259) and place it to corresponding folders in *data/*.
Download FASTA and annotation files (both .gtf and .gff3) for mouse and human genomes from [GENCODE](https://www.gencodegenes.org).
* Versions: Mouse Release M25 (GRCm38.p6), Human release 40 (GRCh38.p13)

### Installation

This repository runs on Python>=3.6 and R version 4.3.1 (2023-06-16) -- "Beagle Scouts". Most of the workflows are wrapped in Snakemake (https://snakemake.readthedocs.io/).

Other dependencies for this software were installed in conda environments and are listed in the respective environment files (.yml), available in each analyses/subfolder.

Environment creation is handled automatically by Snakemake for each workflow (when a Snakefile is available) or can be created manually using the provided environment files, e.g. for the peak calling workflow:
```
# from the peak calling folder
conda create -f env/env.roquin.yml
```

## Contributors

- Giulia Cantini - [@giulic3](https://github.com/giulic3)
- Lambert Moyon - [@lambosaur](https://github.com/lambosaur)

## How to cite
If you use this code in your research, please cite the following paper: