# sctailor-tools
Data processing and analysis tools for scTaILoR-seq

## Description
The tools in this repository are intended to create cell-by-transcript counts tables from raw data generated from scTaILoR-seq experiments.

## Requirements
1) `python` version 3.6.3
2) `pysam` version 0.21.0
3) `wf-single-cell` version 0.1.5 from epi2me-labs using the workflow manager `Nextflow` version 22.04.4 which runs on `Java` version 11.0.2
4) `UMI-tools` version 1.1.2 which was run using an `Anaconda3` version 5.0.1  environment
5) `Isoquant` version 3.3 was also run using an `Anaconda3` version 5.0.1 environment

## Installation
1) Clone this repository
2) Install nextflow
3) Clone `wf-single-cell` workflow
4) Create `UMI-tools` conda environment: `conda env create -f umi_tools_conda_enviroment.yml`
5) Create `Isoquant` conda environment: `conda env create -f isoquant_conda_enviroment.yml`

## Usage
TBD
