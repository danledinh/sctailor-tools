# sctailor-tools
Data processing and analysis tools for scTaILoR-seq

## Description
The tools in this repository are intended to create cell-by-transcript counts tables from raw data generated from scTaILoR-seq experiments.

## Requirements
### Hardware
InDev

### Software
1) `anaconda3` version 5.0.1
2) `python` version 3.6.3
3) `pysam` version 0.21.0
4) `java` version 11.0.2
5) `nextflow` version 22.04.4
8) `wf-single-cell` version 0.1.5 (epi2me-labs)
9) `umi-tools` version 1.1.2 (conda environment)
10) `isoquant` version 3.3 (conda environment)

## Installation
Installation time in the order of minutes.

1) Clone this repository
2) Install listed requirements
3) Clone `wf-single-cell` workflow
4) Create `umi-tools` conda environment: `conda env create -f umi_tools_conda_enviroment.yml`
5) Create `isoquant` conda environment: `conda env create -f isoquant_conda_enviroment.yml`

## Usage
Runtime in the order of hours to days depending on sequencing depth. Test data requires approximately XXX minute when using XXX cores. 

### Example dataset
InDev

### Cell barcode (CB) and unique molecular identifier (UMI) assignment using `wf-single-cell`
InDev

#### Output
InDev

### UMI deduplication using `umi-tools`
InDev

#### Output
InDev

### Transcript detection and quantitation using `isoquant`
InDev

#### Output
InDev
