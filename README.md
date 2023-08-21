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
```
nextflow run epi2me-labs/wf-single-cell -r v0.1.5 \
    -w ${OUTPUT}/${PREFIX}_workspace
    -c ${CONFIG} \
    -profile singularity \
    --max_threads {CORES} \
	  --resources_mm2_max_threads {CORES} \
    --fastq ${FQ} \
    --ref_genome_dir ${REFERENCE} \
    --out_dir ${OUTPUT}

# merge .bam intermediates
cd ${OUTPUT}/bams
samtools merge wf_SC.bam *.bam
samtools index wf_SC.bam
```

#### Output
InDev

### UMI deduplication using `umi-tools`
```
# tag merged bam
umi_tools group -I ${OUTPUT}/wf_SC.bam --group-out=grouped.tsv --output-bam --log=group.log --paired

# keep longest read in each UMI group
python3 dedup_UMI.py ${OUTPUT}/grouped_sorted.bam
```

#### Output
InDev

### Transcript detection and quantitation using `isoquant`
```
isoquant.py --reference ${REFERENCE} --genedb ${GTF} --bam ${OUTPUT}/bamf.bam --data_type (nanopore) -o ${OUTPUT}
```

#### Output
InDev
