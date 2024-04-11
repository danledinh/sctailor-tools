# sctailor-tools
Data processing and analysis tools for scTaILoR-seq

## Description
The tools in this repository are intended to create cell-by-transcript counts tables from raw data generated from scTaILoR-seq experiments.

## Requirements
### Hardware
Data processing was performed on a High-Performance Computing (HPC) platform running `CentOS 7.9`.

### Software
1) `anaconda3` version 5.0.1
2) `python` version 3.6.3
3) `pysam` version 0.21.0
4) `java` version 11.0.2
5) `nextflow` version 22.04.4
8) `wf-single-cell` version 0.1.5 (epi2me-labs)
9) `umi-tools` version 1.1.2 (conda environment)
10) `isoquant` version 3.3 (conda environment)
11) `samtools` version 1.18

## Installation
Installation time in the order of minutes.

1) Clone this repository
2) Install listed requirements
3) Clone `wf-single-cell` workflow
4) Create `umitools` conda environment: `conda env create -f umi_tools_conda_enviroment.yml`
5) Create `isoquant` conda environment: `conda env create -f isoquant_conda_enviroment.yml`

## Usage
Runtime in the order of hours to days depending on sequencing depth. Test data requires approximately 8 hours using 4 cores and 50 GB RAM. 

### Example dataset
Download `LR_3CL_cancer_R1_1.sub1000k.fastq.gz` using SRA-toolkit and BioProject accession PRJNA993664. See tutorial here: https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump. Then, downsample data to 1 million reads using SeqKit (https://bioinf.shenwei.me/seqkit/usage/#sample).

### Cell barcode (CB) and unique molecular identifier (UMI) assignment using `wf-single-cell`
Details found here: https://github.com/epi2me-labs/wf-single-cell
```
nextflow run epi2me-labs/wf-single-cell \
    -r v0.1.5 \
    -w ${OUTPUT}/${PREFIX}_workspace \
    -c ${CONFIG} \
    -profile singularity \
    --max_threads ${CORES} \
    --resources_mm2_max_threads ${CORES} \
    --fastq ${FQ} \
    --ref_genome_dir ${REFERENCE} \
    --out_dir ${OUTPUT}

# merge .bam intermediates
cd ${OUTPUT}/${PREFIX}/bams
samtools merge wf_SC.bam *.bam
samtools index wf_SC.bam
```
#### Output
Intermediate chromosome-specific .bam files are created from the `wf-single-cell` workflow. Then, `samtools` is required to merge these .bam files before UMI group assignment using the `umi-tools` package.

### UMI deduplication using `umi-tools`
Details found here: https://github.com/CGATOxford/UMI-tools
```
# prior to running `umi-tools`, filter out records missing a 12-n.t. UMI sequence.
python3 tidy_UMI.py ${OUTPUT}/${PREFIX}/bams/wf_SC.bam    # output has a .tidy.bam suffix

# load umitools conda env
source activate umitools

# tag merged bam
umi_tools group \
    --output-bam \
    --stdin=${OUTPUT}/${PREFIX}/bams/wf_SC.tidy.bam \
    --stdout=${OUTPUT}/${PREFIX}/bams/wf_SC.grouped.bam \
    --per-cell \
    --per-gene \
    --extract-umi-method=tag \
    --umi-tag=UB \
    --cell-tag=CB \
    --gene-tag=GN

# keep longest read in each UMI group (dedup_UMI.py included in git repo)
### Alternatively, use `umi_tools dedup`
python3 dedup_UMI.py ${OUTPUT}/${PREFIX}/bams/wf_SC.grouped.bam

samtools view \
-N ${OUTPUT}/${PREFIX}/bams/qname_umitools.txt \
-o ${OUTPUT}/${PREFIX}/bams/wf_SC.dedup.bam \
${OUTPUT}/${PREFIX}/bams/wf_SC.grouped.bam

samtools sort \
-o ${OUTPUT}/${PREFIX}/bams/wf_SC.dedup.sorted.bam
${OUTPUT}/${PREFIX}/bams/wf_SC.dedup.bam

samtools index ${OUTPUT}/${PREFIX}/bams/wf_SC.dedup.sorted.bam
```
#### Output
Resultant .bam file contains representative sequence for each UMI group.

### Transcript detection and quantitation using `isoquant`
Details found here: https://github.com/ablab/IsoQuant
```
# load isoquant conda env
source activate isoquant

# run isoquant
isoquant.py \
    --reference ${REFERENCE} \
    --genedb ${GTF} \
    --bam ${OUTPUT}/${PREFIX}/bams/wf_SC.dedup.sorted.bam \
    --data_type nanopore \
    --read_group tag:CB \
    -o ${OUTPUT}/isoquant
```
#### Output
The cell-by-gene and cell-by-transcript tables (`SAMPLE_ID.gene_grouped_counts.tsv` and `SAMPLE_ID.transcript_grouped_counts.tsv`, respectively) were used in downstream analyses. 

### Specialized analysis scripts and functions
See files with `analysis.` prefix.
1) `analysis.merge.py` = Regression-based matrix merge (two or more targeting panels on a single sample) 
2) `analysis.haplotype.py` = Haplotype determination
