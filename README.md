# Evolutionary constraint and deleterious variation across grass genomes

Code associated with:

**Liu et al. — *Deep evolutionary constraint predicts deleterious variation in grass genomes***

## Overview

This repository contains scripts used to analyze evolutionary constraint and putatively deleterious variation in:

* *Phyllostachys edulis*
* *Zea mays* / teosinte
* *Aegilops tauschii*

The study integrates comparative genomics across 94 grass species with population genomic data from the three focal species.

## Workflow

1. Grass phylogeny
2. Genome-wide synteny
3. Whole-genome alignment
4. GERP-based evolutionary constraint
5. Population SNP filtering
6. Deleterious variant annotation
7. ROH, FROH, and individual burden
8. Synteny-based gene comparisons

## Repository structure

```text
jobs/       LSF job scripts
scripts/    analysis and plotting scripts
metadata/   sample and species information
config/     software and path configuration
```

## Manuscript figures

| Figure     | Analysis                               |
| ---------- | -------------------------------------- |
| Figure 1   | Phylogeny and synteny                  |
| Figure 2   | Whole-genome alignment                 |
| Figure 3–4 | Evolutionary constraint                |
| Figure 5   | Putatively deleterious variation       |
| Figure 6   | Heterozygosity and ROH                 |
| Figure 7   | Synteny-based constraint and variation |

## Main thresholds

```text
GERP > 2   evolutionary constraint
GERP > 4   main deleterious-variant analysis
GERP > 6   high-stringency analysis

SIFT < 0.05   deleterious
```

## Major software

ROADIES, treePL, GeneTribe, JCVI, LAST, MULTIZ/ROAST, PHAST, GERP++, BWA-MEM, GATK, VCFtools, PLINK, SIFT 4G, and SnpEff.

Exact parameters are provided in the manuscript Methods.

## Data

Large genomic datasets are not included in this repository. Genome assemblies and population resequencing datasets are described in the manuscript and supplementary tables.

## Citation

Please cite the associated manuscript:

```text
Liu L-M, et al.
Deep evolutionary constraint predicts deleterious variation in grass genomes.
```
