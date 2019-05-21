---
title: Workflows for analyzing single-cell RNA-seq data with R/Bioconductor
author: 
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
- name: Davis J. McCarthy
  affiliation: 
  - &EMBL EMBL European Bioinformatics Institute, Wellcome Genome Campus, Hinxton, Cambridge CB10 1SD, United Kingdom
  - St Vincent's Institute of Medical Research, 41 Victoria Parade, Fitzroy, Victoria 3065, Australia
- name: John C. Marioni
  affiliation: 
  - *CRUK
  - *EMBL
  - Wellcome Trust Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridge CB10 1SA, United Kingdom
date: "2019-05-20"
vignette: >
  %\VignetteIndexEntry{01. Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    toc_float: true
bibliography: ref.bib
---



# Workflow version information

**R version**: R version 3.6.0 Patched (2019-05-02 r76458)

**Bioconductor version**: 3.10

**Package**: 1.9.3

# Motivation

Single-cell RNA sequencing (scRNA-seq) is widely used to measure the genome-wide expression profile of individual cells.
From each cell, mRNA is isolated and reverse transcribed to cDNA for high-throughput sequencing [@stegle2015computational].
This can be done using microfluidics platforms like the Fluidigm C1 [@pollen2014lowcoverage], protocols based on microtiter plates like Smart-seq2 [@picelli2014fulllength], or droplet-based technologies like inDrop [@klein2015droplet;@macosko2015highly].
The number of reads mapped to each gene is then used to quantify its expression in each cell.
Alternatively, unique molecular identifiers (UMIs) can be used to directly measure the number of transcript molecules for each gene [@islam2014quantitative].
Count data are analyzed to detect highly variable genes (HVGs) that drive heterogeneity across cells in a population, to find correlations between genes and cellular phenotypes, or to identify new subpopulations via dimensionality reduction and clustering. 
This provides biological insights at a single-cell resolution that cannot be achieved with conventional bulk RNA sequencing of cell populations.

Strategies for scRNA-seq data analysis differ markedly from those for bulk RNA-seq.
One technical reason is that scRNA-seq data are much noisier than bulk data [@brennecke2013accounting;@marinov2014singlecell].
Reliable capture (i.e., conversion) of transcripts into cDNA for sequencing is difficult with the low quantity of RNA in a single cell.
This increases the frequency of drop-out events where none of the transcripts for a gene are captured.
Dedicated steps are required to deal with this noise during analysis, especially during quality control.
In addition, scRNA-seq data can be used to study cell-to-cell heterogeneity, e.g., to identify new cell subtypes, to characterize differentiation processes, to assign cells into their cell cycle phases, or to identify HVGs driving variability across the population [@vallejos2015basics;@fan2016characterizing;@trapnell2014dynamics].
This is simply not possible with bulk data, meaning that custom methods are required to perform these analyses. 

# scRNA-seq data analysis with Bioconductor

This package contains a set of computational workflows for basic analysis of scRNA-seq data, using software from the open-source Bioconductor project [@huber2015orchestrating].
The workflows start from a count matrix and describe a number of key steps for scRNA-seq data analysis, including:

- quality control to remove problematic cells;
- normalization of cell-specific biases, with and without spike-ins; 
- correction for batch effects; 
- cell cycle phase classification from gene expression data; 
- data exploration to identify putative subpopulations; 
- and finally, HVG and marker gene identification to prioritize interesting genes.

The application of these procedures will be demonstrated on several public scRNA-seq datasets involving immortalized myeloid progenitors, brain cells, haematopoietic stem cells, T-helper cells and mouse embryonic stem cells, generated with a range of experimental protocols and platforms [@lun2017assessing;@wilson2015combined;@zeisel2015brain;@islam2011characterization;@buettner2015computational;@zheng2017massively].
The aim is to provide a variety of modular usage examples that can be applied by readers to construct custom analysis pipelines for their own experiments.

See the *[simpleSingleCell](https://bioconductor.org/packages/3.10/simpleSingleCell)* landing page for links to individual workflows and for instructions on how to install the required packages.
To cite any of these workflows, please refer to http://f1000research.com/articles/5-2122/v2 for instructions.

# Obtaining a count matrix

All of these workflows start from a publicly available count matrix.
For simplicity, we forego a description of the read processing steps required to generate the count matrix, i.e., read alignment and counting into features.
For SMART-seq2 data [@picelli2014fulllength], quantification procedures developed for bulk RNA-seq are generally satisfactory [@love2015rnaseq;@chen2016from].
Users favouring an R-based approach to read alignment and counting might consider using the methods in the *[Rsubread](https://bioconductor.org/packages/3.10/Rsubread)* package [@liao2013subread;@liao2014featurecounts].

Many other scRNA-seq protocols contain bespoke sequence structures that require careful processing:

- Unique molecular identifiers (UMIs) [@islam2014quantitative] are widely used to mitigate the effects of amplification biases.
Reads with the same UMI mapping to the same gene represent a single underlying transcript molecule and only increment the count of that gene by one.
Processing of this data requires extraction of the UMI sequence from each read or read pair, and a method to collapse UMI-based duplicates into a single count [@smith2017umitools].
- Protocols may also use custom cell barcodes to improve multiplexing efficiency beyond that offered by the standard sequencing barcodes (e.g., from Illumina).
This includes data generated from droplet-based experiments [@zheng2017massively] or from very high-throughput plate-based protocols like MARS-seq [@jaitin2014massively].
Processing of this data usually requires a separate step to extract the barcode sequence from each read and to allocate the read to the correct per-cell sequencing library.

For these data sets, the *[scPipe](https://bioconductor.org/packages/3.10/scPipe)* package [@lian2018scpipe] provides an R-based processing pipeline for obtaining a count matrix. 

If spike-in RNA was added, the sequences of the spike-in transcripts can be included as additional FASTA files during genome index building prior to alignment.
Similarly, genomic intervals for both spike-in transcripts and endogenous genes can be concatenated into a single GTF file prior to counting.

# Author information

## Author contributions

A.T.L.L. developed and tested workflows on all datasets.
A.T.L.L. and D.J.M. implemented improvements to the software packages required by the workflow.
J.C.M. provided direction to the software and workflow development.
All authors wrote and approved the final manuscript.

## Competing interests

No competing interests were disclosed.

## Grant information

A.T.L.L. and J.C.M. were supported by core funding from Cancer Research UK (award no. A17197).
D.J.M. was supported by a CJ Martin Fellowship from the National Health and Medical Research Council of Australia.
D.J.M and J.C.M. were also supported by core funding from EMBL.

## Acknowledgements

We would like to thank Antonio Scialdone for helpful discussions, as well as Michael Epstein, James R. Smith and John Wilson-Kanamori for testing the workflow on other datasets.

# References

