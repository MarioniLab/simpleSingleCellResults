---
title: Analyzing single-cell RNA-seq data containing read counts
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
date: "2019-02-28"
vignette: >
  %\VignetteIndexEntry{02. Read count data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---



# Overview

In this workflow, we use a relatively simple dataset [@lun2017assessing] to introduce most of the concepts of scRNA-seq data analysis.
This dataset contains two plates of 416B cells (an immortalized mouse myeloid progenitor cell line), processed using the Smart-seq2 protocol [@picelli2014fulllength].
A constant amount of spike-in RNA from the External RNA Controls Consortium (ERCC) was also added to each cell's lysate prior to library preparation.
High-throughput sequencing was performed and the expression of each gene was quantified by counting the total number of reads mapped to its exonic regions.
Similarly, the quantity of each spike-in transcript was measured by counting the number of reads mapped to the spike-in reference sequences.

Counts for all genes/transcripts in each cell are available from ArrayExpress using the accession number [E-MTAB-5522](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5522).
We download both the count tables (in the "processed files") as well as the metadata file using the *[BiocFileCache](https://bioconductor.org/packages/3.9/BiocFileCache)* package. 
This saves the files to a local cache (`raw_data`) and avoids re-downloading them if they are already present.


```r
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
lun.zip <- bfcrpath(bfc, 
    file.path("https://www.ebi.ac.uk/arrayexpress/files",
        "E-MTAB-5522/E-MTAB-5522.processed.1.zip"))
lun.sdrf <- bfcrpath(bfc, 
    file.path("https://www.ebi.ac.uk/arrayexpress/files",
        "E-MTAB-5522/E-MTAB-5522.sdrf.txt"))
unzip(lun.zip, exdir=tempdir())
```

<!--
We don't use a global cache to ensure that the package build is self-contained on the Bioc build system and elsewhere.
This avoids problems with cross-package behaviour.
It also ensures that a new file is downloaded at every build, which is necessary for checking that the files are current.
We avoid hosting an extra copy of everything (locally) by soft-linking raw_data to ~/.cache/BiocFileCache/.
-->

# Setting up the data

## Loading in the count matrix

Our first task is to load the count matrices into memory.
One matrix was generated for each plate of cells used in the study.
In each matrix, each row represents an endogenous gene or a spike-in transcript, and each column represents a cell.
Subsequently, the count in each entry of the matrix represents the number of reads mapped to a particular gene/transcript in a particular cell.


```r
plate1 <- read.delim(file.path(tempdir(), "counts_Calero_20160113.tsv"), 
    header=TRUE, row.names=1, check.names=FALSE)
plate2 <- read.delim(file.path(tempdir(), "counts_Calero_20160325.tsv"), 
    header=TRUE, row.names=1, check.names=FALSE)

gene.lengths <- plate1$Length # First column is the gene length.
plate1 <- as.matrix(plate1[,-1]) # Discarding gene length (as it is not a cell).
plate2 <- as.matrix(plate2[,-1])
rbind(Plate1=dim(plate1), Plate2=dim(plate2))
```

```
##         [,1] [,2]
## Plate1 46703   96
## Plate2 46703   96
```

We combine the two matrices into a single object for further processing.
This is done after verifying that the genes are in the same order between the two matrices.


```r
stopifnot(identical(rownames(plate1), rownames(plate2)))
all.counts <- cbind(plate1, plate2)
```

For convenience, the count matrix is stored in a `SingleCellExperiment` object from the *[SingleCellExperiment](https://bioconductor.org/packages/3.9/SingleCellExperiment)* package.
This allows different types of row- and column-level metadata to be stored alongside the counts for synchronized manipulation throughout the workflow.


```r
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=all.counts))
rowData(sce)$GeneLength <- gene.lengths
sce
```

```
## class: SingleCellExperiment 
## dim: 46703 192 
## metadata(0):
## assays(1): counts
## rownames(46703): ENSMUSG00000102693 ENSMUSG00000064842 ... SIRV7
##   CBFB-MYH11-mcherry
## rowData names(1): GeneLength
## colnames(192): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1
##   SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ...
##   SLX-11312.N712_S508.H5H5YBBXX.s_8.r_1
##   SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
## colData names(0):
## reducedDimNames(0):
## spikeNames(0):
```

We identify the rows corresponding to ERCC spike-in transcripts from the row names.
We store this information in the `SingleCellExperiment` object for future use.
This is necessary as spike-ins require special treatment in downstream steps such as normalization.


```r
isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
summary(isSpike(sce, "ERCC"))
```

```
##    Mode   FALSE    TRUE 
## logical   46611      92
```

This dataset is slightly unusual in that it contains information from another set of spike-in transcripts, the Spike-In RNA Variants (SIRV) set.
For simplicity, we will only use the ERCC spike-ins in this analysis.
Thus, we must remove the rows corresponding to the SIRV transcripts prior to further analysis, which can be done simply by subsetting the `SingleCellExperiment` object.


```r
is.sirv <- grepl("^SIRV", rownames(sce))
sce <- sce[!is.sirv,] 
summary(is.sirv)
```

```
##    Mode   FALSE    TRUE 
## logical   46696       7
```

**Comments from Aaron:**

- Some feature-counting tools will report mapping statistics in the count matrix (e.g., the number of unaligned or unassigned reads).
While these values can be useful for quality control, they would be misleading if treated as gene expression values.
Thus, they should be removed (or at least moved to the `colData`) prior to further analyses.
- Be aware of using the `^ERCC` regular expression for human data where the row names of the count matrix are gene symbols.
An ERCC gene family actually exists in human annotation, so this would result in incorrect identification of genes as spike-in transcripts.
This problem can be avoided by publishing count matrices with standard identifiers (e.g., Ensembl, Entrez).

## Incorporating cell-based annotation

We load in the metadata for each library/cell from the `sdrf.txt` file.
It is important to check that the rows of the metadata table are in the same order as the columns of the count matrix.
Otherwise, incorrect metadata will be assigned to each cell.


```r
metadata <- read.delim(lun.sdrf, check.names=FALSE, header=TRUE)
m <- match(colnames(sce), metadata[["Source Name"]]) # Enforcing identical order.
stopifnot(all(!is.na(m))) # Checking that nothing's missing.
metadata <- metadata[m,]
head(colnames(metadata))
```

```
## [1] "Source Name"                "Comment[ENA_SAMPLE]"       
## [3] "Comment[BioSD_SAMPLE]"      "Characteristics[organism]" 
## [5] "Characteristics[cell line]" "Characteristics[cell type]"
```

We only retain relevant metadata fields to avoid storing unnecessary information in the `colData` of the `SingleCellExperiment` object.
In particular, we keep the plate of origin (i.e., `block`) and phenotype of each cell.
The second field is relevant as all of the cells contain a _CBFB-MYH11_ oncogene, but the expression of this oncogene is only induced in a subset of the cells.


```r
colData(sce)$Plate <- factor(metadata[["Factor Value[block]"]])
pheno <- metadata[["Factor Value[phenotype]"]]
levels(pheno) <- c("induced", "control")
colData(sce)$Oncogene <- pheno
table(colData(sce)$Oncogene, colData(sce)$Plate)
```

```
##          
##           20160113 20160325
##   induced       48       48
##   control       48       48
```

## Incorporating gene-based annotation

Feature-counting tools typically report genes in terms of standard identifiers from Ensembl or Entrez.
These identifiers are used as they are unambiguous and highly stable.
However, they are difficult to interpret compared to the gene symbols which are more commonly used in the literature.
Given the Ensembl identifiers, we obtain the corresponding gene symbols using annotation packages like *[org.Mm.eg.db](https://bioconductor.org/packages/3.9/org.Mm.eg.db)*.


```r
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))
```

```
## DataFrame with 6 rows and 3 columns
##                    GeneLength            ENSEMBL      SYMBOL
##                     <integer>        <character> <character>
## ENSMUSG00000102693       1070 ENSMUSG00000102693          NA
## ENSMUSG00000064842        110 ENSMUSG00000064842          NA
## ENSMUSG00000051951       6094 ENSMUSG00000051951        Xkr4
## ENSMUSG00000102851        480 ENSMUSG00000102851          NA
## ENSMUSG00000103377       2819 ENSMUSG00000103377          NA
## ENSMUSG00000104017       2233 ENSMUSG00000104017          NA
```

It is often desirable to rename the row names of `sce` to the gene symbols, as these are easier to interpret.
However, this requires some work to account for missing and duplicate symbols.
The code below will replace missing symbols with the Ensembl identifier and concatenate duplicated symbols with the (unique) Ensembl identifiers.


```r
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
head(rownames(sce))
```

```
## [1] "ENSMUSG00000102693" "ENSMUSG00000064842" "Xkr4"              
## [4] "ENSMUSG00000102851" "ENSMUSG00000103377" "ENSMUSG00000104017"
```

We also determine the chromosomal location for each gene using the *[TxDb.Mmusculus.UCSC.mm10.ensGene](https://bioconductor.org/packages/3.9/TxDb.Mmusculus.UCSC.mm10.ensGene)* package.
This will be useful later as several quality control metrics will be computed from rows corresponding to mitochondrial genes.


```r
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
    column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")
```

```
##    Mode   FALSE    TRUE    NA's 
## logical   22428      13   24255
```

Alternatively, annotation from BioMart resources can be directly added to the object using the `getBMFeatureAnnos` function from *[scater](https://bioconductor.org/packages/3.9/scater)*.
This may be more convenient than the approach shown above, but depends on an available internet connection to the BioMart databases.

# Quality control on the cells 

## Defining the quality control metrics

Low-quality cells need to be removed to ensure that technical effects do not distort downstream analysis results.
We use several quality control (QC) metrics:

- The library size is defined as the total sum of counts across all features, i.e., genes and spike-in transcripts.
Cells with small library sizes are of low quality as the RNA has not been efficiently captured (i.e., converted into cDNA and amplified) during library preparation.
- The number of expressed features in each cell is defined as the number of features with non-zero counts for that cell.
Any cell with very few expressed genes is likely to be of poor quality as the diverse transcript population has not been successfully captured.
- The proportion of reads mapped to spike-in transcripts is calculated relative to the library size for each cell.
High proportions are indicative of poor-quality cells, where endogenous RNA has been lost during processing (e.g., due to cell lysis or RNA degradation).
The same amount of spike-in RNA to each cell, so an enrichment in spike-in counts is symptomatic of loss of endogenous RNA.
- In the absence of spike-in transcripts, the proportion of reads mapped to genes in the mitochondrial genome can also be used.
High proportions are indicative of poor-quality cells [@islam2014quantitative;@ilicic2016classification], possibly because of loss of cytoplasmic RNA from perforated cells.
The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.

For each cell, we calculate these quality control metrics using the `calculateQCMetrics` function from the *[scater](https://bioconductor.org/packages/3.9/scater)* package [@mccarthy2017scater].
These are stored in the row- and column-wise metadata of the `SingleCellExperiment` for future reference.


```r
mito <- which(rowData(sce)$CHR=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
head(colnames(colData(sce)), 10)
```

```
##  [1] "Plate"                          "Oncogene"                      
##  [3] "is_cell_control"                "total_features_by_counts"      
##  [5] "log10_total_features_by_counts" "total_counts"                  
##  [7] "log10_total_counts"             "pct_counts_in_top_50_features" 
##  [9] "pct_counts_in_top_100_features" "pct_counts_in_top_200_features"
```

The distributions of these metrics are shown in Figure \@ref(fig:qcplot416b), stratified by oncogene induction status and plate of origin.
The aim is to remove putative low-quality cells that have low library sizes, low numbers of expressed features, and high spike-in (or mitochondrial) proportions.
Such cells can interfere with downstream analyses, e.g., by forming distinct clusters that complicate interpretation of the results.


```r
sce$PlateOnco <- paste0(sce$Oncogene, ".", sce$Plate)
multiplot(
    plotColData(sce, y="total_counts", x="PlateOnco"),
    plotColData(sce, y="total_features_by_counts", x="PlateOnco"),
    plotColData(sce, y="pct_counts_ERCC", x="PlateOnco"),
    plotColData(sce, y="pct_counts_Mt", x="PlateOnco"),
    cols=2)
```

<div class="figure">
<img src="reads_files/figure-html/qcplot416b-1.png" alt="Distributions of various QC metrics for all cells in the 416B dataset. This includes the library sizes, number of expressed genes, and proportion of reads mapped to spike-in transcripts or mitochondrial genes." width="100%"  class="widefigure" />
<p class="caption">(\#fig:qcplot416b)Distributions of various QC metrics for all cells in the 416B dataset. This includes the library sizes, number of expressed genes, and proportion of reads mapped to spike-in transcripts or mitochondrial genes.</p>
</div>

It is also valuable to examine how the QC metrics behave with respect to each other (Figure \@ref(fig:qcbiplot416b)).
Generally, they will be in rough agreement, i.e., cells with low total counts will also have low numbers of expressed features and high ERCC/mitochondrial proportions.
Clear discrepancies may correspond to technical differences between batches of cells (see below) or genuine biological differences in RNA content.


```r
par(mfrow=c(1,3))
plot(sce$total_features_by_counts, sce$total_counts/1e6, xlab="Number of expressed genes",
    ylab="Library size (millions)")
plot(sce$total_features_by_counts, sce$pct_counts_ERCC, xlab="Number of expressed genes",
    ylab="ERCC proportion (%)")
plot(sce$total_features_by_counts, sce$pct_counts_Mt, xlab="Number of expressed genes",
    ylab="Mitochondrial proportion (%)")
```

<div class="figure">
<img src="reads_files/figure-html/qcbiplot416b-1.png" alt="Behaviour of each QC metric compared to the total number of expressed features. Each point represents a cell in the 416B dataset." width="960" />
<p class="caption">(\#fig:qcbiplot416b)Behaviour of each QC metric compared to the total number of expressed features. Each point represents a cell in the 416B dataset.</p>
</div>

## Identifying outliers for each metric 

Picking a threshold for these metrics is not straightforward as their absolute values depend on the experimental protocol.
For example, sequencing to greater depth will lead to more reads and more expressed features, regardless of the quality of the cells.
Similarly, using more spike-in RNA in the protocol will result in higher spike-in proportions.
To obtain an adaptive threshold, we assume that most of the dataset consists of high-quality cells, and identify cells that are outliers for the various QC metrics.

Outliers are defined based on the median absolute deviation (MADs) from the median value of each metric across all cells.
We remove cells with log-library sizes that are more than 3 MADs below the median log-library size.
A log-transformation improves resolution at small values, especially when the MAD of the raw values is comparable to or greater than the median.
We also remove cells where the log-transformed number of expressed genes is 3 MADs below the median value.


```r
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce$PlateOnco)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce$PlateOnco)
```

The `batch=` argument ensures that outliers are identified _within_ each level of the specified plate/oncogene factor.
This allows `isOutlier()` to accommodate systematic differences in the QC metrics across plates (Figure \@ref(fig:qcplot416b)),
which can arise due to technical differences in processing (e.g., differences in sequencing depth) rather than any changes in quality.
The same reasoning applies to the oncogene induction status, where induced cells may have naturally fewer expressed genes for biological reasons.
Failing to account for these systematic differences would inflate the MAD estimate and compromise the removal of low-quality cells.

We identify outliers for the proportion-based metrics in a similar manner.
Here, no transformation is required as we are identifying large outliers, for which the distinction should be fairly clear on the raw scale.
We do not need to use the mitochondrial proportions as we already have the spike-in proportions (which serve a similar purpose) for this dataset.
This avoids potential issues arising from genuine differences in mitochondrial content between cell types that may confound outlier identification.


```r
spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher",
    batch=sce$PlateOnco)
```

Subsetting by column will retain only the high-quality cells that pass each filter described above.
We examine the number of cells removed by each filter as well as the total number of retained cells.
Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data quality.


```r
keep <- !(libsize.drop | feature.drop | spike.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=sum(keep))
```

```
##   ByLibSize ByFeature BySpike Remaining
## 1         5         4       6       183
```

We then subset the `SingleCellExperiment` object to retain only the putative high-quality cells.
We also save the original object to file for later use.


```r
sce$PassQC <- keep
saveRDS(sce, file="416B_preQC.rds")
sce <- sce[,keep]
dim(sce)
```

```
## [1] 46696   183
```

**Comments from Aaron:**

- See [this section](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/qc.html#assumptions-of-outlier-identification) for a more detailed discussion of the assumptions underlying the outlier-based detection of low-quality cells.
- `isOutlier()` will also return the exact filter thresholds for each metric (within each batch, if `batch=` is specified).
These may be useful for checking whether the automatically selected thresholds are appropriate.

    
    ```r
    attr(libsize.drop, "thresholds")
    ```
    
    ```
    ##        control.20160113 control.20160325 induced.20160113 induced.20160325
    ## lower            674264         419448.6         499131.8         461501.1
    ## higher              Inf              Inf              Inf              Inf
    ```
    
    ```r
    attr(spike.drop, "thresholds")
    ```
    
    ```
    ##        control.20160113 control.20160325 induced.20160113 induced.20160325
    ## lower              -Inf             -Inf             -Inf             -Inf
    ## higher          8.99581         8.105749         15.50477         12.71858
    ```

# Classification of cell cycle phase 

We use the prediction method described by @scialdone2015computational to classify cells into cell cycle phases based on the gene expression data.
Using a training dataset, the sign of the difference in expression between two genes was computed for each pair of genes.
Pairs with changes in the sign across cell cycle phases were chosen as markers.
Cells in a test dataset can then be classified into the appropriate phase, based on whether the observed sign for each marker pair is consistent with one phase or another.

This approach is implemented in the `cyclone` function from the *[scran](https://bioconductor.org/packages/3.9/scran)* package.
The package contains a pre-trained set of marker pairs for mouse data, which we can load in the the `readRDS` function.
We use the Ensembl identifiers for each gene in our dataset to match up with the names in the pre-trained set of gene pairs.


```r
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
    package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
```

The `cyclone` result for each cell in the HSC dataset is shown in Figure \@ref(fig:phaseplot416b).
Each cell is assigned a score for each phase, with a higher score corresponding to a higher probability that the cell is in that phase.
We focus on the G1 and G2/M scores as these are the most informative for classification.


```r
plot(assignments$score$G1, assignments$score$G2M, 
    xlab="G1 score", ylab="G2/M score", pch=16)
```

<div class="figure">
<img src="reads_files/figure-html/phaseplot416b-1.png" alt="Cell cycle phase scores from applying the pair-based classifier on the 416B dataset. Each point represents a cell, plotted according to its scores for G1 and G2/M phases." width="100%" />
<p class="caption">(\#fig:phaseplot416b)Cell cycle phase scores from applying the pair-based classifier on the 416B dataset. Each point represents a cell, plotted according to its scores for G1 and G2/M phases.</p>
</div>

Cells are classified as being in G1 phase if the G1 score is above 0.5 and greater than the G2/M score; 
    in G2/M phase if the G2/M score is above 0.5 and greater than the G1 score; 
    and in S phase if neither score is above 0.5.
Here, the vast majority of cells are classified as being in G1 phase.
We save these assignments into the `SingleCellExperiment` object for later use.


```r
sce$phases <- assignments$phases
table(sce$phases)
```

```
## 
##  G1 G2M   S 
##  98  62  23
```

Pre-trained classifiers are available in *[scran](https://bioconductor.org/packages/3.9/scran)* for human and mouse data. 
While the mouse classifier used here was trained on data from embryonic stem cells, it is still accurate for other cell types [@scialdone2015computational].
This may be due to the conservation of the transcriptional program associated with the cell cycle [@bertoli2013control;@conboy2007cell].
The pair-based method is also a non-parametric procedure that is robust to most technical differences between datasets.

__Comments from Aaron:__

- To remove confounding effects due to cell cycle phase, we can filter the cells to only retain those in a particular phase (usually G1) for downstream analysis.
Alternatively, if a non-negligible number of cells are in other phases, we can use the assigned phase as a blocking factor.
This protects against cell cycle effects without discarding information, and will be discussed later in more detail.
- The classifier may not be accurate for data that are substantially different from those used in the training set, e.g., due to the use of a different protocol.
In such cases, users can construct a custom classifier from their own training data using the `sandbag` function.
This will also be necessary for other model organisms where pre-trained classifiers are not available.
- Do not filter out low-abundance genes before applying `cyclone`.
Even if a gene is not expressed in *any* cell, it may still be useful for classification if it is phase-specific.
Its lack of expression relative to other genes will still yield informative pairs, and filtering them out would reduce power.

# Examining gene-level expression metrics

## Inspecting the most highly expressed genes

We examine the identities of the most highly expressed genes (Figure \@ref(fig:topgene416b)).
This should generally be dominated by constitutively expressed transcripts, such as those for ribosomal or mitochondrial proteins.
The presence of other classes of features may be cause for concern if they are not consistent with expected biology.
For example, a top set containing many spike-in transcripts suggests that too much spike-in RNA was added during library preparation, while the absence of ribosomal proteins and/or the presence of their pseudogenes are indicative of suboptimal alignment.


```r
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotHighestExprs(sce, n=50) + fontsize
```

<div class="figure">
<img src="reads_files/figure-html/topgene416b-1.png" alt="Percentage of total counts assigned to the top 50 most highly-abundant features in the 416B dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature." width="100%"  class="widefigure" />
<p class="caption">(\#fig:topgene416b)Percentage of total counts assigned to the top 50 most highly-abundant features in the 416B dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature.</p>
</div>

## Filtering out low-abundance genes

Low-abundance genes are problematic as zero or near-zero counts do not contain much information for reliable statistical inference [@bourgon2010independent].
These genes typically do not provide enough evidence to reject the null hypothesis during testing, yet they still increase the severity of the multiple testing correction.
In addition, the discreteness of the counts may interfere with statistical procedures, e.g., by compromising the accuracy of continuous approximations.
Thus, low-abundance genes are often removed in many RNA-seq analysis pipelines before the application of downstream methods.

The "optimal" choice of filtering strategy depends on the downstream application.
A more aggressive filter is usually required to remove discreteness (e.g., for normalization) compared to that required for removing underpowered tests.
For hypothesis testing, the filter statistic should also be independent of the test statistic under the null hypothesis.
Thus, we (or the relevant function) will filter at each step as needed, rather than applying a single filter for the entire analysis.

Several metrics can be used to define low-abundance genes.
The most obvious is the average count for each gene, computed across all cells in the dataset.
We calculate this using the `calcAverage()` function, which also performs some adjustment for library size differences between cells. 
We typically observe a peak of moderately expressed genes following a plateau of lowly expressed genes (Figure \@ref(fig:abhist416b)).


```r
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
    xlab=expression(Log[10]~"average count"))
```

<div class="figure">
<img src="reads_files/figure-html/abhist416b-1.png" alt="Histogram of log-average counts for all genes in the 416B dataset." width="100%" />
<p class="caption">(\#fig:abhist416b)Histogram of log-average counts for all genes in the 416B dataset.</p>
</div>

A minimum threshold can be applied to this value to filter out genes that are lowly expressed.
The example below demonstrates how we _could_ remove genes with average counts less than 1.
The number of `TRUE` values in `demo.keep` corresponds to the number of retained rows/genes after filtering.


```r
demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)
```

```
##    Mode   FALSE    TRUE 
## logical   33490   13206
```

We also examine the number of cells that express each gene.
This is closely related to the average count for most genes, as expression in many cells will result in a higher average (Figure \@ref(fig:nexprshist416b)).
Genes expressed in very few cells are often uninteresting as they are driven by amplification artifacts (though they may also also arise from rare populations).
We could then remove genes that are expressed in fewer than _n_ cells.


```r
num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
    xlab=expression(Log[10]~"average count"))
```

<div class="figure">
<img src="reads_files/figure-html/nexprshist416b-1.png" alt="The number of cells expressing each gene in the 416B dataset, plotted against the log-average count. Intensity of colour corresponds to the number of genes at any given location." width="100%" />
<p class="caption">(\#fig:nexprshist416b)The number of cells expressing each gene in the 416B dataset, plotted against the log-average count. Intensity of colour corresponds to the number of genes at any given location.</p>
</div>

As mentioned above, we will apply these filters at each step rather than applying them globally by subsetting `sce`.
This ensures that the most appropriate filter is used in each application.
Nonetheless, we remove genes that are not expressed in any cell to reduce computational work in downstream steps. 
Such genes provide no information and would be removed by any filtering strategy.


```r
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)
```

```
##    Mode   FALSE    TRUE 
## logical   22833   23863
```

# Normalization of cell-specific biases

## Using the deconvolution method to deal with zero counts

Read counts are subject to differences in capture efficiency and sequencing depth between cells [@stegle2015computational].
Normalization is required to eliminate these cell-specific biases prior to downstream quantitative analyses.
This is often done by assuming that most genes are not differentially expressed (DE) between cells.
Any systematic difference in count size across the non-DE majority of genes between two cells is assumed to represent bias and is removed by scaling.
More specifically, "size factors" are calculated that represent the extent to which counts should be scaled in each library.

Size factors can be computed with several different approaches, e.g., using the `estimateSizeFactorsFromMatrix` function in the *[DESeq2](https://bioconductor.org/packages/3.9/DESeq2)* package [@anders2010differential;@love2014moderated], or with the `calcNormFactors` function [@robinson2010scaling] in the *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* package.
However, single-cell data can be problematic for these bulk data-based methods due to the dominance of low and zero counts.
To overcome this, we pool counts from many cells to increase the size of the counts for accurate size factor estimation [@lun2016pooling].
Pool-based size factors are then "deconvolved" into cell-based factors for normalization of each cell's expression profile.


```r
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.3464  0.7117  0.9098  1.0000  1.1396  3.5695
```

The size factors are well-correlated with the library sizes for all cells (Figure \@ref(fig:normplot416b)).
This suggests that most of the systematic differences between cells are driven by differences in capture efficiency or sequencing depth.
Any DE between cells would yield a non-linear trend between the total count and size factor, and/or increased scatter around the trend.
We observe some evidence of this after oncogene induction, where the size factors after induction are systematically lower.
This is consistent with composition biases [@robinson2010scaling] introduced by upregulation of genes after induction.


```r
plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
    xlab="Library size (millions)", ylab="Size factor",
    col=c("red", "black")[sce$Oncogene], pch=16)
legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,
    legend=levels(sce$Oncogene))
```

<div class="figure">
<img src="reads_files/figure-html/normplot416b-1.png" alt="Size factors from deconvolution, plotted against library sizes for all cells in the 416B dataset. Axes are shown on a log-scale. Wild-type cells are shown in black and oncogene-induced cells are shown in red." width="100%" />
<p class="caption">(\#fig:normplot416b)Size factors from deconvolution, plotted against library sizes for all cells in the 416B dataset. Axes are shown on a log-scale. Wild-type cells are shown in black and oncogene-induced cells are shown in red.</p>
</div>

__Comments from Aaron:__

- While the deconvolution approach is robust to the high frequency of zeroes in scRNA-seq data, it will eventually fail if too many counts are zero.
This manifests as negative size factors that are obviously nonsensical.
To avoid this, the `computeSumFactors` function will automatically remove low-abundance genes prior to the calculation of size factors.
Genes with a library size-adjusted average count below a specified threshold (`min.mean`) are ignored.
For read count data, the default value of 1 is usually satisfactory.
- Cell-based QC should always be performed prior to normalization, to remove cells with very low numbers of expressed genes.
Otherwise, `computeSumFactors()` may yield negative size factors for low-quality cells.
This is because too many zeroes are present in these cells, reducing the effectiveness of pooling to eliminate zeroes.
See `?computeSumFactors` for more details on how the function resolves negative size factors.
- The `sizes` argument can be used to specify the number of pool sizes to use to compute the size factors.
More `sizes` yields more precise estimates at the cost of some computational time and memory.
In general, all `sizes` should be above 20 cells to ensure that there are sufficient non-zero expression values in each pool.
The total number of cells should also be at least 100 for effective pooling.
- For highly heterogeneous datasets, it is advisable to perform a rough clustering of the cells to weaken the non-DE assumption.
This can be done with the `quickCluster()` function and the results passed to `computeSumFactors()` via the `clusters` argument.
We demonstrate this approach with a larger data set in the next [workflow](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/umis.html#normalization-of-cell-specific-biases).

## Computing separate size factors for spike-in transcripts

Size factors computed from the counts for endogenous genes are usually not appropriate for normalizing the counts for spike-in transcripts.
Consider an experiment without library quantification, i.e., the amount of cDNA from each library is _not_ equalized prior to pooling and multiplexed sequencing.
Here, cells containing more RNA have greater counts for endogenous genes and thus larger size factors to scale down those counts.
However, the same amount of spike-in RNA is added to each cell during library preparation.
This means that the counts for spike-in transcripts are not subject to the effects of RNA content.
Attempting to normalize the spike-in counts with the gene-based size factors will lead to over-normalization and incorrect quantification of expression.
Similar reasoning applies in cases where library quantification is performed. 
For a constant total amount of cDNA, any increases in endogenous RNA content will suppress the coverage of spike-in transcripts.
As a result, the bias in the spike-in counts will be opposite to that captured by the gene-based size factor.

To ensure normalization is performed correctly, we compute a separate set of size factors for the spike-in set.
For each cell, the spike-in-specific size factor is defined as the total count across all transcripts in the spike-in set.
This assumes that none of the spike-in transcripts are differentially expressed, which is reasonable given that the same amount and composition of spike-in RNA should have been added to each cell [@lun2017assessing].
(See below for a more detailed discussion on spike-in normalization.)
These size factors are stored in a separate field of the `SingleCellExperiment` object by setting `general.use=FALSE` in `computeSpikeFactors`.
This ensures that they will only be used with the spike-in transcripts but not the endogenous genes.


```r
sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)
```

## Applying the size factors to normalize gene expression

The count data are used to compute normalized log-expression values for use in downstream analyses.
Each value is defined as the log~2~-ratio of each count to the size factor for the corresponding cell, after adding a prior count of 1 to avoid undefined values at zero counts.
Division by the size factor ensures that any cell-specific biases are removed.
If spike-in-specific size factors are present in `sce`, they will be automatically applied to normalize the spike-in transcripts separately from the endogenous genes. 


```r
sce <- normalize(sce)
```

The log-transformation is useful as it means that any differences in the values directly represent log~2~-fold changes in expression between cells.
This is usually more relevant than the absolute differences in coverage, which need to be interpreted in the context of the overall abundance.
The log-transformation also provides some measure of variance stabilization [@law2014voom], so that high-abundance genes with large variances do not dominate downstream analyses.
The computed values are stored as an `"logcounts"` matrix in addition to the other assay elements.



# Modelling the technical noise in gene expression

Variability in the observed expression values across genes can be driven by genuine biological heterogeneity or uninteresting technical noise. 
To distinguish between these two possibiltiies, we need to model the technical component of the variance of the expression values for each gene.
We do so using the set of spike-in transcripts, which were added in the same quantity to each cell.
Thus, the spike-in transcripts should exhibit no biological variability, i.e., any variance in their counts should be technical in origin.

We use the `trendVar()` function to fit a mean-dependent trend to the variances of the log-expression values for the spike-in transcripts.
We set `block=` to block on the plate of origin for each cell, to ensure that technical differences between plates do not inflate the variances.
Given the mean abundance of a gene, the fitted value of the trend is then used as an estimate of the technical component for that gene.
The biological component of the variance is finally calculated by subtracting the technical component from the total variance of each gene with the `decomposeVar` function.


```r
var.fit <- trendVar(sce, parametric=TRUE, block=sce$Plate,
    loess.args=list(span=0.3))
var.out <- decomposeVar(sce, var.fit)
head(var.out)
```

```
## DataFrame with 6 rows and 6 columns
##                                   mean               total
##                              <numeric>           <numeric>
## ENSMUSG00000103377 0.00807160215928894   0.011921865486065
## ENSMUSG00000103147  0.0346526072192529  0.0722196162535234
## ENSMUSG00000103161 0.00519472222570747 0.00493857699521053
## ENSMUSG00000102331   0.018666093059853   0.032923591860573
## ENSMUSG00000102948   0.059057000132083  0.0881371257735823
## Rp1                 0.0970243712569606    0.45233813529556
##                                    bio               tech            p.value
##                              <numeric>          <numeric>          <numeric>
## ENSMUSG00000103377 -0.0238255786088717 0.0357474440949367                  1
## ENSMUSG00000103147 -0.0812680860584481  0.153487702311972  0.999999999992144
## ENSMUSG00000103161 -0.0180705438722202 0.0230091208674307                  1
## ENSMUSG00000102331 -0.0497487337065681  0.082672325567141  0.999999999998056
## ENSMUSG00000102948  -0.173441452696662  0.261578578470245                  1
## Rp1                 0.0226096722909625  0.429728463004597 0.0354980966384924
##                                  FDR
##                            <numeric>
## ENSMUSG00000103377                 1
## ENSMUSG00000103147                 1
## ENSMUSG00000103161                 1
## ENSMUSG00000102331                 1
## ENSMUSG00000102948                 1
## Rp1                0.153727758280855
```

We visually inspect the trend to confirm that it corresponds to the spike-in variances (Figure \@ref(fig:hvgplot416b))). 
The wave-like shape is typical of the mean-variance trend for log-expression values.
A linear increase in the variance is observed as the mean increases from zero, as larger variances are possible when the counts increase.
At very high abundances, the effect of sampling noise decreases due to the law of large numbers, resulting in a decrease in the variance.


```r
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
```

<div class="figure">
<img src="reads_files/figure-html/hvgplot416b-1.png" alt="Variance of normalized log-expression values for each gene in the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)." width="100%" />
<p class="caption">(\#fig:hvgplot416b)Variance of normalized log-expression values for each gene in the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red).</p>
</div>

We check the distribution of expression values for the genes with the largest biological components.
This ensures that the variance estimate is not driven by one or two outlier cells (Figure \@ref(fig:hvgvioplot416b)).


```r
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize
```

<div class="figure">
<img src="reads_files/figure-html/hvgvioplot416b-1.png" alt="Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the 416B dataset. Each point represents the log-expression value in a single cell." width="100%" />
<p class="caption">(\#fig:hvgvioplot416b)Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the 416B dataset. Each point represents the log-expression value in a single cell.</p>
</div>

**Comments from Aaron:**

- In practice, trend fitting is complicated by the small number of spike-in transcripts and the uneven distribution of their abundances.
See [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/var.html#details-of-trend-fitting-parameters) for more details on how to refine the fit.
- In the absence of spike-ins, users can set `use.spikes=FALSE` to fit a trend to the variances of the endogenous genes (see [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/var.html#wwhen-spike-ins-are-unavailable)).
Alternatively, we can create a trend based on the assumption of Poisson technical noise, as described [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/tenx.html#modelling-the-mean-variance-trend).
- Negative biological components are often obtained from `decomposeVar`. 
These are intuitively meaningless as it is impossible for a gene to have total variance below technical noise.
Nonetheless, such values occur due to imprecise estimation of the total variance, especially for low numbers of cells.
- `decomposeVar` also yields _p_-values that can be used to define HVGs at a specific threshold for the false discovery rate (FDR).
We will discuss this in more detail later, as formal detection of HVGs is not necessary for feature selection during data exploration.

# Removing the batch effect

As previously mentioned, the data were collected on two plates.
Small uncontrollable differences in processing between plates can result in a batch effect, i.e., systematic differences in expression between cells on different plates.
Such differences are not interesting and can be removed by applying the `removeBatchEffect()` function from the *[limma](https://bioconductor.org/packages/3.9/limma)* package [@ritchie2015limma].
This removes the effect of the plate of origin while accounting for the (interesting) effect of oncogene induction.


```r
library(limma)
assay(sce, "corrected") <- removeBatchEffect(logcounts(sce), 
    design=model.matrix(~sce$Oncogene), batch=sce$Plate)
assayNames(sce)
```

```
## [1] "counts"    "logcounts" "corrected"
```

Manual batch correction is necessary for downstream procedures that are not model-based, e.g., clustering and most forms of dimensionality reduction.
However, if an analysis method can accept a design matrix, blocking on nuisance factors in the design matrix is preferable to using `removeBatchEffect()`.
This is because the latter does not account for the loss of residual degrees of freedom, nor the uncertainty of estimation of the blocking factor terms.

**Comments from Aaron:**

- `removeBatchEffect()` performs a linear regression and sets the coefficients corresponding to the blocking factors to zero.
This is effective provided that the population composition within each batch is known (and supplied as `design=`) or identical across batches.
Such an assumption is reasonable for this dataset, involving a homogeneous cell line population on both plates.
However, in most scRNA-seq applications, the factors of variation are not identical across batches and not known in advance.
This motivates the use of more sophisticated batch correction methods such as `mnnCorrect()`.

# Denoising expression values using PCA

Once the technical noise is modelled, we can use principal components analysis (PCA) to remove random technical noise.
Consider that each cell represents a point in the high-dimensional expression space, where the spread of points represents the total variance.
PCA identifies axes in this space that capture as much of this variance as possible.
Each axis is a principal component (PC), where any early PC will explain more of the variance than a later PC.

We assume that biological processes involving co-regulated groups of genes will account for the most variance in the data.
If this is the case, this process should be represented by one or more of the earlier PCs.
In contrast, random technical noise affects each gene independently and will be represented by later PCs.
The `denoisePCA()` function removes later PCs until the total discarded variance is equal to the sum of technical components for all genes used in the PCA.


```r
sce <- denoisePCA(sce, technical=var.out, assay.type="corrected")
dim(reducedDim(sce, "PCA")) 
```

```
## [1] 183  24
```

The function returns a `SingleCellExperiment` object containing the PC scores for each cell in the `reducedDims` slot.
The aim is to eliminate technical noise and enrich for biological signal in the retained PCs.
This improves resolution of the underlying biology during downstream procedures such as clustering.

__Comments from Aaron:__

- `denoisePCA()` will only use genes that have positive biological components, i.e., variances greater than the fitted trend.
This guarantees that the total technical variance to be discarded will not be greater than the total variance in the data.
- For the `technical=` argument, the function will also accept the trend function directly (i.e., `var.fit$trend`) or a vector of technical components per gene.
Here, we supply the `DataFrame` from `decomposeVar()` to allow the function to adjust for the loss of residual degrees of freedom after batch correction.
Specifically, the variance in the batch-corrected matrix is slightly understated, requiring some rescaling of the technical components to compensate.
- No filtering is performed on abundance here, which ensures that PCs corresponding to rare subpopulations can still be detected. 
Discreteness is less of an issue as low-abundance genes also have lower variance, thus reducing their contribution to the PCA.
- It is also possible to obtain a low-rank approximation of the original expression matrix, capturing the variance equivalent to the retained PCs.
This is useful for denoising prior to downstream procedures that require gene-wise expression values.


```r
sce2 <- denoisePCA(sce, technical=var.fit$trend, 
    assay.type="corrected", value="lowrank") 
assayNames(sce2)
```

```
## [1] "counts"    "logcounts" "corrected" "lowrank"
```



# Visualizing data in low-dimensional space

## With PCA

We visualize the relationships between cells by constructing pairwise PCA plots for the first three components (Figure \@ref(fig:pcaplot416b-onco)).
Cells with similar expression profiles should be located close together in the plot, while dissimilar cells should be far apart.
In this case, we observe a clear separation of cells based on the oncogene induction status, consistent with the expected effects on the transcriptome.


```r
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
    colour_by="Oncogene") + fontsize
```

<div class="figure">
<img src="reads_files/figure-html/pcaplot416b-onco-1.png" alt="Pairwise PCA plots of the first three PCs in the 416B dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured according to oncogene induction status." width="864" />
<p class="caption">(\#fig:pcaplot416b-onco)Pairwise PCA plots of the first three PCs in the 416B dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured according to oncogene induction status.</p>
</div>

By comparison, we observe no clear separation of cells by batch (Figure \@ref(fig:pcaplot416b-batch)).
This indicates that our batch correction step using `removeBatchEffect()` was successful.


```r
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
    colour_by="Plate") + fontsize
```

<div class="figure">
<img src="reads_files/figure-html/pcaplot416b-batch-1.png" alt="Pairwise PCA plots of the first three PCs in the 416B dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured according to the plate of origin." width="864" />
<p class="caption">(\#fig:pcaplot416b-batch)Pairwise PCA plots of the first three PCs in the 416B dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured according to the plate of origin.</p>
</div>

Note that `plotReducedDim()` will use the PCA results that were already stored in `sce` by `denoisePCA()`.
This allows us to rapidly generate new plots with different aesthetics, without repeating the entire PCA computation. 
Similarly, `plotPCA()` will use existing results if they are available in the `SingleCellExperiment`, and will recalculate them otherwise.
Users should set `rerun=TRUE` to forcibly recalculate the PCs in the presence of existing results.

__Comments from Aaron:__

- For each visualization method, additional cell-specific information can be incorporated into the size or shape of each point.
This is done using the `size_by=` and `shape_by=` arguments in most plotting functions.
- More components can be shown but these are usually less informative as they explain less of the variance. 
They are also often more difficult to interpret as they are defined to be orthogonal to earlier PCs (and thus dependent on what is detected in those PCs).

## With _t_-SNE

Another widely used approach for dimensionality reduction is the _t_-stochastic neighbour embedding (_t_-SNE) method [@van2008visualizing].
_t_-SNE tends to work better than PCA for separating cells in more diverse populations.
This is because the former can directly capture non-linear relationships in high-dimensional space, whereas the latter must represent them on linear axes.
However, this improvement comes at the cost of more computational effort and requires the user to consider parameters such as the random seed and perplexity (see comments).

We demonstrate the generation of _t_-SNE plots in Figure \@ref(fig:tsneplot416b) using the `plotTSNE()` function.
We set `use_dimred="PCA"` to perform the _t_-SNE on the low-rank approximation of the data, allowing the algorithm to take advantage of the previous denoising step.


```r
set.seed(100)
out5 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=5),
    colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 5")

set.seed(100)
out10 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=10),
    colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 10")

set.seed(100)
out20 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=20),
    colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 20")

multiplot(out5, out10, out20, cols=3)
```

<div class="figure">
<img src="reads_files/figure-html/tsneplot416b-1.png" alt="_t_-SNE plots constructed from the denoised PCs in the 416B dataset, using a range of perplexity values. Each point represents a cell, coloured according to its oncogene induction status. Bars represent the coordinates of the cells on each axis." width="1440" />
<p class="caption">(\#fig:tsneplot416b)_t_-SNE plots constructed from the denoised PCs in the 416B dataset, using a range of perplexity values. Each point represents a cell, coloured according to its oncogene induction status. Bars represent the coordinates of the cells on each axis.</p>
</div>

_t_-SNE is a stochastic method, so users should run the algorithm several times to ensure that the results are representative.
Scripts should set a seed (via the `set.seed()` command) to ensure that the chosen results are reproducible.
It is also advisable to test different settings of the "perplexity" parameter as this will affect the distribution of points in the low-dimensional space.

Here, we call `runTSNE()` with a perplexity of 20 to store the _t_-SNE results inside our `SingleCellExperiment` object.
This avoids repeating the calculations whenever we want to create a new plot with `plotTSNE()`, as the stored results will be used instead.
Again, users can set `rerun=TRUE` to force recalculation in the presence of stored results.


```r
set.seed(100)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=20)
reducedDimNames(sce)
```

```
## [1] "PCA"  "TSNE"
```

There are many other dimensionality reduction techniques that we do not consider here but could also be used, e.g., multidimensional scaling, diffusion maps.
These have their own advantages and disadvantages -- for example, diffusion maps (see `plotDiffusionMap`) place cells along a continuous trajectory and are suited for visualizing graduated processes like differentiation [@angerer2016destiny].

__Comments from Aaron:__

- A good guide on how to interpret _t_-SNE plots can be found at http://distill.pub/2016/misread-tsne/.
This demonstrates how distances between clusters in the 2-dimensional embedding have little meaning, as does the apparent "size" (i.e., spread) of the clusters.

# Clustering cells into putative subpopulations

## Defining cell clusters from expression data

The denoised log-expression values are used to cluster cells into putative subpopulations.
Specifically, we perform hierarchical clustering on the Euclidean distances between cells, using Ward's criterion to minimize the total variance within each cluster.
This yields a dendrogram that groups together cells with similar expression patterns across the chosen genes.


```r
pcs <- reducedDim(sce, "PCA")
my.dist <- dist(pcs)
my.tree <- hclust(my.dist, method="ward.D2")
```

Clusters are explicitly defined by applying a dynamic tree cut [@langfelder2008defining] to the dendrogram.
This exploits the shape of the branches in the dendrogram to refine the cluster definitions, and is more appropriate than `cutree` for complex dendrograms.
Greater control of the empirical clusters can be obtained by manually specifying `cutHeight` in `cutreeDynamic`.
We also set `minClusterSize` to a lower value than the default of 20, to avoid spurious aggregation of distant small clusters.


```r
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), 
    minClusterSize=10, verbose=0))
```

We examine the distribution of cells in each cluster with respect to known factors.
Each cluster is comprised of cells from both batches, indicating that the clustering is not driven by a batch effect.
Differences in the composition of each cluster are observed with respect to `Oncogene`, consistent with a biological effect of oncogene induction.


```r
table(my.clusters, sce$Plate)
```

```
##            
## my.clusters 20160113 20160325
##           1       41       39
##           2       19       20
##           3       16       11
##           4       10       14
##           5        5        8
```

```r
table(my.clusters, sce$Oncogene)
```

```
##            
## my.clusters induced control
##           1      80       0
##           2       0      39
##           3       0      27
##           4       0      24
##           5      13       0
```

We visualize the cluster assignments for all cells on the _t_-SNE plot in Figure \@ref(fig:tsnecluster416b).
Adjacent cells are generally assigned to the same cluster, indicating that the clustering procedure was applied correctly.


```r
sce$cluster <- factor(my.clusters)
plotTSNE(sce, colour_by="cluster") + fontsize
```

<div class="figure">
<img src="reads_files/figure-html/tsnecluster416b-1.png" alt="_t_-SNE plot of the denoised PCs of the 416B dataset. Each point represents a cell and is coloured according to the cluster identity to which it was assigned." width="100%" />
<p class="caption">(\#fig:tsnecluster416b)_t_-SNE plot of the denoised PCs of the 416B dataset. Each point represents a cell and is coloured according to the cluster identity to which it was assigned.</p>
</div>

We check the separatedness of the clusters using the silhouette width (Figure \@ref(fig:silhouette416b)).
Cells with large positive silhouette widths are closer to other cells in the _same_ cluster than to cells in _different_ clusters.
Conversely, cells with negative widths are closer to other clusters than to other cells in the cluster to which it was assigned.
Each cluster would ideally contain many cells with large positive widths, indicating that it is well-separated from other clusters.


```r
library(cluster)
clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
sil <- silhouette(my.clusters, dist = my.dist)
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
    border=sil.cols, col=sil.cols, do.col.sort=FALSE) 
```

<div class="figure">
<img src="reads_files/figure-html/silhouette416b-1.png" alt="Barplot of silhouette widths for cells in each cluster. Each cluster is assigned a colour and cells with positive widths are coloured according to the colour of its assigned cluster. Any cell with a negative width is coloured according to the colour of the cluster that it is closest to. The average width for all cells in each cluster is shown, along with the average width for all cells in the dataset." width="100%" />
<p class="caption">(\#fig:silhouette416b)Barplot of silhouette widths for cells in each cluster. Each cluster is assigned a colour and cells with positive widths are coloured according to the colour of its assigned cluster. Any cell with a negative width is coloured according to the colour of the cluster that it is closest to. The average width for all cells in each cluster is shown, along with the average width for all cells in the dataset.</p>
</div>

The silhouette width can be used to determine the parameter values that maximize the separation between clusters.
For example, we could vary the cut height or splitting depth in `cutreeDynamic` to maximize the average silhouette width across all cells.
This usually provides a satisfactory initial clustering for further examination.
However, keep in mind that the granularity of clustering is much like the magnification on a microscope.
Different views of the data can be obtained with different granularities, some of which may be suboptimal on measures of separation.
Users should not fixate on the clustering with the greatest separation if it does not provide the desired granularity for a particular biological question.

Most cells have relatively small silhouette positive widths in Figure \@ref(fig:silhouette416b), indicating that the separation between clusters is weak.
This may be symptomatic of over-clustering where clusters that are clearly defined on oncogene induction status are further split into subsets that are less well separated.
Nonetheless, we will proceed with the current clustering scheme in `my.clusters`, as it provides reasonable partitions for further characterization of heterogeneity.



**Comments from Aaron:**

- An alternative clustering strategy is to use a matrix of distances derived from correlations (e.g., as in `quickCluster`).
This is more robust to noise and normalization errors, but is also less sensitive to subtle changes in the expression profiles. 
- Both Ward's criterion and complete linkage yield spherical, compact clusters.
In particular, complete linkage favours the formation of clusters with the same diameter.
This may be desirable in some cases but is less appropriate when subpopulations differ in their variance.
Thus, we typically use Ward's criterion for our initial clustering. 
Of course, it is simple (and recommended) to try other approaches provided that some assessment is performed, e.g., using the silhouette width.

## Detecting marker genes between clusters

Once putative subpopulations are identified by clustering, we can identify marker genes for each cluster using the `findMarkers` function.
This performs Welch $t$-tests on the log-expression values for every gene and between every pair of clusters [@soneson2018bias].
The aim is to test for DE in each cluster compared to the others while blocking on uninteresting factors such as the plate of origin (see [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/de.html#blocking-on-uninteresting-factors-of-variation) for details).
The top DE genes are likely to be good candidate markers as they can effectively distinguish between cells in different clusters.


```r
markers <- findMarkers(sce, my.clusters, block=sce$Plate)
```

For each cluster, the DE results of the relevant comparisons are consolidated into a single output table.
This allows a set of marker genes to be easily defined by taking the top DE genes from each pairwise comparison between clusters.
For example, to construct a marker set for cluster 1 from the top 10 genes of each comparison, one would filter `marker.set` to retain rows with `Top` less than or equal to 10.
Other statistics are also reported for each gene, including the adjusted $p$-values (see below) and the log-fold changes relative to every other cluster.




```r
marker.set <- markers[["1"]]
head(marker.set, 10)
```

```
## DataFrame with 10 rows and 7 columns
##              Top              p.value                  FDR           logFC.2
##        <integer>            <numeric>            <numeric>         <numeric>
## Aurkb          1 6.65863261824492e-75 1.58362259559721e-70 -7.37163350569914
## Tk1            1 6.41442387689833e-64 3.81385607660686e-60 -4.92754579848056
## Myh11          1 4.95086465211539e-49 9.81220116843835e-46  4.42159318376489
## Cdca8          1 2.22334940769686e-46  3.5251945975503e-43 -6.84273526334783
## Ccna2          2 1.16841222330534e-68 1.38941739534356e-64  -7.3079756458424
## Rrm2           2 1.48873842020081e-56 5.05809512109088e-53 -5.52120322191947
## Cks1b          2 3.83636139889977e-39 2.40105745131667e-36  -6.6792118199827
## Pirb           2 1.83893803490065e-34 6.15992440620318e-32  5.25803749166673
## Pimreg         3 7.41004737119548e-68 5.87443855430467e-64 -7.30454210816126
## Pclaf          3  8.9651722100101e-51 2.13218690670669e-47 -5.60087985621813
##                  logFC.3            logFC.4            logFC.5
##                <numeric>          <numeric>          <numeric>
## Aurkb  -6.72799345321135  -1.95039440976238  -6.91128802096519
## Tk1    -7.74065113926215  -3.53749565362853  -4.63516273649457
## Myh11   4.30812918035861   4.45235717737968    1.0413149433198
## Cdca8  -4.88595092732349  -2.43821402084038  -7.12791471326961
## Ccna2   -6.9676852052539  -2.46589325823089  -7.12692721843331
## Rrm2   -7.94685699818751  -3.19173143688883  -5.42878091964762
## Cks1b  -5.92137181826573  -4.37146346518812    -6.214592473138
## Pirb    5.18195596569259   5.87631057587633 0.0704964218555272
## Pimreg -5.91099762335684 -0.874660676141792  -7.01798853404623
## Pclaf  -7.56997893312346  -2.36631043435985  -5.16956927698937
```



We save the list of candidate marker genes for further examination.


```r
write.table(marker.set, file="416B_marker_1.tsv", sep="\t", 
    quote=FALSE, col.names=NA)
```

We visualize the expression profiles of the top candidates to verify that the DE signature is robust (Figure \@ref(fig:heatmapmarker416b)). 
Most of the top markers have strong and consistent up- or downregulation in cells of cluster 1 compared to some or all of the other clusters.
A cursory examination of the heatmap indicates that cluster 1 contains oncogene-induced cells with strong downregulation of DNA replication and cell cycle genes.
This is consistent with the potential induction of senescence as an anti-tumorigenic response [@wajapeyee2010senescence].
A more comprehensive investigation of the function of these markers can be performed with gene set enrichment analyses, e.g., using `kegga` or `goana` from *[limma](https://bioconductor.org/packages/3.9/limma)*.


```r
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=top.markers, columns=order(sce$cluster), 
    colour_columns_by=c("cluster", "Plate", "Oncogene"),
    cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5)) 
```

<div class="figure">
<img src="reads_files/figure-html/heatmapmarker416b-1.png" alt="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 1 in the 416B dataset. Column colours represent the cluster to which each cell is assigned, the plate of origin or the oncogene induction status of each cell, as indicated by the legend." width="960" />
<p class="caption">(\#fig:heatmapmarker416b)Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 1 in the 416B dataset. Column colours represent the cluster to which each cell is assigned, the plate of origin or the oncogene induction status of each cell, as indicated by the legend.</p>
</div>

Many of the markers in Figure \@ref(fig:heatmapmarker416b) are not uniquely up- or downregulated in the chosen cluster.
Testing for unique DE tends to be too stringent as it overlooks important genes that are expressed in two or more clusters.
For example, in a mixed population of CD4^+^-only, CD8^+^-only, double-positive and double-negative T cells, neither _Cd4_ or _Cd8_ would be detected as subpopulation-specific markers because each gene is expressed in two subpopulations.
With our approach, both of these genes will be picked up as candidate markers as they will be DE between at least one pair of subpopulations.
A combination of markers can then be chosen to characterize a subpopulation, which is more flexible than trying to find uniquely DE genes.

We strongly recommend selecting some markers for use in validation studies with an independent replicate population of cells.
The aim is to identify a corresponding subset of cells that express the upregulated markers and do not express the downregulated markers.
Ideally, a different technique for quantifying expression would also be used during validation, e.g., fluorescent _in situ_ hybridisation or quantitative PCR.
This confirms that the subpopulation genuinely exists and is not an artifact of scRNA-seq or the computational analysis.

__Comments from Aaron:__

- By setting `direction="up"`, `findMarkers` will only return genes that are upregulated in each cluster compared to the others.
This is convenient in highly heterogeneous populations to focus on genes that can immediately identify each cluster.
While lack of expression may also be informative, it is less useful for positive identification.
- `findMarkers` can also be directed to find genes that are DE between the chosen cluster and _all_ other clusters.
This should be done by setting `pval.type="all"`, which defines the p-value for each gene as the maximum value across all pairwise comparisons involving the chosen cluster.
Combined with `direction="up"`, this can be used to identify unique markers for each cluster.
However, this is sensitive to overclustering, as unique marker genes will no longer exist if a cluster is split into two smaller subclusters.
- It must be stressed that the (adjusted) _p_-values computed here cannot be properly interpreted as measures of significance.
This is because the clusters have been empirically identified from the data, see [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/de.html#misinterpretation-of-de-$p$-values).

# Concluding remarks

Once the basic analysis is completed, it is often useful to save the `SingleCellExperiment` object to file with the `saveRDS` function.
The object can then be easily restored into new R sessions using the `readRDS` function.
This allows further work to be conducted without having to repeat all of the processing steps described above.


```r
saveRDS(file="416B_data.rds", sce)
```

A variety of methods are available to perform more complex analyses on the processed expression data.
For example, cells can be ordered in pseudotime (e.g., for progress along a differentiation pathway) with *[monocle](https://bioconductor.org/packages/3.9/monocle)* [@trapnell2014dynamics] or *[TSCAN](https://bioconductor.org/packages/3.9/TSCAN)* [@ji2016tscan]; 
cell-state hierarchies can be characterized with the *[sincell](https://bioconductor.org/packages/3.9/sincell)* package [@julia2015sincell];
and oscillatory behaviour can be identified using *[Oscope](https://bioconductor.org/packages/3.9/Oscope)* [@leng2015oscope].
HVGs can be used in gene set enrichment analyses to identify biological pathways and processes with heterogeneous activity, using packages designed for bulk data like *[topGO](https://bioconductor.org/packages/3.9/topGO)* or with dedicated single-cell methods like *[scde](https://bioconductor.org/packages/3.9/scde)* [@fan2016characterizing].
Full descriptions of these analyses are outside the scope of this workflow, so interested users are advised to consult the relevant documentation.



All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org).
The specific version numbers of the packages used are shown below, along with the version of the R installation.


```r
sessionInfo()
```

```
## R Under development (unstable) (2019-02-19 r76128)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.5 LTS
## 
## Matrix products: default
## BLAS: /home/cri.camres.org/lun01/Software/R/trunk/lib/libRblas.so
## LAPACK: /home/cri.camres.org/lun01/Software/R/trunk/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] cluster_2.0.7-1                       
##  [2] dynamicTreeCut_1.63-1                 
##  [3] limma_3.39.12                         
##  [4] scran_1.11.20                         
##  [5] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0
##  [6] GenomicFeatures_1.35.7                
##  [7] scater_1.11.11                        
##  [8] ggplot2_3.1.0                         
##  [9] org.Mm.eg.db_3.7.0                    
## [10] AnnotationDbi_1.45.0                  
## [11] SingleCellExperiment_1.5.2            
## [12] SummarizedExperiment_1.13.0           
## [13] DelayedArray_0.9.8                    
## [14] BiocParallel_1.17.15                  
## [15] matrixStats_0.54.0                    
## [16] Biobase_2.43.1                        
## [17] GenomicRanges_1.35.1                  
## [18] GenomeInfoDb_1.19.2                   
## [19] IRanges_2.17.4                        
## [20] S4Vectors_0.21.10                     
## [21] BiocGenerics_0.29.1                   
## [22] BiocFileCache_1.7.0                   
## [23] dbplyr_1.3.0                          
## [24] knitr_1.21                            
## [25] BiocStyle_2.11.0                      
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] RColorBrewer_1.1-2       progress_1.2.0          
##  [5] httr_1.4.0               tools_3.6.0             
##  [7] R6_2.4.0                 irlba_2.3.3             
##  [9] KernSmooth_2.23-15       vipor_0.4.5             
## [11] DBI_1.0.0                lazyeval_0.2.1          
## [13] colorspace_1.4-0         withr_2.1.2             
## [15] processx_3.2.1           tidyselect_0.2.5        
## [17] gridExtra_2.3            prettyunits_1.0.2       
## [19] bit_1.1-14               curl_3.3                
## [21] compiler_3.6.0           BiocNeighbors_1.1.12    
## [23] rtracklayer_1.43.1       labeling_0.3            
## [25] bookdown_0.9             scales_1.0.0            
## [27] callr_3.1.1              rappdirs_0.3.1          
## [29] stringr_1.4.0            digest_0.6.18           
## [31] Rsamtools_1.99.2         rmarkdown_1.11          
## [33] XVector_0.23.0           pkgconfig_2.0.2         
## [35] htmltools_0.3.6          highr_0.7               
## [37] rlang_0.3.1              RSQLite_2.1.1           
## [39] DelayedMatrixStats_1.5.2 dplyr_0.8.0.1           
## [41] RCurl_1.95-4.11          magrittr_1.5            
## [43] BiocSingular_0.99.12     simpleSingleCell_1.7.17 
## [45] GenomeInfoDbData_1.2.0   Matrix_1.2-16           
## [47] Rcpp_1.0.0               ggbeeswarm_0.6.0        
## [49] munsell_0.5.0            viridis_0.5.1           
## [51] edgeR_3.25.3             stringi_1.3.1           
## [53] yaml_2.2.0               zlibbioc_1.29.0         
## [55] Rtsne_0.15               plyr_1.8.4              
## [57] grid_3.6.0               blob_1.1.1              
## [59] crayon_1.3.4             lattice_0.20-38         
## [61] Biostrings_2.51.2        cowplot_0.9.4           
## [63] hms_0.4.2                locfit_1.5-9.1          
## [65] ps_1.3.0                 pillar_1.3.1            
## [67] igraph_1.2.4             reshape2_1.4.3          
## [69] codetools_0.2-16         biomaRt_2.39.2          
## [71] XML_3.98-1.17            glue_1.3.0              
## [73] evaluate_0.13            BiocManager_1.30.4      
## [75] gtable_0.2.0             purrr_0.3.0             
## [77] assertthat_0.2.0         xfun_0.5                
## [79] rsvd_1.0.0               viridisLite_0.3.0       
## [81] pheatmap_1.0.12          tibble_2.0.1            
## [83] GenomicAlignments_1.19.1 beeswarm_0.2.3          
## [85] memoise_1.1.0            statmod_1.4.30
```

# References

