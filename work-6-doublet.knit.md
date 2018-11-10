---
title: Detecting doublets in single-cell RNA-seq data
author:
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
date: "2018-11-10"
vignette: >
  %\VignetteIndexEntry{08. Detecting doublets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---



# Overview

In single-cell RNA sequencing (scRNA-seq) experiments, doublets are artifactual libraries generated from two cells.
They typically arise due to errors in cell sorting or capture, especially in droplet-based protocols [@zheng2017massively] involving thousands of cells.
Doublets are obviously undesirable when the aim is to characterize populations at the _single_-cell level.
In particular, they can incorrectly suggest the existence of intermediate populations or transitory states that not actually exist.
Thus, it is desirable to remove doublet libraries so that they do not compromise interpretation of the results.

Several experimental strategies are available for doublet removal.
One approach exploits natural genetic variation when pooling cells from multiple donor individuals [@kang2018multiplexed].
Doublets can be identified as libraries with allele combinations that do not exist in any single donor.
Another approach is to mark a subset of cells (e.g., all cells from one sample) with an antibody conjugated to a different oligonucleotide [@stoeckius2017hashing].
Upon pooling, libraries that are observed to have different oligonucleotides are considered to be doublets and removed.
These approaches can be highly effective but rely on experimental information that may not be available.

A more general approach is to infer doublets from the expression profiles alone [@dahlin2018single].
In this workflow, we will describe two purely computational approaches for detecting doublets from scRNA-seq data.
The main difference between These two methods is whether or not they need cluster information beforehand.
Both are implemented in the *[scran](https://bioconductor.org/packages/3.9/scran)* package from the open-source Bioconductor project [@huber2015orchestrating].
We will demonstrate the use of these methods on data from a droplet-based scRNA-seq study of the mouse mammary gland [@bach2017differentiation],
available from NCBI GEO with the accession code GSE106273.


```r
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
base.path <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2834nnn/GSM2834500/suppl"
barcode.fname <- bfcrpath(bfc, file.path(base.path, 
    "GSM2834500%5FG%5F1%5Fbarcodes%2Etsv%2Egz"))
gene.fname <- bfcrpath(bfc, file.path(base.path,
    "GSM2834500%5FG%5F1%5Fgenes%2Etsv%2Egz"))
counts.fname <- bfcrpath(bfc, file.path(base.path,
    "GSM2834500%5FG%5F1%5Fmatrix%2Emtx%2Egz"))
```

# Preparing the data

## Reading in the counts

We create a `SingleCellExperiment` object from the count matrix.
The files have been modified from the _CellRanger_ output, so we have to manually load them in rather than using `read10xCounts()`.


```r
library(scater)
library(Matrix)
gene.info <- read.table(gene.fname, stringsAsFactors=FALSE)
colnames(gene.info) <- c("Ensembl", "Symbol")
sce <- SingleCellExperiment(
    list(counts=as(readMM(counts.fname), "dgCMatrix")), 
    rowData=gene.info, 
    colData=DataFrame(Barcode=readLines(barcode.fname))
)
```

We put some more meaningful information in the row and column names.
Note the use of `uniquifyFeatureNames()` to generate unique row names from gene symbols.


```r
rownames(sce) <- uniquifyFeatureNames(
    rowData(sce)$Ensembl, rowData(sce)$Symbol)
colnames(sce) <- sce$Barcode
sce
```

```
## class: SingleCellExperiment 
## dim: 27998 2915 
## metadata(0):
## assays(1): counts
## rownames(27998): Xkr4 Gm1992 ... Vmn2r122 CAAA01147332.1
## rowData names(2): Ensembl Symbol
## colnames(2915): AAACCTGAGGATGCGT-1 AAACCTGGTAGTAGTA-1 ...
##   TTTGTCATCCTTAATC-1 TTTGTCATCGAACGGA-1
## colData names(1): Barcode
## reducedDimNames(0):
## spikeNames(0):
```

We add some genomic location annotation for downstream use.


```r
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
chrloc <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keytype="GENEID", 
    keys=rowData(sce)$Ensembl, column="CDSCHROM")
rowData(sce)$Chr <- chrloc
```

## Quality control

We compute quality control (QC) metrics using the `calculateQCMetrics()` function from the *[scater](https://bioconductor.org/packages/3.9/scater)* package [@mccarthy2017scater].


```r
is.mito <- rowData(sce)$Chr == "chrM"
summary(is.mito)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical   21767      13    6218
```

```r
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(is.mito)))
```

We remove cells that are outliers for any of these metrics, as previously discussed.
Note that some quality control was already performed by the authors of the original study, so relatively few cells are discarded here.


```r
low.lib <- isOutlier(sce$total_counts, log=TRUE, nmads=3, type="lower")
low.nexprs <- isOutlier(sce$total_features_by_counts, log=TRUE, nmads=3, type="lower")
high.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")
discard <- low.lib | low.nexprs | high.mito
DataFrame(LowLib=sum(low.lib), LowNum=sum(low.nexprs), HighMito=sum(high.mito), 
    Discard=sum(discard), Kept=sum(!discard))
```

```
## DataFrame with 1 row and 5 columns
##      LowLib    LowNum  HighMito   Discard      Kept
##   <integer> <integer> <integer> <integer> <integer>
## 1         0         0       143       143      2772
```

We then subset the `SingleCellExperiment` object to remove these low-quality cells.


```r
sce <- sce[,!discard]
```

## Normalization for cell-specific biases

We apply the pooling method [@lun2016pooling] to compute size factors for scaling normalization of cell-specific biases.


```r
library(scran)
set.seed(1000)
sce <- simpleSumFactors(sce, min.mean=0.1, approximate=TRUE)
summary(sizeFactors(sce))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2734  0.5168  0.7479  1.0000  1.2069 11.0371
```

We then compute log-normalized expression values for downstream use.
This data set does not contain spike-in transcripts so separate normalization with `computeSpikeFactors()` is not required.


```r
sce <- normalize(sce)
assayNames(sce)
```

```
## [1] "counts"    "logcounts"
```

## Modelling and removing noise

As we have no spike-ins, we model technical noise using the `makeTechTrend()` function.


```r
tech.trend <- makeTechTrend(x=sce)
fit <- trendVar(sce, use.spikes=FALSE)
plot(fit$mean, fit$var, pch=16, 
    xlab="Mean log-expression",
    ylab="Variance of log-expression")
curve(tech.trend(x), add=TRUE, col="red")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-6-doublet_files/figure-html/varplot-1.png" alt="Variance of the log-expression values as a function of the mean log-expression in the mammary gland data set. Each point represents a gene, and the red line corresponds to Poisson variance." width="100%" />
<p class="caption">(\#fig:varplot)Variance of the log-expression values as a function of the mean log-expression in the mammary gland data set. Each point represents a gene, and the red line corresponds to Poisson variance.</p>
</div>

We use `denoisePCA()` to choose the number of principal components (PCs) to retain based on the technical noise per gene.
We need to set the seed for reproducibility when `approximate=TRUE`, due to the use of randomized methods from *[irlba](https://bioconductor.org/packages/3.9/irlba)*.


```r
set.seed(12345)
sce <- denoisePCA(sce, technical=tech.trend, approximate=TRUE)
ncol(reducedDim(sce))
```

```
## [1] 15
```

## Clustering into subpopulations

We cluster cells into putative subpopulations using `buildSNNGraph()` [@xu2015identification].
We use a higher `k` to increase connectivity and reduce the granularity of the clustering.


```r
snn.gr <- buildSNNGraph(sce, use.dimred="PCA", k=25)
sce$Cluster <- factor(igraph::cluster_walktrap(snn.gr)$membership)
table(sce$Cluster)
```

```
## 
##   1   2   3   4   5   6   7   8   9  10  11 
## 891 470 775  24  79  87 297  25  52  33  39
```

We visualize the clustering on a _t_-SNE plot [@van2008visualizing].
Figure \@ref(fig:tsneclust) shows that there are a number of well-separated clusters as well as some more inter-related clusters.


```r
set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA")
plotTSNE(sce, colour_by="Cluster")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-6-doublet_files/figure-html/tsneclust-1.png" alt="t-SNE plot of the mammary gland data set. Each point is a cell coloured according to its assigned cluster identity." width="100%" />
<p class="caption">(\#fig:tsneclust)t-SNE plot of the mammary gland data set. Each point is a cell coloured according to its assigned cluster identity.</p>
</div>

# Doublet detection with clusters

The `doubletCluster()` function will identify clusters that have intermediate expression profiles of two other clusters [@bach2017differentiation].
Specifically, it will examine every possible triplet of clusters consisting of a query cluster and its two "parents".
It will then compute a number of statistics:

- The number of genes (`N`) that are differentially expressed in the same direction in the query cluster compared to _both_ of the parent clusters.
Such genes would be unique markers for the query cluster and provide evidence against the null hypothesis, i.e., that the query cluster consists of doublets from the two parents.
Clusters with few unique genes are more likely to be doublets.
- The ratio of the median library size in each parent to the median library size in the query (`lib.size*`).
Doublet libraries are generated from a larger initial pool of RNA compared to libraries for single cells, and thus the former should have larger library sizes.
Library size ratios much greater than unity are inconsistent with a doublet identity for the query.
- The proportion of cells in the query cluster should also be reasonable - typically less than 5% of all cells, depending on how many cells were loaded onto the 10X Genomics device.


```r
dbl.out <- doubletCluster(sce, sce$Cluster)
dbl.out
```

```
## DataFrame with 11 rows and 9 columns
##        source1     source2         N        best              p.value
##    <character> <character> <integer> <character>            <numeric>
## 6            2           1         2         Fau    0.027481774634694
## 5            4           2        31         Ptn 1.51830177144041e-11
## 7            6           4        78       Cotl1 7.84265653798545e-09
## 3            8           1       120        Pigr 4.85238637379419e-25
## 2           10           5       136      Igfbp5 2.50430821433591e-10
## 4            9           8       140        Gde1 9.30194891490462e-12
## 8           11           4       197       Epcam 4.90533186350571e-20
## 9           11           4       270    AF251705   3.296606994399e-24
## 10          11           4       300       Fabp4 2.70725398963721e-32
## 11          10           8       388         Dcn 4.93706079643116e-32
## 1            6           4       575       Cpeb3 4.35646192474353e-20
##            lib.size1         lib.size2                prop
##            <numeric>         <numeric>           <numeric>
## 6  0.542929292929293   1.1996632996633  0.0313852813852814
## 5  0.961305817470201  1.14748265433197  0.0284992784992785
## 7   1.63658906185425 0.744386279101805   0.107142857142857
## 3  0.662732784755305  1.54309224772629    0.27958152958153
## 2   1.25767441860465 0.871472868217054    0.16955266955267
## 4   1.16785416859443  1.13278430646803 0.00865800865800866
## 8  0.882698905407613 0.882780591406633 0.00901875901875902
## 9  0.856192060850963 0.856271293875287  0.0187590187590188
## 10 0.666050295857988 0.666111932938856  0.0119047619047619
## 11  1.50138811771238  1.13288913566537  0.0140692640692641
## 1  0.833567218635981 0.379139769856862   0.321428571428571
##                                    all.pairs
##                              <DataFrameList>
## 6         2:1:2:...,3:2:3:...,7:5:10:...,...
## 5     4:2:31:...,10:2:81:...,8:2:130:...,...
## 7     6:4:78:...,8:6:146:...,9:6:238:...,...
## 3   8:1:120:...,8:6:161:...,11:1:213:...,...
## 2  10:5:136:...,11:5:186:...,8:5:228:...,...
## 4   9:8:140:...,9:5:156:...,10:9:157:...,...
## 8   11:4:197:...,4:3:240:...,6:4:240:...,...
## 9  11:4:270:...,10:4:338:...,8:4:366:...,...
## 10 11:4:300:...,11:8:329:...,4:2:336:...,...
## 11 10:8:388:...,10:9:403:...,9:5:431:...,...
## 1    6:4:575:...,5:3:730:...,8:6:816:...,...
```

Examination of the output of `doubletCluster()` indicates that cluster 6 has the fewest unique genes and library sizes that are comparable to or greater than its parents.
We see that every gene detected in this cluster is also expressed in either of the two proposed parent clusters (Figure \@ref(fig:heatclust)).


```r
markers <- findMarkers(sce, sce$Cluster, direction="up")
dbl.markers <- markers[["6"]]
chosen <- rownames(dbl.markers)[dbl.markers$Top <= 10]
plotHeatmap(sce, columns=order(sce$Cluster), colour_columns_by="Cluster", 
    features=chosen, cluster_cols=FALSE, center=TRUE, symmetric=TRUE, 
    zlim=c(-5, 5), show_colnames=FALSE)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-6-doublet_files/figure-html/heatclust-1.png" alt="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 6 in the mammary gland dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend." width="100%" />
<p class="caption">(\#fig:heatclust)Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 6 in the mammary gland dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend.</p>
</div>



Closer examination of some known markers suggests that the offending cluster consists of doublets of basal cells (_Acta2_) and alveolar cells (_Csn2_) (Figure \@ref(fig:markerexprs)).
Indeed, no cell type is known to strongly express both of these genes at the same time, which supports the hypothesis that this cluster consists solely of doublets.
Of course, it is possible that this cluster represents an entirely novel cell type, though the presence of doublets provides a more sober explanation for its expression profile.


```r
plotExpression(sce, features=c("Acta2", "Csn2"), 
    x="Cluster", colour_by="Cluster")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-6-doublet_files/figure-html/markerexprs-1.png" alt="Distribution of log-normalized expression values for _Acta2_ and _Csn2_ in each cluster. Each point represents a cell." width="960" />
<p class="caption">(\#fig:markerexprs)Distribution of log-normalized expression values for _Acta2_ and _Csn2_ in each cluster. Each point represents a cell.</p>
</div>

The strength of `doubletCluster()` lies in its simplicity and ease of interpretation.
Suspect clusters can be quickly flagged for further investigation, based on the metrics returned by the function.
However, it is obviously dependent on the quality of the clustering.
Clusters that are too coarse will fail to separate doublets from other cells, while clusters that are too fine will complicate interpretation.

**Comments from Aaron:**

- The output of `doubletClusters()` should be treated as a prioritization of "high-risk" clusters that require more careful investigation.
We do not recommend using a fixed threshold on any of the metrics to identify doublet clusters.
Any appropriate threshold for the metrics will depend on the quality of the clustering and the biological context.
- The pair of parents for each query cluster are identified solely on the basis of the lowest `N`.
This means that any `lib.size*` above unity is not definitive evidence against a doublet identity for a query cluster.
It is possible for the "true" parent clusters to be adjacent to the detected parents but with slightly higher `N`.
If this occurs, inspect the `all.pairs` field for statistics on all possible parent pairs for a given query cluster. 
- Clusters with few cells are implicitly more likely to be detected as doublets.
This is because they will have less power to detect DE genes and thus the value of `N` is more likely to be small.
Fortunately for us, this is a desirable effect as doublets should be rare in a properly performed scRNA-seq experiment.

# Doublet detection by simulation

## Background

The other doublet detection strategy involves _in silico_ simulation of doublets from the single-cell expression profiles [@dahlin2018single].
This is performed using the `doubletCells()` function from *[scran](https://bioconductor.org/packages/3.9/scran)*, which will:

1. Simulate thousands of doublets by adding together two randomly chosen single-cell profiles.
2. For each original cell, compute the density of simulated doublets in the surrounding neighbourhood.
3. For each original cell, compute the density of other observed cells in the neighbourhood.
4. Return the ratio between the two densities as a "doublet score" for each cell.

This approach assumes that the simulated doublets are good approximations for real doublets.
The use of random selection accounts for the relative abundances of different subpopulations, which affect the likelihood of their involvement in doublets;
and the calculation of a ratio avoids high scores for non-doublet cells in highly abundant subpopulations.
We see the function in action below:


```r
set.seed(100)
dbl.dens <- doubletCells(sce, approximate=TRUE)
summary(dbl.dens)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00000 0.02052 0.06882 0.15756 0.13988 7.90837
```

The highest doublet scores are concentrated in a single cluster of cells in the centre of Figure \@ref(fig:denstsne).


```r
sce$DoubletScore <- dbl.dens
plotTSNE(sce, colour_by="DoubletScore")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-6-doublet_files/figure-html/denstsne-1.png" alt="t-SNE plot of the mammary gland data set. Each point is a cell coloured according to its doublet density." width="100%" />
<p class="caption">(\#fig:denstsne)t-SNE plot of the mammary gland data set. Each point is a cell coloured according to its doublet density.</p>
</div>

From the clustering information, we see that the affected cells belong to the same cluster that was identified using `doubletCluster()` (Figure \@ref(fig:densclust)). 


```r
plotColData(sce, x="Cluster", y="DoubletScore", colour_by="Cluster")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-6-doublet_files/figure-html/densclust-1.png" alt="Distribution of doublet scores for each cluster in the mammary gland data set. Each point is a cell." width="100%" />
<p class="caption">(\#fig:densclust)Distribution of doublet scores for each cluster in the mammary gland data set. Each point is a cell.</p>
</div>

The advantage of `doubletCells()` is that it does not depend on clusters, reducing the sensitivity of the results to clustering quality.
The downside is that it requires some strong assumptions about how doublets form, such as the combining proportions and the sampling from pure subpopulations.
In particular, `doubletCells()` treats the library size of each cell as an accurate proxy for its total RNA content.
If this is not true, the simulation will not combine expression profiles from different cells in the correct proportions.
This means that the simulated doublets will be systematically shifted away from the real doublets, resulting in doublet scores that are too low.

It should also be stressed that the interpretation of the doublet scores is relative.
It is difficult to generally define a fixed threshold above which libraries are to be considered doublets.
One strategy is to define a threshold at a percentile corresponding to the expected proportion of doublets, based on experimental knowledge of the number of loaded cells.
Another approach is to use these scores in the context of cluster annotation, where clusters with higher-than-average doublet scores are considered suspect.

**Comments from Aaron:**

- To speed up the density calculations, `doubletCells()` will perform a PCA on the log-expression matrix. 
When `approximate=TRUE`, methods from the *[irlba](https://CRAN.R-project.org/package=irlba)* package are used to perform a fast approximate PCA.
This involves randomization so it is necessary to call `set.seed()` to ensure that results are reproducible.
- In some cases, we can improve the clarity of the result by setting `force.match=TRUE` in the `doubletCells()` call.
This will forcibly match each simulated doublet to the nearest neighbouring cells in the original data set.
Any systematic differences between simulated and real doublets will be removed, provided that the former are close enough to the latter to identify the correct nearest neighbours.
This overcomes some of issues related to RNA content but is a rather aggressive strategy that may incorrectly inflate the reported doublet scores.
- The issue of unknown combining proportions can be solved completely if spike-in information is available, e.g., in plate-based protocols.
This will provide an accurate estimate of the total RNA content of each cell.
To this end, size factors from `computeSpikeFactors()` (see [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/xtra-2-spike)) can be supplied to the `doubletCells()` function via the `size.factors.content=` argument.
This will use the spike-in size factors to scale the contribution of each cell to a doublet library.

# Concluding remarks

Doublet detection procedures should only be applied to libraries generated in the same experimental batch.
It is obviously impossible for doublets to form between two cells that were captured separately.
Thus, some understanding of the experimental design is required prior to the use of the above functions.
This avoids unnecessary concerns about the validity of clusters when they cannot possibly consist of doublets.

It is also difficult to interpret doublet predictions in data containing cellular trajectories.
By definition, cells in the middle of a trajectory are always intermediate between other cells and are liable to be incorrectly detected as doublets.
Some protection is provided by the non-linear nature of many real trajectories, which reduces the risk of simulated doublets coinciding with real cells in `doubletCells()`.
One can also put more weight on the relative library sizes in `doubletCluster()` instead of relying on `N`, 
under the assumption that sudden spikes in RNA content are unlikely in a continuous biological process.

# Session information

We save the `SingleCellExperiment` object with its associated data to file for future use.


```r
saveRDS(sce, file="mammary.rds")
```

All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org).
The specific version numbers of the packages used are shown below, along with the version of the R installation.


```r
sessionInfo()
```

```
## R Under development (unstable) (2018-11-02 r75535)
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
##  [1] scran_1.11.1                          
##  [2] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0
##  [3] GenomicFeatures_1.35.1                
##  [4] AnnotationDbi_1.45.0                  
##  [5] Matrix_1.2-15                         
##  [6] scater_1.11.2                         
##  [7] ggplot2_3.1.0                         
##  [8] SingleCellExperiment_1.5.0            
##  [9] SummarizedExperiment_1.13.0           
## [10] DelayedArray_0.9.0                    
## [11] BiocParallel_1.17.1                   
## [12] matrixStats_0.54.0                    
## [13] Biobase_2.43.0                        
## [14] GenomicRanges_1.35.0                  
## [15] GenomeInfoDb_1.19.0                   
## [16] IRanges_2.17.1                        
## [17] S4Vectors_0.21.1                      
## [18] BiocGenerics_0.29.1                   
## [19] bindrcpp_0.2.2                        
## [20] BiocFileCache_1.7.0                   
## [21] dbplyr_1.2.2                          
## [22] knitr_1.20                            
## [23] BiocStyle_2.11.0                      
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] RColorBrewer_1.1-2       progress_1.2.0          
##  [5] httr_1.3.1               rprojroot_1.3-2         
##  [7] dynamicTreeCut_1.63-1    tools_3.6.0             
##  [9] backports_1.1.2          irlba_2.3.3             
## [11] R6_2.3.0                 HDF5Array_1.11.0        
## [13] vipor_0.4.5              DBI_1.0.0               
## [15] lazyeval_0.2.1           colorspace_1.3-2        
## [17] withr_2.1.2              tidyselect_0.2.5        
## [19] gridExtra_2.3            prettyunits_1.0.2       
## [21] bit_1.1-14               curl_3.2                
## [23] compiler_3.6.0           BiocNeighbors_1.1.1     
## [25] labeling_0.3             rtracklayer_1.43.0      
## [27] bookdown_0.7             scales_1.0.0            
## [29] rappdirs_0.3.1           stringr_1.3.1           
## [31] digest_0.6.18            Rsamtools_1.35.0        
## [33] rmarkdown_1.10           XVector_0.23.0          
## [35] pkgconfig_2.0.2          htmltools_0.3.6         
## [37] highr_0.7                limma_3.39.1            
## [39] rlang_0.3.0.1            RSQLite_2.1.1           
## [41] DelayedMatrixStats_1.5.0 bindr_0.1.1             
## [43] dplyr_0.7.8              RCurl_1.95-4.11         
## [45] magrittr_1.5             GenomeInfoDbData_1.2.0  
## [47] Rcpp_1.0.0               ggbeeswarm_0.6.0        
## [49] munsell_0.5.0            Rhdf5lib_1.5.0          
## [51] viridis_0.5.1            edgeR_3.25.0            
## [53] stringi_1.2.4            yaml_2.2.0              
## [55] zlibbioc_1.29.0          Rtsne_0.15              
## [57] rhdf5_2.27.0             plyr_1.8.4              
## [59] grid_3.6.0               blob_1.1.1              
## [61] crayon_1.3.4             lattice_0.20-38         
## [63] cowplot_0.9.3            Biostrings_2.51.1       
## [65] hms_0.4.2                locfit_1.5-9.1          
## [67] pillar_1.3.0             igraph_1.2.2            
## [69] reshape2_1.4.3           biomaRt_2.39.2          
## [71] XML_3.98-1.16            glue_1.3.0              
## [73] evaluate_0.12            BiocManager_1.30.3      
## [75] gtable_0.2.0             purrr_0.2.5             
## [77] assertthat_0.2.0         xfun_0.4                
## [79] viridisLite_0.3.0        pheatmap_1.0.10         
## [81] tibble_1.4.2             GenomicAlignments_1.19.0
## [83] beeswarm_0.2.3           memoise_1.1.0           
## [85] statmod_1.4.30
```

# References
