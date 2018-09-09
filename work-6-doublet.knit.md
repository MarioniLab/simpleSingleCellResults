---
title: Detecting doublets in single-cell RNA-seq data
author:
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
date: "2018-09-09"
vignette: >
  %\VignetteIndexEntry{8. Detecting doublets in scRNA-seq data}
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
Both are implemented in the *[scran](https://bioconductor.org/packages/3.8/scran)* package from the open-source Bioconductor project [@huber2015orchestrating].
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

We compute quality control (QC) metrics using the `calculateQCMetrics()` function from the *[scater](https://bioconductor.org/packages/3.8/scater)* package [@mccarthy2017scater].


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

We apply the deconvolution method with pre-clustering [@lun2016pooling] to compute size factors for scaling normalization of cell-specific biases.


```r
library(scran)
set.seed(1000)
clusters <- quickCluster(sce, method="igraph", min.mean=0.1)
table(clusters)
```

```
## clusters
##    1    2    3 
##  919 1045  808
```

```r
sce <- computeSumFactors(sce, clusters=clusters, min.mean=0.1)
summary(sizeFactors(sce))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2688  0.5310  0.7672  1.0000  1.1987 10.6856
```

We then compute log-normalized expression values for downstream use.
This data set does not contain spike-in transcripts so separate normalziation with `computeSpikeFactors()` is not required.


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
We need to set the seed for reproducibility when `approximate=TRUE`, due to the use of randomized methods from *[irlba](https://bioconductor.org/packages/3.8/irlba)*.


```r
set.seed(12345)
sce <- denoisePCA(sce, technical=tech.trend, approximate=TRUE)
ncol(reducedDim(sce))
```

```
## [1] 14
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
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14 
## 475 411 398  54 378 271  84  74 289 219  33  22  25  39
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

The `doubletCluster()` function will identify clusters that have intermediate expression profiles of two other clusters [@bach2018differentiation].
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
## DataFrame with 14 rows and 9 columns
##         source1     source2         N        best              p.value
##     <character> <character> <integer> <character>            <numeric>
## 7             8           5         0       Fabp3   0.0516239051106075
## 10           13           2        13        Xist 1.73868501814697e-29
## 8            12           1        28         Ptn 1.10432261800312e-14
## 5            12           7        49       Cotl1 7.90164177238926e-08
## 6             7           2        56      Sec61b  1.3329887878903e-07
## ...         ...         ...       ...         ...                  ...
## 2            12           6       182       Cdc20 1.47756925944611e-19
## 13           14          12       198        Gpx3  1.1282711558589e-19
## 4            14          12       270        C1qb 9.41841774529012e-49
## 11           14          12       299       Fabp4  2.7072539896372e-32
## 14           13          11       388         Dcn 4.93706079643113e-32
##             lib.size1         lib.size2                prop
##             <numeric>         <numeric>           <numeric>
## 7   0.456830005034402 0.574593052525592  0.0303030303030303
## 10  0.677851605758582  1.51240310077519   0.079004329004329
## 8     0.9924694645973  1.18853889246028  0.0266955266955267
## 5           0.7890625  1.74036214953271   0.136363636363636
## 6   0.509882775733721 0.584281680499701  0.0977633477633478
## ...               ...               ...                 ...
## 2   0.395657904371385  1.71150325840228   0.148268398268398
## 13  0.882698905407613 0.882780591406633 0.00901875901875902
## 4   0.856192060850963 0.856271293875287  0.0194805194805195
## 11  0.666050295857988 0.666111932938856  0.0119047619047619
## 14   1.13288913566537  1.50138811771238  0.0140692640692641
##                                         all.pairs
##                                   <DataFrameList>
## 7               8:5:0:...,6:1:1:...,9:1:3:...,...
## 10        13:2:13:...,12:2:26:...,12:3:37:...,...
## 8        12:1:28:...,11:1:97:...,13:1:143:...,...
## 5         12:7:49:...,12:9:55:...,13:9:96:...,...
## 6        7:2:56:...,10:9:145:...,10:7:146:...,...
## ...                                           ...
## 2      12:6:182:...,10:9:248:...,10:8:267:...,...
## 13   14:12:198:...,14:10:200:...,12:3:202:...,...
## 4   14:12:270:...,12:11:331:...,13:12:371:...,...
## 11   14:12:299:...,14:13:329:...,12:1:338:...,...
## 14     13:11:388:...,11:4:406:...,8:4:434:...,...
```

Examination of the output of `doubletCluster()` indicates that cluster 7 has the fewest unique genes and library sizes that are comparable to or greater than its parents.
We see that every gene detected in this cluster is also expressed in either of the two proposed parent clusters (Figure \@ref(fig:heatclust)).


```r
markers <- findMarkers(sce, sce$Cluster, direction="up")
dbl.markers <- markers[["7"]]
chosen <- rownames(dbl.markers)[dbl.markers$Top <= 10]
plotHeatmap(sce, columns=order(sce$Cluster), colour_columns_by="Cluster", 
    features=chosen, cluster_cols=FALSE, center=TRUE, symmetric=TRUE, 
    zlim=c(-5, 5), show_colnames=FALSE)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-6-doublet_files/figure-html/heatclust-1.png" alt="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 7 in the 416B dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend." width="100%" />
<p class="caption">(\#fig:heatclust)Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 7 in the 416B dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend.</p>
</div>



Closer examination of some known markers suggests that the offending cluster consists of doublets of basal cells (_Acta2_) and alveolar cells (_Csn2_) (Figure \@ref(fig:markerexprs))).
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
This is performed using the `doubletCells()` function from *[scran](https://bioconductor.org/packages/3.8/scran)*, which will:

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
## 0.00000 0.01771 0.07163 0.15626 0.13794 7.55793
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

- In some cases, we can improve the clarity of the result by setting `force.match=TRUE` in the `doubletCells()` call.
This will forcibly match each simulated doublet to the nearest neighbouring cells in the original data set.
Any systematic differences between simulated and real doublets will be removed, provided that the former are close enough to the latter to identify the correct nearest neighbours.
This overcomes some of issues related to RNA content but is a rather aggressive strategy that may incorrectly inflate the reported doublet scores.
- The issue of unknown combining proportions can be solved completely if spike-in information is available, e.g., in plate-based protocols.
This will provide an accurate estimate of the total RNA content of each cell.
To this end, size factors from `computeSpikeFactors()` (see [here](https://bioconductor.org/packages/3.8/simpleSingleCell/vignettes/xtra-2-spike)) can be supplied to the `doubletCells()` function via the `size.factors.content=` argument.
This will use the spike-in size factors to scale the contribution of each cell to a doublet library.

# Session information


```r
saveRDS(sce, file="mammary.rds")
```


```r
sessionInfo()
```

```
## R version 3.5.0 Patched (2018-04-30 r74679)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.5 LTS
## 
## Matrix products: default
## BLAS: /home/cri.camres.org/lun01/Software/R/R-3-5-branch/lib/libRblas.so
## LAPACK: /home/cri.camres.org/lun01/Software/R/R-3-5-branch/lib/libRlapack.so
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
##  [1] scran_1.9.20                          
##  [2] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0
##  [3] GenomicFeatures_1.33.2                
##  [4] AnnotationDbi_1.43.1                  
##  [5] Matrix_1.2-14                         
##  [6] scater_1.9.20                         
##  [7] ggplot2_3.0.0                         
##  [8] SingleCellExperiment_1.3.10           
##  [9] SummarizedExperiment_1.11.6           
## [10] DelayedArray_0.7.37                   
## [11] BiocParallel_1.15.11                  
## [12] matrixStats_0.54.0                    
## [13] Biobase_2.41.2                        
## [14] GenomicRanges_1.33.13                 
## [15] GenomeInfoDb_1.17.1                   
## [16] IRanges_2.15.17                       
## [17] S4Vectors_0.19.19                     
## [18] BiocGenerics_0.27.1                   
## [19] bindrcpp_0.2.2                        
## [20] BiocFileCache_1.5.5                   
## [21] dbplyr_1.2.2                          
## [22] knitr_1.20                            
## [23] BiocStyle_2.9.6                       
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] RColorBrewer_1.1-2       progress_1.2.0          
##  [5] httr_1.3.1               rprojroot_1.3-2         
##  [7] dynamicTreeCut_1.63-1    tools_3.5.0             
##  [9] backports_1.1.2          irlba_2.3.2             
## [11] R6_2.2.2                 HDF5Array_1.9.15        
## [13] vipor_0.4.5              DBI_1.0.0               
## [15] lazyeval_0.2.1           colorspace_1.3-2        
## [17] withr_2.1.2              tidyselect_0.2.4        
## [19] gridExtra_2.3            prettyunits_1.0.2       
## [21] bit_1.1-14               compiler_3.5.0          
## [23] labeling_0.3             rtracklayer_1.41.5      
## [25] bookdown_0.7             scales_1.0.0            
## [27] rappdirs_0.3.1           stringr_1.3.1           
## [29] digest_0.6.16            Rsamtools_1.33.5        
## [31] rmarkdown_1.10           XVector_0.21.3          
## [33] pkgconfig_2.0.2          htmltools_0.3.6         
## [35] highr_0.7                limma_3.37.4            
## [37] rlang_0.2.2              RSQLite_2.1.1           
## [39] DelayedMatrixStats_1.3.8 bindr_0.1.1             
## [41] dplyr_0.7.6              RCurl_1.95-4.11         
## [43] magrittr_1.5             GenomeInfoDbData_1.1.0  
## [45] Rcpp_0.12.18             ggbeeswarm_0.6.0        
## [47] munsell_0.5.0            Rhdf5lib_1.3.3          
## [49] viridis_0.5.1            edgeR_3.23.3            
## [51] stringi_1.2.4            yaml_2.2.0              
## [53] zlibbioc_1.27.0          Rtsne_0.13              
## [55] rhdf5_2.25.9             plyr_1.8.4              
## [57] grid_3.5.0               blob_1.1.1              
## [59] crayon_1.3.4             lattice_0.20-35         
## [61] cowplot_0.9.3            Biostrings_2.49.1       
## [63] hms_0.4.2                locfit_1.5-9.1          
## [65] pillar_1.3.0             igraph_1.2.2            
## [67] kmknn_0.99.16            reshape2_1.4.3          
## [69] biomaRt_2.37.6           XML_3.98-1.16           
## [71] glue_1.3.0               evaluate_0.11           
## [73] BiocManager_1.30.2       gtable_0.2.0            
## [75] purrr_0.2.5              assertthat_0.2.0        
## [77] xfun_0.3                 viridisLite_0.3.0       
## [79] pheatmap_1.0.10          tibble_1.4.2            
## [81] GenomicAlignments_1.17.3 beeswarm_0.2.3          
## [83] memoise_1.1.0            statmod_1.4.30
```

# References