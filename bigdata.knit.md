---
title: Scalable analyses for big scRNA-seq data with Bioconductor 
author:
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
date: "2019-04-27"
vignette: >
  %\VignetteIndexEntry{11. Scalability for big data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---





# Overview

Advances in single-cell RNA sequencing (scRNA-seq) technologies have increased the number of cells that can be assayed in routine experiments.
For effective data analysis, the computational methods need to scale with the increasing size of scRNA-seq data sets.
Scalability requires greater use of parallelization, out-of-memory representations and fast approximate algorithms to process data efficiently.
Fortunately, this is easy to achieve within the Bioconductor ecosystem.
This workflow will discuss how to tune the previous analysis pipelines for greater speed to handle large scRNA-seq data sets.

# Out of memory representations

As we have previously discussed, the count matrix is the central structure around which our analyses are based.
In the previous workflows, this has been held fully in memory as a dense `matrix` or as a sparse `dgCMatrix`.
Howevever, in-memory representations may not be feasible for very large data sets, especially on machines with limited memory.
For example, the 1.3 million brain cell data set from 10X Genomics [@zheng2017massively] would require over 100 GB of RAM to hold as a `matrix` and around 30 GB as a `dgCMatrix`.
This makes it challenging to investigate the data on anything less than a high-performance computing system.

The obvious solution is to use a file-backed matrix representation where the data are held on disk and subsets are retrieved into memory as requested.
While a number of implementations of file-backed matrices are available (e.g., *[bigmemory](https://CRAN.R-project.org/package=bigmemory)*, *[matter](https://bioconductor.org/packages/3.9/matter)*),
we will be using the implementation from the *[HDF5Array](https://bioconductor.org/packages/3.9/HDF5Array)* package.
This uses the popular HDF5 format as the underlying data store, which provides a measure of standardization and portability across systems. 
We demonstrate with a subset of 20,000 cells from the 1.3 million brain cell data set, as provided by the *[TENxBrainData](https://bioconductor.org/packages/3.9/TENxBrainData)* package^[We could instead obtain the full-sized data set by using `TENxBrainData()`, but we will use the smaller data set here for demonstration purposes.].


```r
library(TENxBrainData)
sce <- TENxBrainData20k() # downloads once and caches it for future use.
sce
```

```
## class: SingleCellExperiment 
## dim: 27998 20000 
## metadata(0):
## assays(1): counts
## rownames: NULL
## rowData names(2): Ensembl Symbol
## colnames: NULL
## colData names(4): Barcode Sequence Library Mouse
## reducedDimNames(0):
## spikeNames(0):
```

Examination of the `SingleCellExperiment` object indicates that the count matrix is a `HDF5Matrix`.
From a comparison of the memory usage, it is clear that this matrix object is simply a stub that points to the much larger HDF5 file that actually contains the data.
This avoids the need for large RAM availability during analyses.


```r
counts(sce)
```

```
## <27998 x 20000> HDF5Matrix object of type "integer":
##              [,1]     [,2]     [,3]     [,4] ... [,19997] [,19998] [,19999]
##     [1,]        0        0        0        0   .        0        0        0
##     [2,]        0        0        0        0   .        0        0        0
##     [3,]        0        0        0        0   .        0        0        0
##     [4,]        0        0        0        0   .        0        0        0
##     [5,]        0        0        0        0   .        0        0        0
##      ...        .        .        .        .   .        .        .        .
## [27994,]        0        0        0        0   .        0        0        0
## [27995,]        0        0        0        1   .        0        2        0
## [27996,]        0        0        0        0   .        0        1        0
## [27997,]        0        0        0        0   .        0        0        0
## [27998,]        0        0        0        0   .        0        0        0
##          [,20000]
##     [1,]        0
##     [2,]        0
##     [3,]        0
##     [4,]        0
##     [5,]        0
##      ...        .
## [27994,]        0
## [27995,]        0
## [27996,]        0
## [27997,]        0
## [27998,]        0
```

```r
object.size(counts(sce))
```

```
## 2144 bytes
```

```r
file.info(path(counts(sce)))$size
```

```
## [1] 76264332
```

Manipulation of the count matrix will generally result in the creation of a `DelayedArray` (from the *[DelayedArray](https://bioconductor.org/packages/3.9/DelayedArray)* package).
This stores delayed operations in the matrix object, to be executed when the modified matrix values are realized for use in calculations. 
The use of delayed operations avoids the need to write the modified values to a new file at every operation, which would unnecessarily require time-consuming disk I/O.


```r
tmp <- counts(sce)
tmp <- log2(tmp + 1)
tmp
```

```
## <27998 x 20000> DelayedMatrix object of type "double":
##              [,1]     [,2]     [,3] ... [,19999] [,20000]
##     [1,]        0        0        0   .        0        0
##     [2,]        0        0        0   .        0        0
##     [3,]        0        0        0   .        0        0
##     [4,]        0        0        0   .        0        0
##     [5,]        0        0        0   .        0        0
##      ...        .        .        .   .        .        .
## [27994,]        0        0        0   .        0        0
## [27995,]        0        0        0   .        0        0
## [27996,]        0        0        0   .        0        0
## [27997,]        0        0        0   .        0        0
## [27998,]        0        0        0   .        0        0
```

Many functions described in the previous workflows are capable of accepting `HDF5Matrix` objects^[If you find one that is not, please contact the maintainers.].
This is powered by the availability of common methods for all matrix representations (e.g., subsetting, combining, methods from *[DelayedMatrixStats](https://bioconductor.org/packages/3.9/DelayedMatrixStats)*)
as well as representation-agnostic C++ code using *[beachmat](https://bioconductor.org/packages/3.9/beachmat)* [@lun2018beachmat].
For example, we compute quality control (QC) metrics below with the same `calculateQCMetrics()` function that we used in the other workflows.


```r
library(scater)
sce <- calculateQCMetrics(sce, compact=TRUE) # compacting for clean output.
sce$scater_qc
```

```
## DataFrame with 20000 rows and 2 columns
##       is_cell_control                            all
##             <logical>                    <DataFrame>
## 1               FALSE 1546:3.18949031369937:3060:...
## 2               FALSE  1694:3.2291697025391:3500:...
## 3               FALSE 1613:3.20790353038605:3092:...
## 4               FALSE 2050:3.31196566036837:4420:...
## 5               FALSE 1813:3.25863728272408:3771:...
## ...               ...                            ...
## 19996           FALSE 2050:3.31196566036837:4431:...
## 19997           FALSE 2704:3.43216726944259:6988:...
## 19998           FALSE 2988:3.47552591503928:8749:...
## 19999           FALSE 1711:3.23350376034113:3842:...
## 20000           FALSE  945:2.97589113640179:1775:...
```

Needless to say, data access from file-backed representations is slower than that from in-memory representations (assuming the latter is not moved into swap space).
The time spent retrieving data from disk is an unavoidable cost of memory efficiency.

**Comments from Aaron:**

- By default, file locking is necessary for reading from HDF5 files via the *[Rhdf5lib](https://bioconductor.org/packages/3.9/Rhdf5lib)* library, but this may be disabled on some file systems.
Users can set the `HDF5_USE_FILE_LOCKING` environment variable to `FALSE` to avoid this requirement.

# Parallelization

In many Bioconductor packages, different parallelization mechanisms are easily tested through the *[BiocParallel](https://bioconductor.org/packages/3.9/BiocParallel)* framework.
We construct a `BiocParallelParam` object that specifies the type of parallelization that we wish to use.
For example, we might use forking^[Not available on Windows.] across 2 cores:


```r
bpp <- MulticoreParam(2)
bpp
```

```
## class: MulticoreParam
##   bpisup: FALSE; bpnworkers: 2; bptasks: 0; bpjobname: BPJOB
##   bplog: FALSE; bpthreshold: INFO; bpstopOnError: TRUE
##   bpRNGseed: ; bptimeout: 2592000; bpprogressbar: FALSE
##   bpexportglobals: TRUE
##   bplogdir: NA
##   bpresultdir: NA
##   cluster type: FORK
```

Another approach would be to distribute jobs across a network of computers:


```r
bpp <- SnowParam(5)
bpp
```

```
## class: SnowParam
##   bpisup: FALSE; bpnworkers: 5; bptasks: 0; bpjobname: BPJOB
##   bplog: FALSE; bpthreshold: INFO; bpstopOnError: TRUE
##   bpRNGseed: ; bptimeout: 2592000; bpprogressbar: FALSE
##   bpexportglobals: TRUE
##   bplogdir: NA
##   bpresultdir: NA
##   cluster type: SOCK
```

High-performance computing systems typically use job schedulers across a cluster of compute nodes.
We can distribute jobs via the scheduler using the `BatchtoolsParam` class.
The example below assumes a SLURM cluster, though the settings can be easily^[In general. Some fiddling may be required, depending on the idiosyncrasies of the cluster set-up.] configured for a particular system (see [here](https://bioconductor.org/packages/3.9/BiocParallel/vignettes/BiocParallel_BatchtoolsParam.pdf) for details).


```r
bpp <- BatchtoolsParam(10, cluster="slurm",
	resources=list(walltime=20000, memory=8000, ncpus=1))
```

Once we have defined the parallelization mechanism, we can pass the `BiocParallelParam` object to the function that we wish to run.
This will instruct the function to run operations in parallel where it is allowed to (as defined by the developer).
Different functions may parallelize operations across cells, or genes, or batches of data, depending on what is most appropriate.
In the example below, we parallelize the QC calculations (across cells) using two cores:


```r
alt <- calculateQCMetrics(sce, BPPARAM=MulticoreParam(2), compact=TRUE)
```

This yields the same result as the single-core calculation, but faster.




```r
all.equal(alt, sce) 
```

```
## [1] TRUE
```

**Comments from Aaron:**

- Efficiently combining parallelization with file-backed matrix representations is likely to require systems that support parallel I/O.

# Fast approximate algorithms

## Nearest neighbours searching

Identification of neighbouring cells in PC or expression space is a common procedure that is used in many functions, e.g., `buildSNNGraph()`, `doubletCells()`.
The default is to favour accuracy over speed by using an exact nearest neighbour search, implemented with the k-means for k-nearest neighbours algorithm [@wang2012fast]. 
However, for large data sets, it may be preferable to use a faster approximate approach.
The *[BiocNeighbors](https://bioconductor.org/packages/3.9/BiocNeighbors)* framework makes it easy to switch between search options.

To demonstrate, we will use the PBMC data from the [previous workflow](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/tenx.html):


```r
sce.pbmc <- readRDS("pbmc_data.rds")
```

We had [previously](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/tenx.html#clustering-with-graph-based-methods) generated a shared nearest neighbor graph with an exact neighbour search.
We repeat this below using an approximate search, implemented using the [Annoy](https://github.com/spotify/Annoy) algorithm.
This involves constructing a `BiocNeighborParam` object to specify the search algorithm, and passing it to the `buildSNNGraph()` function.


```r
library(scran)
library(BiocNeighbors)
snn.gr <- buildSNNGraph(sce.pbmc, BNPARAM=AnnoyParam(), use.dimred="PCA")
```

The results from the exact and approximate searches are consistent with most clusters from the former re-appearing in the latter.
This suggests that the inaccuracy from the approximation can be largely ignored.
However, if the approximation was unacceptable, it would be simple to switch back to an exact algorithm by altering `BNPARAM`.


```r
clusters <- igraph::cluster_walktrap(snn.gr)
table(Exact=sce.pbmc$Cluster, Approx=clusters$membership)
```

```
##      Approx
## Exact   1   2   3   4   5   6   7   8   9  10  11  12  13  14
##    1    1   0   0   0  13 732   0   0   1   0  70   1   0   0
##    2  495   0   0  16   0   0   0  25   0   0   0   0   0   0
##    3    0 124   0   0   0   0   0   0   0   0   0   0   0   0
##    4    1   0   0 527   0   0  25   0   0   0   0   0   0   0
##    5    0   0 515   0   0   0   0   0   0   0   0   0   0   0
##    6    0   0   0   0   0   0 173   0   0   0   0   0   0   0
##    7    3   0   0   0   0   0   0 833   0   0   0   0   0   0
##    8    0   0   0   0  37   3   0   1   0   0   0   0   0   0
##    9    0   0   0   0   0   0   0   0  45   0   0   0   0   0
##    10   0   0   0   0   0   0   0   6   0 143   0   0   0   0
##    11   0   0   0   0   0   0   0   0   0   0   0  81   0   0
##    12   0   0   0   0   0   0   0   0   0   0   0   0  17   0
##    13   0   0   0   0   0   0   0   0   0   0   0   0   0  37
```

**Comments from Aaron:**

- The neighbour search algorithms are interoperable with *[BiocParallel](https://bioconductor.org/packages/3.9/BiocParallel)*, so it is straightforward to parallelize the search for greater speed.

## Principal components analysis

We have already introduced the `IrlbaParam()` function, which creates an `IrlbaParam` parameter object from the *[BiocSingular](https://bioconductor.org/packages/3.9/BiocSingular)* package.
If passed to a function via the `BSPARAM=` argument, this instructs the function to use fast approximate methods from *[irlba](https://CRAN.R-project.org/package=irlba)* to perform SVD or PCA.
However, this is not the only SVD strategy that is provided in the *[BiocSingular](https://bioconductor.org/packages/3.9/BiocSingular)* framework.
For example, one could perform randomized SVD^[Via the *[rsvd](https://CRAN.R-project.org/package=rsvd)* package.] to perform an approximate SVD.


```r
library(BiocSingular)

# As the name suggests, it is random, so we need to set the seed.
set.seed(999)    
r.out <- BiocSingular::runPCA(t(logcounts(sce.pbmc)), rank=20, 
    BSPARAM=RandomParam(deferred=TRUE))
str(r.out)
```

```
## List of 3
##  $ sdev    : num [1:20] 10.54 7.03 5.54 4.3 3.54 ...
##  $ rotation: num [1:33694, 1:20] 1.28e-16 -2.31e-17 -1.47e-17 -9.36e-05 -8.34e-05 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
##  $ x       : num [1:3925, 1:20] -17.25 -16.11 9.57 9.4 -8.15 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
```

```r
# For comparison:    
i.out <- BiocSingular::runPCA(t(logcounts(sce.pbmc)), rank=20, 
    BSPARAM=IrlbaParam(fold=Inf))
str(i.out)
```

```
## List of 3
##  $ sdev    : num [1:20] 10.54 7.03 5.54 4.3 3.54 ...
##  $ rotation: num [1:33694, 1:20] -1.17e-19 3.36e-19 -5.18e-20 -9.36e-05 -8.35e-05 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
##  $ x       : num [1:3925, 1:20] -17.25 -16.11 9.57 9.4 -8.15 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : NULL
##   .. ..$ : chr [1:20] "PC1" "PC2" "PC3" "PC4" ...
```

Users can often obtain considerable speed-ups through a careful choice of `BSPARAM=` that considers the structure of the underlying matrix data.
For example, setting `deferred=TRUE` will defer the centering during matrix multiplications to avoid loss of sparsity.
Another consideration would be whether the data are stored on file, in which case one may prefer an algorithm that performs fewer reads at the expense of more computations.

# Concluding remarks

All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org).
The specific version numbers of the packages used are shown below, along with the version of the R installation.


```r
sessionInfo()
```

```
## R Under development (unstable) (2019-04-11 r76379)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.2 LTS
## 
## Matrix products: default
## BLAS:   /home/luna/Software/R/trunk/lib/libRblas.so
## LAPACK: /home/luna/Software/R/trunk/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] BiocSingular_0.99.18        BiocNeighbors_1.1.13       
##  [3] scran_1.11.27               scater_1.11.16             
##  [5] ggplot2_3.1.1               TENxBrainData_1.3.0        
##  [7] HDF5Array_1.11.12           rhdf5_2.27.19              
##  [9] SingleCellExperiment_1.5.2  SummarizedExperiment_1.13.0
## [11] DelayedArray_0.9.9          BiocParallel_1.17.19       
## [13] matrixStats_0.54.0          Biobase_2.43.1             
## [15] GenomicRanges_1.35.1        GenomeInfoDb_1.19.3        
## [17] IRanges_2.17.5              S4Vectors_0.21.24          
## [19] BiocGenerics_0.29.2         knitr_1.22                 
## [21] BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6                  bit64_0.9-7                  
##  [3] httr_1.4.0                    dynamicTreeCut_1.63-1        
##  [5] tools_3.7.0                   R6_2.4.0                     
##  [7] irlba_2.3.3                   vipor_0.4.5                  
##  [9] DBI_1.0.0                     lazyeval_0.2.2               
## [11] colorspace_1.4-1              withr_2.1.2                  
## [13] tidyselect_0.2.5              gridExtra_2.3                
## [15] processx_3.3.0                bit_1.1-14                   
## [17] curl_3.3                      compiler_3.7.0               
## [19] bookdown_0.9                  scales_1.0.0                 
## [21] callr_3.2.0                   rappdirs_0.3.1               
## [23] stringr_1.4.0                 digest_0.6.18                
## [25] rmarkdown_1.12                XVector_0.23.2               
## [27] pkgconfig_2.0.2               htmltools_0.3.6              
## [29] limma_3.39.18                 dbplyr_1.4.0                 
## [31] rlang_0.3.4                   RSQLite_2.1.1                
## [33] shiny_1.3.2                   DelayedMatrixStats_1.5.2     
## [35] dplyr_0.8.0.1                 RCurl_1.95-4.12              
## [37] magrittr_1.5                  simpleSingleCell_1.7.21      
## [39] GenomeInfoDbData_1.2.1        Matrix_1.2-17                
## [41] Rcpp_1.0.1                    ggbeeswarm_0.6.0             
## [43] munsell_0.5.0                 Rhdf5lib_1.5.4               
## [45] viridis_0.5.1                 edgeR_3.25.7                 
## [47] stringi_1.4.3                 yaml_2.2.0                   
## [49] zlibbioc_1.29.0               plyr_1.8.4                   
## [51] BiocFileCache_1.7.10          AnnotationHub_2.15.15        
## [53] grid_3.7.0                    blob_1.1.1                   
## [55] dqrng_0.2.0                   promises_1.0.1               
## [57] ExperimentHub_1.9.3           crayon_1.3.4                 
## [59] lattice_0.20-38               beachmat_1.99.8              
## [61] locfit_1.5-9.1                ps_1.3.0                     
## [63] pillar_1.3.1                  igraph_1.2.4.1               
## [65] codetools_0.2-16              glue_1.3.1                   
## [67] evaluate_0.13                 BiocManager_1.30.4           
## [69] httpuv_1.5.1                  gtable_0.3.0                 
## [71] purrr_0.3.2                   assertthat_0.2.1             
## [73] xfun_0.6                      rsvd_1.0.0                   
## [75] mime_0.6                      xtable_1.8-4                 
## [77] later_0.8.0                   viridisLite_0.3.0            
## [79] tibble_2.1.1                  AnnotationDbi_1.45.1         
## [81] beeswarm_0.2.3                memoise_1.1.0                
## [83] statmod_1.4.30                interactiveDisplayBase_1.21.0
```

# References 
