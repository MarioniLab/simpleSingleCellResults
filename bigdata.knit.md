---
title: Scalable analyses for big scRNA-seq data with Bioconductor 
author:
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
date: "2019-01-04"
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
##   bptimeout: 2592000; bpprogressbar: FALSE; bpexportglobals: TRUE
##   bpRNGseed: 
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
##   bptimeout: 2592000; bpprogressbar: FALSE; bpexportglobals: TRUE
##   bpRNGseed: 
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

# Approximate nearest neighbours searching

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
##    1    0 508   1   0   0   0   0   2   0   0   0   0   0   0
##    2   53   1   0   0   1   0   0   0   0   0   0   0   0   0
##    3    0   0   0 205   0   0   0   0   1   0   0   0   0   0
##    4    2   0   0   0 726   1   0   0   0   0   0   0   0   0
##    5    0   0 545   0   0   0   0   0   0   0   0   0   0   0
##    6    0   0   1   0   0 523   0   0   0   0   0   0   0   0
##    7    0   0   0   0   0   0 127   0   0   0   0   0   0   0
##    8    0  37   0   0   0   0   0 786   0   0   0   0   0   0
##    9    0  45   0   0   0   0   0   1   0   0 108   0   0   0
##    10   0   0   0   0   0   0   0   0  40   0   0   0   0   0
##    11   0   0   0   0   3   0   0   0   0  59   0   0   0   0
##    12   0   0   0   0   0   0   0   0   0   0   0  85   0   0
##    13   0   0   0   0   0   0   0   0   0   0   0   0   0  14
##    14   0   0   0   0   0   0   0   0   0   0   0   0  46   0
```

**Comments from Aaron:**

- The neighbour search algorithms are interoperable with *[BiocParallel](https://bioconductor.org/packages/3.9/BiocParallel)*, so it is straightforward to parallelize the search for greater speed.

# Concluding remarks

All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org).
The specific version numbers of the packages used are shown below, along with the version of the R installation.


```r
sessionInfo()
```

```
## R Under development (unstable) (2018-12-07 r75787)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: OS X El Capitan 10.11.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] BiocNeighbors_1.1.7         scran_1.11.12              
##  [3] scater_1.11.5               ggplot2_3.1.0              
##  [5] TENxBrainData_1.3.0         HDF5Array_1.11.10          
##  [7] rhdf5_2.27.4                SingleCellExperiment_1.5.1 
##  [9] SummarizedExperiment_1.13.0 DelayedArray_0.9.5         
## [11] BiocParallel_1.17.3         matrixStats_0.54.0         
## [13] Biobase_2.43.0              GenomicRanges_1.35.1       
## [15] GenomeInfoDb_1.19.1         IRanges_2.17.3             
## [17] S4Vectors_0.21.8            BiocGenerics_0.29.1        
## [19] knitr_1.21                  BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6                  bit64_0.9-7                  
##  [3] httr_1.4.0                    dynamicTreeCut_1.63-1        
##  [5] tools_3.6.0                   R6_2.3.0                     
##  [7] vipor_0.4.5                   DBI_1.0.0                    
##  [9] lazyeval_0.2.1                colorspace_1.3-2             
## [11] withr_2.1.2                   tidyselect_0.2.5             
## [13] gridExtra_2.3                 processx_3.2.1               
## [15] bit_1.1-14                    curl_3.2                     
## [17] compiler_3.6.0                bookdown_0.9                 
## [19] scales_1.0.0                  callr_3.1.1                  
## [21] stringr_1.3.1                 digest_0.6.18                
## [23] rmarkdown_1.11                XVector_0.23.0               
## [25] pkgconfig_2.0.2               htmltools_0.3.6              
## [27] limma_3.39.3                  rlang_0.3.0.1                
## [29] RSQLite_2.1.1                 shiny_1.2.0                  
## [31] DelayedMatrixStats_1.5.0      bindr_0.1.1                  
## [33] dplyr_0.7.8                   RCurl_1.95-4.11              
## [35] magrittr_1.5                  simpleSingleCell_1.7.10      
## [37] GenomeInfoDbData_1.2.0        Matrix_1.2-15                
## [39] Rcpp_1.0.0                    ggbeeswarm_0.6.0             
## [41] munsell_0.5.0                 Rhdf5lib_1.5.1               
## [43] viridis_0.5.1                 edgeR_3.25.2                 
## [45] stringi_1.2.4                 yaml_2.2.0                   
## [47] zlibbioc_1.29.0               plyr_1.8.4                   
## [49] AnnotationHub_2.15.3          grid_3.6.0                   
## [51] blob_1.1.1                    promises_1.0.1               
## [53] ExperimentHub_1.9.0           crayon_1.3.4                 
## [55] lattice_0.20-38               locfit_1.5-9.1               
## [57] ps_1.3.0                      pillar_1.3.1                 
## [59] igraph_1.2.2                  codetools_0.2-16             
## [61] glue_1.3.0                    evaluate_0.12                
## [63] BiocManager_1.30.4            httpuv_1.4.5.1               
## [65] gtable_0.2.0                  purrr_0.2.5                  
## [67] assertthat_0.2.0              xfun_0.4                     
## [69] mime_0.6                      xtable_1.8-3                 
## [71] later_0.7.5                   viridisLite_0.3.0            
## [73] tibble_1.4.2                  AnnotationDbi_1.45.0         
## [75] beeswarm_0.2.3                memoise_1.1.0                
## [77] bindrcpp_0.2.2                statmod_1.4.30               
## [79] interactiveDisplayBase_1.21.0
```

# References 
