---
title: Normalizing single-cell RNA-seq data using spike-in information
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
date: "2018-11-12"
vignette: >
  %\VignetteIndexEntry{07. Spike-in normalization}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---



# Motivation 

Scaling normalization strategies for scRNA-seq data can be broadly divided into two classes.
The first class assumes that there exists a subset of genes that are not DE between samples, as described in the previous workflows.
The second class uses the fact that the same amount of spike-in RNA was added to each cell [@lun2017assessing].
Differences in the coverage of the spike-in transcripts can only be due to cell-specific biases, e.g., in capture efficiency or sequencing depth.
Scaling normalization is then applied to equalize spike-in coverage across cells.

The choice between these two normalization strategies depends on the biology of the cells and the features of interest.
If the majority of genes are expected to be DE and there is no reliable house-keeping set, spike-in normalization may be the only option for removing cell-specific biases.
Spike-in normalization should also be used if differences in the total RNA content of individual cells are of interest.
In any particular cell, an increase in the amount of endogenous RNA will not increase spike-in coverage (with or without library quantification).
Thus, the former will not be represented as part of the bias in the latter, which means that the effects of total RNA content on expression will not be removed upon scaling.
With non-DE normalization, an increase in RNA content will systematically increase the expression of all genes in the non-DE subset, such that it will be treated as bias and removed.

# Setting up the data

## Obtaining the dataset

We demonstrate the use of spike-in normalization on a dataset involving different cell types -- namely, mouse embryonic stem cells (mESCs) and mouse embryonic fibroblasts (MEFs) [@islam2011characterization].
The count table was obtained from the NCBI Gene Expression Omnibus (GEO) as a supplementary file using the accession number [GSE29087](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29087).


```r
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask=FALSE)
islam.fname <- bfcrpath(bfc, file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
    "GSE29nnn/GSE29087/suppl/GSE29087_L139_expression_tab.txt.gz"))
```

We load the counts into R, using `colClasses` to speed up `read.table` by pre-defining the type of each column.
We also specify the rows corresponding to spike-in transcripts.


```r
library(SingleCellExperiment)
counts <- read.table(islam.fname,
    colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), 
    rep("integer", 96)), skip=6, sep='\t', row.names=1)

is.spike <- grep("SPIKE", rownames(counts)) 
sce.islam <- SingleCellExperiment(list(counts=as.matrix(counts)))
isSpike(sce.islam, "spike") <- is.spike
dim(sce.islam)
```

```
## [1] 22936    96
```

## Applying quality control

We perform some quality control to remove low-quality cells using the `calculateQCMetrics` function.
Outliers are identified within each cell type to avoid issues with systematic differences in the metrics between cell types.
The negative control wells do not contain any cells and are useful for quality control (as they _should_ manifest as outliers for the various metrics), but need to be removed prior to downstream analysis.


```r
library(scater)
sce.islam <- calculateQCMetrics(sce.islam)
sce.islam$grouping <- rep(c("mESC", "MEF", "Neg"), c(48, 44, 4))

libsize.drop <- isOutlier(sce.islam$total_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce.islam$grouping)
feature.drop <- isOutlier(sce.islam$total_features_by_counts, nmads=3, type="lower", 
    log=TRUE, batch=sce.islam$grouping)
spike.drop <- isOutlier(sce.islam$pct_counts_spike, nmads=3, type="higher", 
    batch=sce.islam$grouping)
    
sce.islam <- sce.islam[,!(libsize.drop | feature.drop | 
    spike.drop | sce.islam$grouping=="Neg")]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=ncol(sce.islam))
```

```
##   ByLibSize ByFeature BySpike Remaining
## 1         4         6      12        77
```

# Calculating spike-in size factors 

We apply the `computeSpikeFactors` method to estimate size factors for all cells.
This method computes the total count over all spike-in transcripts in each cell, and calculates size factors to equalize the total spike-in count across cells. 
Here, we set `general.use=TRUE` as we intend to apply the spike-in factors to all counts.


```r
library(scran)
sce.islam <- computeSpikeFactors(sce.islam, general.use=TRUE)
head(sizeFactors(sce.islam))    
```

```
## [1] 1.1486524 1.1274936 0.4285218 1.1011014 0.6646450 1.5932121
```

```r
head(sizeFactors(sce.islam, "spike")) # same as general size factors.
```

```
## [1] 1.1486524 1.1274936 0.4285218 1.1011014 0.6646450 1.5932121
```

Running `normalize` will use the spike-in-based size factors to compute normalized log-expression values.
Unlike the previous analyses, we do not have to define separate size factors for the spike-in transcripts.
This is because the relevant factors are already being used for all genes and spike-in transcripts when `general.use=TRUE`.
(The exception is if the experiment uses multiple spike-in sets that behave differently and need to be normalized separately.)


```r
sce.islam <- normalize(sce.islam)
```

For comparison, we also compute the deconvolution size factors for this data set [@lun2016pooling].


```r
deconv.sf <- computeSumFactors(sce.islam, sf.out=TRUE, 
    cluster=sce.islam$grouping)
head(deconv.sf)
```

```
## [1] 0.06424166 0.12172399 0.15547626 0.12876863 0.05680651 0.07888561
```

We observe a negative correlation between the two sets of size factors (Figure \@ref(fig:normplotspikemef)).
This is because MEFs contain more endogenous RNA, which reduces the relative spike-in coverage in each library (thereby decreasing the spike-in size factors) but increases the coverage of endogenous genes (thus increasing the deconvolution size factors).
If the spike-in size factors were applied to the counts, the expression values in MEFs would be scaled up while expression in mESCs would be scaled down.
However, the opposite would occur if deconvolution size factors were used.


```r
colours <- c(mESC="red", MEF="grey")
plot(sizeFactors(sce.islam), deconv.sf, col=colours[sce.islam$grouping], pch=16, 
    log="xy", xlab="Size factor (spike-in)", ylab="Size factor (deconvolution)")
legend("bottomleft", col=colours, legend=names(colours), pch=16)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-2-spike_files/figure-html/normplotspikemef-1.png" alt="Size factors from spike-in normalization, plotted against the size factors from deconvolution for all cells in the mESC/MEF dataset. Axes are shown on a log-scale, and cells are coloured according to their identity. Deconvolution size factors were computed with small pool sizes owing to the low number of cells of each type." width="100%" />
<p class="caption">(\#fig:normplotspikemef)Size factors from spike-in normalization, plotted against the size factors from deconvolution for all cells in the mESC/MEF dataset. Axes are shown on a log-scale, and cells are coloured according to their identity. Deconvolution size factors were computed with small pool sizes owing to the low number of cells of each type.</p>
</div>

Whether or not total RNA content is relevant -- and thus, the choice of normalization strategy -- depends on the biological hypothesis. 
In the HSC and brain analyses, variability in total RNA across the population was treated as noise and removed by non-DE normalization.
This may not always be appropriate if total RNA is associated with a biological difference of interest.
For example, @islam2011characterization observe a 5-fold difference in total RNA between mESCs and MEFs.
Similarly, the total RNA in a cell changes across phases of the cell cycle [@buettner2015computational].
Spike-in normalization will preserve these differences in total RNA content such that the corresponding biological groups can be easily resolved in downstream analyses.

__Comments from Aaron:__

- Setting `sf.out=TRUE` will directly return the size factors, rather than a `SingleCellExperiment` object containing those factors.
This is more convenient when only the size factors are required for further analysis.

# Concluding remarks 

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
## [14] GenomicRanges_1.35.1                  
## [15] GenomeInfoDb_1.19.0                   
## [16] IRanges_2.17.1                        
## [17] S4Vectors_0.21.3                      
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

