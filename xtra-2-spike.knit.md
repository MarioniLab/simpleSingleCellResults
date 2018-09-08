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
date: "2018-09-08"
vignette: >
  %\VignetteIndexEntry{7. Spike-in normalization of scRNA-seq data}
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
deconv.sf <- computeSumFactors(sce.islam, sf.out=TRUE, cluster=sce.islam$grouping)
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

- We only use genes with average counts greater than 1 (as specified in `min.mean`) to compute the deconvolution size factors.
This avoids problems with discreteness as mentioned in our previous uses of `computeSumFactors`.
- Setting `sf.out=TRUE` will directly return the size factors, rather than a `SingleCellExperiment` object containing those factors.
This is more convenient when only the size factors are required for further analysis.

# Concluding remarks 

All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org).
The specific version numbers of the packages used are shown below, along with the version of the R installation.


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
##  [1] scran_1.9.20                bindrcpp_0.2.2             
##  [3] BiocFileCache_1.5.5         dbplyr_1.2.2               
##  [5] scRNAseq_1.7.0              edgeR_3.23.3               
##  [7] limma_3.37.4                scater_1.9.20              
##  [9] ggplot2_3.0.0               SingleCellExperiment_1.3.10
## [11] SummarizedExperiment_1.11.6 DelayedArray_0.7.37        
## [13] BiocParallel_1.15.11        matrixStats_0.54.0         
## [15] Biobase_2.41.2              GenomicRanges_1.33.13      
## [17] GenomeInfoDb_1.17.1         IRanges_2.15.17            
## [19] S4Vectors_0.19.19           BiocGenerics_0.27.1        
## [21] knitr_1.20                  BiocStyle_2.9.6            
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0         colorspace_1.3-2        
##   [3] mvoutlier_2.0.9          class_7.3-14            
##   [5] modeltools_0.2-22        rio_0.5.10              
##   [7] dynamicTreeCut_1.63-1    mclust_5.4.1            
##   [9] rprojroot_1.3-2          XVector_0.21.3          
##  [11] pls_2.7-0                cvTools_0.3.2           
##  [13] bit64_0.9-7              flexmix_2.3-14          
##  [15] mvtnorm_1.0-8            splines_3.5.0           
##  [17] sROC_0.1-2               robustbase_0.93-2       
##  [19] robCompositions_2.0.8    kernlab_0.9-27          
##  [21] cluster_2.0.7-1          HDF5Array_1.9.15        
##  [23] httr_1.3.1               BiocManager_1.30.2      
##  [25] rrcov_1.4-4              compiler_3.5.0          
##  [27] backports_1.1.2          assertthat_0.2.0        
##  [29] Matrix_1.2-14            lazyeval_0.2.1          
##  [31] htmltools_0.3.6          tools_3.5.0             
##  [33] igraph_1.2.2             gtable_0.2.0            
##  [35] glue_1.3.0               GenomeInfoDbData_1.1.0  
##  [37] reshape2_1.4.3           dplyr_0.7.6             
##  [39] rappdirs_0.3.1           Rcpp_0.12.18            
##  [41] carData_3.0-1            trimcluster_0.1-2.1     
##  [43] cellranger_1.1.0         zCompositions_1.1.1     
##  [45] sgeostat_1.0-27          kmknn_0.99.16           
##  [47] fpc_2.1-11.1             DelayedMatrixStats_1.3.8
##  [49] lmtest_0.9-36            xfun_0.3                
##  [51] laeken_0.4.6             stringr_1.3.1           
##  [53] openxlsx_4.1.0           statmod_1.4.30          
##  [55] DEoptimR_1.0-8           zlibbioc_1.27.0         
##  [57] MASS_7.3-50              zoo_1.8-3               
##  [59] scales_1.0.0             VIM_4.7.0               
##  [61] hms_0.4.2                rhdf5_2.25.9            
##  [63] RColorBrewer_1.1-2       yaml_2.2.0              
##  [65] curl_3.2                 memoise_1.1.0           
##  [67] NADA_1.6-1               gridExtra_2.3           
##  [69] RSQLite_2.1.1            reshape_0.8.7           
##  [71] stringi_1.2.4            highr_0.7               
##  [73] pcaPP_1.9-73             e1071_1.7-0             
##  [75] boot_1.3-20              zip_1.0.0               
##  [77] truncnorm_1.0-8          prabclus_2.2-6          
##  [79] rlang_0.2.2              pkgconfig_2.0.2         
##  [81] bitops_1.0-6             evaluate_0.11           
##  [83] lattice_0.20-35          purrr_0.2.5             
##  [85] Rhdf5lib_1.3.3           bindr_0.1.1             
##  [87] bit_1.1-14               tidyselect_0.2.4        
##  [89] GGally_1.4.0             plyr_1.8.4              
##  [91] magrittr_1.5             bookdown_0.7            
##  [93] R6_2.2.2                 DBI_1.0.0               
##  [95] pillar_1.3.0             haven_1.1.2             
##  [97] foreign_0.8-71           withr_2.1.2             
##  [99] survival_2.42-6          abind_1.4-5             
## [101] RCurl_1.95-4.11          sp_1.3-1                
## [103] nnet_7.3-12              tibble_1.4.2            
## [105] crayon_1.3.4             car_3.0-2               
## [107] rmarkdown_1.10           viridis_0.5.1           
## [109] locfit_1.5-9.1           grid_3.5.0              
## [111] readxl_1.1.0             data.table_1.11.4       
## [113] blob_1.1.1               forcats_0.3.0           
## [115] diptest_0.75-7           vcd_1.4-4               
## [117] digest_0.6.16            munsell_0.5.0           
## [119] beeswarm_0.2.3           viridisLite_0.3.0       
## [121] vipor_0.4.5
```

# References
