---
title: Elaborating on cell-based quality control in single-cell RNA-seq data
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
date: "2019-04-13"
vignette: >
  %\VignetteIndexEntry{06. Quality control details}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---





# Overview

Low-quality cells can often yield misleading results in downstream analyses, by:

- forming their own distinct cluster(s), complicating interpretation of the results.
This can be most obviously driven by increased mitochondrial proportions or enrichment for nuclear RNAs after cell damage.
However, very small libraries can also form their own clusters due to shifts in the mean upon log-transformation.
- containing genes that appear to be strongly "upregulated" due to the presence of very small size factors.
This is most problematic if some transcripts are present at constant levels in the ambient solution for all cells (i.e., wells or droplets).
Small counts will then be greatly inflated upon normalization with these size factors.
- containing genes that appear to be strongly "downregulated" due to the loss of RNA upon cell damage.
This seems most pronounced with ribosomal protein genes, though other cytoplasmic transcripts are likely to be affected.
- distorting the characterization of population heterogeneity during variance estimation or principal components analysis.
The first few principal components will capture differences in quality rather than biology, reducing the effectiveness of dimensionality reduction.
Similarly, genes with the largest variances will be driven by differences between low- and high-quality cells.

As such, we need to remove these cells at the start of the analysis.
Recall that we were defining low-quality cells as those with outlier values for various quality control (QC) metrics,
using the `isOutlier()` and `calculateQCMetrics()` functions from the *[scater](https://bioconductor.org/packages/3.9/scater)* package [@mccarthy2017scater].
Here, we will examine some of the reasoning behind the outlier-based QC in more detail.

# Assumptions of outlier identification 

An outlier-based definition for low-quality cells assumes that most cells are of high quality.
This is usually reasonable and can be experimentally supported in some situations by visually checking that the cells are intact, e.g., on the microwell plate.
Another assumption is that the QC metrics are independent on the biological state of each cell.
This ensures that any outlier values for these metrics are driven by technical factors rather than biological processes.
Thus, removing cells based on the metrics will not misrepresent the biology in downstream analyses.

The second assumption is most likely to be violated in highly heterogeneous cell populations.
For example, some cell types may naturally have less RNA or express fewer genes than other cell types.
Such cell types are more likely to be considered outliers and removed, even if they are of high quality.
The use of the MAD mitigates this problem by accounting for biological variability in the QC metrics.
A heterogeneous population should have higher variability in the metrics among high-quality cells, increasing the MAD and reducing the chance of incorrectly removing particular cell types (at the cost of reducing power to remove low-quality cells).
Nonetheless, filtering based on outliers may not be appropriate in extreme cases where one cell type is very different from the others.

Systematic differences in the QC metrics can be handled to some extent using the `batch=` argument in the `isOutlier()` function.
For example, setting `batch` to the plate of origin will identify outliers within each level of `batch`, using plate-specific median and MAD estimates.
This is obviously useful for accommodating known differences in experimental processing, e.g., sequencing at different depth or different amounts of added spike-in RNA. 
We can also include biological factors in `batch`, if those factors could result in systematically fewer expressed genes or lower RNA content.
However, this is not applicable in experiments where the factors are not known in advance.

# Checking for discarded cell types

## In the 416B data set 

We can diagnose loss of distinct cell types during QC by looking for differences in gene expression between the discarded and retained cells.
To demonstrate, we compute the average count across the discarded and retained pools in the 416B data set.


```r
library(SingleCellExperiment)
sce.full.416b <- readRDS("416B_preQC.rds")

library(scater)
lost <- calcAverage(counts(sce.full.416b)[,!sce.full.416b$PassQC])
kept <- calcAverage(counts(sce.full.416b)[,sce.full.416b$PassQC])
```

If the discarded pool is enriched for a certain cell type, we should observe increased expression of the corresponding marker genes.
No systematic upregulation of genes is apparent in the discarded pool in Figure \@ref(fig:discardplot416b), 
indicating that the QC step did not inadvertently filter out a cell type in the 416B dataset.


```r
# Avoid loss of points where either average is zero.
capped.lost <- pmax(lost, min(lost[lost>0]))
capped.kept <- pmax(kept, min(kept[kept>0]))

plot(capped.lost, capped.kept, xlab="Average count (discarded)", 
    ylab="Average count (retained)", log="xy", pch=16)
is.spike <- isSpike(sce.full.416b)
points(capped.lost[is.spike], capped.kept[is.spike], col="red", pch=16)
is.mito <- rowData(sce.full.416b)$is_feature_control_Mt
points(capped.lost[is.mito], capped.kept[is.mito], col="dodgerblue", pch=16)
```

<div class="figure">
<img src="qc_files/figure-html/discardplot416b-1.png" alt="Average counts across all discarded and retained cells in the 416B dataset. Each point represents a gene, with spike-in and mitochondrial transcripts in red and blue respectively." width="100%" />
<p class="caption">(\#fig:discardplot416b)Average counts across all discarded and retained cells in the 416B dataset. Each point represents a gene, with spike-in and mitochondrial transcripts in red and blue respectively.</p>
</div>

We examine this more closely by computing log-fold changes between the average counts of the two pools.
The `predFC` function stabilizes the log-fold change estimates by adding a prior count to the average of each pool.
We only examine the log-fold changes rather than formally testing for differential expression, as we are not interested in penalizing intra-pool heterogeneity.


```r
library(edgeR)
coefs <- predFC(cbind(lost, kept), design=cbind(1, c(1, 0)))[,2]
info <- data.frame(logFC=coefs, Lost=lost, Kept=kept, 
    row.names=rownames(sce.full.416b))
head(info[order(info$logFC, decreasing=TRUE),], 20)
```

```
##                       logFC      Lost        Kept
## ENSMUSG00000104647 6.844237  7.515034 0.000000000
## Nmur1              6.500909  5.909533 0.000000000
## Retn               6.250501 10.333172 0.196931018
## Fut9               6.096356  4.696897 0.010132032
## ENSMUSG00000102352 6.077614  9.393793 0.206637368
## ENSMUSG00000102379 6.029758  4.244690 0.000000000
## 1700101I11Rik      5.828821  4.483094 0.039199404
## Gm4952             5.698380  6.580862 0.172999108
## ENSMUSG00000106680 5.670156  3.389236 0.005234611
## ENSMUSG00000107955 5.554616  5.268508 0.132532601
## Gramd1c            5.446975  4.435342 0.103783669
## Jph3               5.361082  4.696897 0.139188080
## ENSMUSG00000092418 5.324462  3.395752 0.056931488
## 1700029I15Rik      5.316226  8.199510 0.394588776
## Pih1h3b            5.307439  2.546814 0.000000000
## ENSMUSG00000097176 5.275459  2.541927 0.003772842
## Olfr456            5.093383  2.186536 0.000000000
## ENSMUSG00000103731 5.017303  3.315016 0.107148909
## Klhdc8b            4.933215 19.861036 1.635081878
## ENSMUSG00000082449 4.881422  1.878759 0.000000000
```

Again, no obvious cell type markers are present in the top set of genes upregulated in the discarded pool.
The magnitude of the log-fold changes is less important, attributable to imprecision with few cells in the discarded pool.
Large log-fold changes can also be driven by enrichment or depletion of mitochondrial, ribosomal protein or nuclear genes upon cell damage.

## In the PBMC data set

For comparison, we consider the PBMC data set in which we previously identified a platelet population
(see the [previous workflow](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/tenx.html#marker-gene-detection)).
Recall that we relied on the use of the `emptyDrops()` method from the *[DropletUtils](https://bioconductor.org/packages/3.9/DropletUtils)* package to retain the platelets.
In contrast, if we had used a naive threshold on the total unique molecular identifier (UMI) count, we would have removed this population during the cell calling step.


```r
sce.pbmc <- readRDS("pbmc_data.rds")
wrong.keep <- sce.pbmc$total_counts >= 1000

lost <- calcAverage(counts(sce.pbmc)[,!wrong.keep])
kept <- calcAverage(counts(sce.pbmc)[,wrong.keep])
```

The presence of a distinct population in the discarded pool manifests in Figure \@ref(fig:discardplotpbmc) as a shift to the bottom-right for a number of genes.
This includes _PF4_, _PPBP_ and _SDPR_ that are strongly upregulated in the platelets.


```r
# Avoid loss of points where either average is zero.
capped.lost <- pmax(lost, min(lost[lost>0]))
capped.kept <- pmax(kept, min(kept[kept>0]))

plot(capped.lost, capped.kept, xlab="Average count (discarded)", 
    ylab="Average count (retained)", log="xy", pch=16)
platelet <- c("PF4", "PPBP", "SDPR")
points(capped.lost[platelet], capped.kept[platelet], col="orange", pch=16)
```

<div class="figure">
<img src="qc_files/figure-html/discardplotpbmc-1.png" alt="Average counts across all discarded and retained cells in the PBMC dataset, after using a more stringent filter on the total UMI count. Each point represents a gene, with platelet-related genes highlighted in orange." width="100%" />
<p class="caption">(\#fig:discardplotpbmc)Average counts across all discarded and retained cells in the PBMC dataset, after using a more stringent filter on the total UMI count. Each point represents a gene, with platelet-related genes highlighted in orange.</p>
</div>

These platelet-specific genes are also present among the top set of positive log-fold changes. 


```r
coefs <- predFC(cbind(lost, kept), design=cbind(1, c(1, 0)))[,2]
info <- data.frame(logFC=coefs, Lost=lost, Kept=kept, 
    row.names=rownames(sce.pbmc))
head(info[order(info$logFC, decreasing=TRUE),], 20)
```

```
##              logFC      Lost       Kept
## PF4       6.668019 4.3696131 0.17566713
## PPBP      6.548488 5.0020215 0.27309672
## HIST1H2AC 6.344409 3.2233183 0.14555421
## GNG11     6.218699 2.5968564 0.10122003
## SDPR      6.003830 2.2143681 0.09835611
## TUBB1     5.664169 1.7026339 0.09063991
## CLU       5.516478 1.3641616 0.05605478
## ACRBP     5.379732 1.2330620 0.05481902
## NRGN      5.172992 1.3592957 0.13054332
## RGS18     5.073244 1.7061891 0.25175703
## MAP3K7CL  4.953282 1.0281334 0.08985204
## SPARC     4.700745 0.6778018 0.02508210
## MMD       4.669388 0.7698856 0.06420673
## PGRMC1    4.563096 0.7579898 0.08169815
## CMTM5     4.220330 0.4618661 0.01657007
## TSC22D1   4.172856 0.5298668 0.05969562
## HRAT92    4.169966 0.4354403 0.01154032
## GP9       4.144189 0.4764224 0.03734688
## CTSA      4.039108 0.8570307 0.27227760
## F13A1     4.011858 0.4597065 0.05299097
```

## Avoiding loss of cell types

If cell types are being incorrectly discarded, the most direct solution is to relax the QC filters by increasing `nmads=` in the `isOutlier()` calls.
We can also avoid filtering on metrics that are associated with genuine biological differences between cell types.
The most extreme approach would be to not perform any QC filtering at all, thus guaranteeing that all cell types in the data are retained.
However, this obviously comes with an increased risk of retaining more low-quality damaged cells.
Such cells will cause problems in downstream analyses as discussed above, which motivates the use of a more strict filter (at least on the first pass) in our workflows.

As an aside, it is worth mentioning that the true technical quality of a cell may be correlated with its type.
(This differs from a correlation between the cell type and the QC metrics, as the latter are our imperfect proxies for quality.)
This can arise if some cell types are not amenable to dissociation or microfluidics handling during the scRNA-seq protocol.
In such cases, it is possible to correctly discard an entire cell type during QC if all of its members are damaged.
Indeed, concerns over the computational removal of cell types during QC are probably minor compared to losses in the experimental protocol.

# Alternative approaches to quality control

## Using fixed thresholds

One alternative strategy is to set pre-defined thresholds on each QC metric.
For example, we might remove all cells with library sizes below 100000 and numbers of expressed genes below 4000.
This avoids any assumptions associated with the use of outliers to identify low-quality cells.
However, it generally requires considerable experience to determine appropriate thresholds for each experimental protocol and biological system.
For example, thresholds for read count-based data are simply not applicable for UMI-based data, and vice versa.
Indeed, even with the same protocol and system, the appropriate threshold can vary from run to run due to the vagaries of RNA capture and sequencing.

## Using PCA-based outliers

Another strategy is to perform a principal components analysis (PCA) based on the quality metrics for each cell, e.g., the total number of reads, the total number of features and the proportion of mitochondrial or spike-in reads.
Outliers on a PCA plot may be indicative of low-quality cells that have aberrant technical properties compared to the (presumed) majority of high-quality cells.
This is demonstrated below on a brain cell dataset from @tasic2016adult, using functions from *[scater](https://bioconductor.org/packages/3.9/scater)*.


```r
# Obtaining the dataset.
library(scRNAseq)
data(allen)

# Setting up the data.
sce.allen <- as(allen, "SingleCellExperiment")
assayNames(sce.allen) <- "counts"
isSpike(sce.allen, "ERCC") <- grep("ERCC", rownames(sce.allen))

# Computing the QC metrics and running PCA.
library(scater)
sce.allen <- calculateQCMetrics(sce.allen)
sce.allen <- runPCA(sce.allen, use_coldata=TRUE, detect_outliers=TRUE)
table(sce.allen$outlier)
```

```
## 
## FALSE  TRUE 
##   374     5
```

Methods like PCA-based outlier detection and support vector machines can provide more power to distinguish low-quality cells from high-quality counterparts [@ilicic2016classification].
This is because they are able to detect subtle patterns across many quality metrics simultaneously. 
However, this comes at some cost to interpretability, as the reason for removing a given cell may not always be obvious.
Users interested in the more sophisticated approaches are referred to the *[scater](https://bioconductor.org/packages/3.9/scater)* and *[cellity](https://bioconductor.org/packages/3.9/cellity)* packages.

## Using the gene expression profiles

For completeness, we note that outliers can also be identified from the gene expression profiles, rather than QC metrics.
We consider this to be a risky strategy as it can remove high-quality cells in rare populations.
Even if subpopulations are explicitly captured with a mixture model, removal of outlier cells will simply reinforce the existing model.
This may be misleading if it understates the biological heterogeneity in each population.

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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] scRNAseq_1.9.0              edgeR_3.25.4               
##  [3] limma_3.39.15               scater_1.11.16             
##  [5] ggplot2_3.1.1               SingleCellExperiment_1.5.2 
##  [7] SummarizedExperiment_1.13.0 DelayedArray_0.9.9         
##  [9] BiocParallel_1.17.18        matrixStats_0.54.0         
## [11] Biobase_2.43.1              GenomicRanges_1.35.1       
## [13] GenomeInfoDb_1.19.3         IRanges_2.17.4             
## [15] S4Vectors_0.21.23           BiocGenerics_0.29.2        
## [17] knitr_1.22                  BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0         colorspace_1.4-1        
##   [3] mvoutlier_2.0.9          modeltools_0.2-22       
##   [5] class_7.3-15             rio_0.5.16              
##   [7] mclust_5.4.3             XVector_0.23.2          
##   [9] pls_2.7-1                BiocNeighbors_1.1.13    
##  [11] cvTools_0.3.2            flexmix_2.3-15          
##  [13] mvtnorm_1.0-10           ranger_0.11.2           
##  [15] codetools_0.2-16         splines_3.7.0           
##  [17] sROC_0.1-2               robustbase_0.93-4       
##  [19] robCompositions_2.0.10   kernlab_0.9-27          
##  [21] cluster_2.0.8            BiocManager_1.30.4      
##  [23] rrcov_1.4-7              compiler_3.7.0          
##  [25] assertthat_0.2.1         Matrix_1.2-17           
##  [27] lazyeval_0.2.2           BiocSingular_0.99.16    
##  [29] htmltools_0.3.6          tools_3.7.0             
##  [31] rsvd_1.0.0               gtable_0.3.0            
##  [33] glue_1.3.1               GenomeInfoDbData_1.2.1  
##  [35] dplyr_0.8.0.1            Rcpp_1.0.1              
##  [37] carData_3.0-2            trimcluster_0.1-2.1     
##  [39] cellranger_1.1.0         zCompositions_1.2.0     
##  [41] sgeostat_1.0-27          fpc_2.1-11.1            
##  [43] DelayedMatrixStats_1.5.2 lmtest_0.9-36           
##  [45] xfun_0.6                 laeken_0.5.0            
##  [47] stringr_1.4.0            ps_1.3.0                
##  [49] openxlsx_4.1.0           irlba_2.3.3             
##  [51] DEoptimR_1.0-8           zlibbioc_1.29.0         
##  [53] MASS_7.3-51.4            zoo_1.8-5               
##  [55] scales_1.0.0             VIM_4.8.0               
##  [57] hms_0.4.2                RColorBrewer_1.1-2      
##  [59] yaml_2.2.0               curl_3.3                
##  [61] NADA_1.6-1               gridExtra_2.3           
##  [63] reshape_0.8.8            stringi_1.4.3           
##  [65] pcaPP_1.9-73             simpleSingleCell_1.7.20 
##  [67] e1071_1.7-1              boot_1.3-22             
##  [69] zip_2.0.1                truncnorm_1.0-8         
##  [71] prabclus_2.2-7           rlang_0.3.4             
##  [73] pkgconfig_2.0.2          bitops_1.0-6            
##  [75] evaluate_0.13            lattice_0.20-38         
##  [77] purrr_0.3.2              processx_3.3.0          
##  [79] tidyselect_0.2.5         GGally_1.4.0            
##  [81] plyr_1.8.4               magrittr_1.5            
##  [83] bookdown_0.9             R6_2.4.0                
##  [85] pillar_1.3.1             haven_2.1.0             
##  [87] foreign_0.8-71           withr_2.1.2             
##  [89] survival_2.44-1.1        abind_1.4-5             
##  [91] RCurl_1.95-4.12          sp_1.3-1                
##  [93] nnet_7.3-12              tibble_2.1.1            
##  [95] crayon_1.3.4             car_3.0-2               
##  [97] rmarkdown_1.12           viridis_0.5.1           
##  [99] locfit_1.5-9.1           grid_3.7.0              
## [101] readxl_1.3.1             data.table_1.12.2       
## [103] callr_3.2.0              forcats_0.4.0           
## [105] diptest_0.75-7           vcd_1.4-4               
## [107] digest_0.6.18            munsell_0.5.0           
## [109] beeswarm_0.2.3           viridisLite_0.3.0       
## [111] vipor_0.4.5
```

# References

