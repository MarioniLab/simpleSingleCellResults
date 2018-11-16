---
title: Further strategies for analyzing single-cell RNA-seq data
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
date: "2018-11-16"
vignette: >
  %\VignetteIndexEntry{12. Further analysis strategies}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---



# Overview

Here, we describe a few additional analyses that can be performed with single-cell RNA sequencing data.
This includes detection of significant correlations between genes
and regressing out the effect of cell cycle from the gene expression matrix.

# Identifying correlated gene pairs with Spearman's rho

scRNA-seq data is commonly used to identify correlations between the expression profiles of different genes.
This is quantified by computing Spearman's rho, which accommodates non-linear relationships in the expression values.
Non-zero correlations between pairs of genes provide evidence for their co-regulation.
However, the noise in the data requires some statistical analysis to determine whether a correlation is significantly non-zero.

To demonstrate, we use the `correlatePairs` function to identify significant correlations between the various histocompatability antigens in the haematopoietic stem cell (HSC) dataset [@wilson2015combined].
The significance of each correlation is determined using a permutation test.
For each pair of genes, the null hypothesis is that the expression profiles of two genes are independent.
Shuffling the profiles and recalculating the correlation yields a null distribution that is used to obtain a _p_-value for each observed correlation value [@phipson2010permutation].


```r
library(scran)
sce.hsc <- readRDS("hsc_data.rds")    

set.seed(100)
var.cor <- correlatePairs(sce.hsc, subset.row=grep("^H2-", rownames(sce.hsc)))
head(var.cor)
```

```
## DataFrame with 6 rows and 6 columns
##         gene1       gene2               rho              p.value
##   <character> <character>         <numeric>            <numeric>
## 1      H2-Ab1      H2-Eb1 0.497634203104049   1.999998000002e-06
## 2       H2-Aa      H2-Ab1 0.488479262672811   1.999998000002e-06
## 3       H2-D1       H2-K1 0.412280566558266  3.3999966000034e-05
## 4       H2-Aa      H2-Eb1  0.41029237242421  3.7999962000038e-05
## 5      H2-Ab1     H2-DMb1 0.359662777615092 0.000455999544000456
## 6       H2-Q6       H2-Q7 0.339981196923693  0.00098999901000099
##                    FDR   limited
##              <numeric> <logical>
## 1 0.000434999565000435      TRUE
## 2 0.000434999565000435      TRUE
## 3  0.00413249586750413     FALSE
## 4  0.00413249586750413     FALSE
## 5   0.0396719603280397     FALSE
## 6   0.0717749282250718     FALSE
```

Correction for multiple testing across many gene pairs is performed by controlling the FDR at 5%.


```r
sig.cor <- var.cor$FDR <= 0.05
summary(sig.cor)
```

```
##    Mode   FALSE    TRUE 
## logical     430       5
```

We can also compute correlations between specific pairs of genes, or between all pairs between two distinct sets of genes.
The example below computes the correlation between _Fos_ and _Jun_, which dimerize to form the AP-1 transcription factor [@angel1991role].


```r
correlatePairs(sce.hsc, subset.row=cbind("Fos", "Jun"))
```

```
## DataFrame with 1 row and 6 columns
##         gene1       gene2               rho            p.value
##   <character> <character>         <numeric>          <numeric>
## 1         Fos         Jun 0.466855724920241 1.999998000002e-06
##                  FDR   limited
##            <numeric> <logical>
## 1 1.999998000002e-06      TRUE
```

Examination of the expression profiles in Figure \@ref(fig:fosjuncorplot) confirms the presence of a modest correlation between these two genes.


```r
library(scater)
plotExpression(sce.hsc, features="Fos", x="Jun")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-4-misc_files/figure-html/fosjuncorplot-1.png" alt="Expression of _Fos_ plotted against the expression of _Jun_ for all cells in the HSC dataset." width="100%" />
<p class="caption">(\#fig:fosjuncorplot)Expression of _Fos_ plotted against the expression of _Jun_ for all cells in the HSC dataset.</p>
</div>

The use of `correlatePairs` is primarily intended to identify correlated gene pairs for validation studies.
Obviously, non-zero correlations do not provide evidence for a direct regulatory interaction, let alone specify causality.
To construct regulatory networks involving many genes, we suggest using dedicated packages such as *[WCGNA](https://CRAN.R-project.org/package=WCGNA)*.

__Comments from Aaron:__

- We suggest only computing correlations between a subset of genes of interest, known either _a priori_ or empirically defined, e.g., as HVGs.
Computing correlations across all genes will take too long; unnecessarily increase the severity of the multiple testing correction; 
and may prioritize strong but uninteresting correlations, e.g., between tightly co-regulated house-keeping genes.
- The `correlatePairs` function can also return gene-centric output by setting `per.gene=TRUE`.
This calculates a combined _p_-value [@simes1986improved] for each gene that indicates whether it is significantly correlated to any other gene.
From a statistical perspective, this is a more natural approach to correcting for multiple testing when genes, rather than pairs of genes, are of interest.
- The `Limited` field indicates whether the _p_-value was lower-bounded by the number of permutations.
If this is `TRUE` for any non-significant gene at the chosen FDR threshold, consider increasing the number of permutations to improve power.

# Blocking on the cell cycle phase

Cell cycle phase is usually uninteresting in studies focusing on other aspects of biology.
However, the effects of cell cycle on the expression profile can mask other effects and interfere with the interpretation of the results.
This cannot be avoided by simply removing cell cycle marker genes, as the cell cycle can affect a substantial number of other transcripts [@buettner2015computational].
Rather, more sophisticated strategies are required, one of which is demonstrated below using data from a study of T Helper 2 (T~H~2) cells [@mahata2014singlecell].


```r
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
mahata.fname <- bfcrpath(bfc, 
    "http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S7.xlsx")
```

@buettner2015computational have already applied quality control and normalized the data, so we can use them directly as log-expression values (accessible as Supplementary Data 1 of https://dx.doi.org/10.1038/nbt.3102).


```r
library(readxl)
incoming <- as.data.frame(read_excel(mahata.fname, sheet=1))
rownames(incoming) <- incoming[,1]
incoming <- incoming[,-1]
incoming <- incoming[,!duplicated(colnames(incoming))] # Remove duplicated genes.
sce.th2 <- SingleCellExperiment(list(logcounts=t(incoming)))
```

We empirically identify the cell cycle phase using the pair-based classifier in `cyclone`.
The majority of cells in Figure \@ref(fig:phaseplotth2) seem to lie in G1 phase, with small numbers of cells in the other phases.


```r
library(org.Mm.eg.db)
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce.th2), keytype="SYMBOL", column="ENSEMBL")

set.seed(100)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
    package="scran"))
assignments <- cyclone(sce.th2, mm.pairs, gene.names=ensembl, assay.type="logcounts")

plot(assignments$score$G1, assignments$score$G2M, 
    xlab="G1 score", ylab="G2/M score", pch=16)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-4-misc_files/figure-html/phaseplotth2-1.png" alt="Cell cycle phase scores from applying the pair-based classifier on the T~H~2 dataset, where each point represents a cell." width="100%" />
<p class="caption">(\#fig:phaseplotth2)Cell cycle phase scores from applying the pair-based classifier on the T~H~2 dataset, where each point represents a cell.</p>
</div>

We can block directly on the phase scores in downstream analyses.
This is more graduated than using a strict assignment of each cell to a specific phase, as the magnitude of the score considers the uncertainty of the assignment.
The phase covariates in the design matrix will absorb any phase-related effects on expression such that they will not affect estimation of the effects of other experimental factors.
Users should also ensure that the phase score is not confounded with other factors of interest.
For example, model fitting is not possible if all cells in one experimental condition are in one phase, and all cells in another condition are in a different phase.


```r
design <- model.matrix(~ G1 + G2M, assignments$score)
fit.block <- trendVar(sce.th2, design=design, parametric=TRUE, use.spikes=NA)
dec.block <- decomposeVar(sce.th2, fit.block)

library(limma)
sce.th2.block <- sce.th2
assay(sce.th2.block, "corrected") <- removeBatchEffect(
    logcounts(sce.th2), covariates=design[,-1])

sce.th2.block <- denoisePCA(sce.th2.block, technical=dec.block, 
    assay.type="corrected")
dim(reducedDim(sce.th2.block, "PCA"))
```

```
## [1] 81  5
```

The result of blocking on `design` is visualized with some PCA plots in Figure \@ref(fig:pcaplotth2).
Before removal, the distribution of cells along the first two principal components is strongly associated with their G1 and G2/M scores.
This is no longer the case after removal, which suggests that the cell cycle effect has been mitigated.


```r
sce.th2$G1score <- sce.th2.block$G1score <- assignments$score$G1
sce.th2$G2Mscore <- sce.th2.block$G2Mscore <- assignments$score$G2M

# Without blocking on phase score.
fit <- trendVar(sce.th2, parametric=TRUE, use.spikes=NA) 
sce.th2 <- denoisePCA(sce.th2, technical=fit$trend)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
out <- plotReducedDim(sce.th2, use_dimred="PCA", ncomponents=2, colour_by="G1score", 
    size_by="G2Mscore") + fontsize + ggtitle("Before removal")

# After blocking on the phase score.
out2 <- plotReducedDim(sce.th2.block, use_dimred="PCA", ncomponents=2, 
    colour_by="G1score", size_by="G2Mscore") + fontsize + 
    ggtitle("After removal")
multiplot(out, out2, cols=2)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-4-misc_files/figure-html/pcaplotth2-1.png" alt="PCA plots before (left) and after (right) removal of the cell cycle effect in the T~H~2 dataset. Each cell is represented by a point with colour and size determined by the G1 and G2/M scores, respectively." width="1152" />
<p class="caption">(\#fig:pcaplotth2)PCA plots before (left) and after (right) removal of the cell cycle effect in the T~H~2 dataset. Each cell is represented by a point with colour and size determined by the G1 and G2/M scores, respectively.</p>
</div>

As an aside, this dataset contains cells at various stages of differentiation [@mahata2014singlecell].
This is an ideal use case for diffusion maps which perform dimensionality reduction along a continuous process.
In Figure \@ref(fig:diffusionth2), cells are arranged along a trajectory in the low-dimensional space.
The first diffusion component is likely to correspond to T~H~2 differentiation, given that a key regulator _Gata3_ [@zhu2006gata3] changes in expression from left to right.


```r
plotDiffusionMap(sce.th2.block, colour_by="Gata3",
    run_args=list(use_dimred="PCA", sigma=25)) + fontsize
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-4-misc_files/figure-html/diffusionth2-1.png" alt="A diffusion map for the T~H~2 dataset, where each cell is coloured by its expression of _Gata3_. A larger `sigma` is used compared to the default value to obtain a smoother plot." width="100%" />
<p class="caption">(\#fig:diffusionth2)A diffusion map for the T~H~2 dataset, where each cell is coloured by its expression of _Gata3_. A larger `sigma` is used compared to the default value to obtain a smoother plot.</p>
</div>

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
##  [1] readxl_1.1.0                          
##  [2] DelayedMatrixStats_1.5.0              
##  [3] gdata_2.18.0                          
##  [4] R.utils_2.7.0                         
##  [5] R.oo_1.22.0                           
##  [6] R.methodsS3_1.7.1                     
##  [7] scRNAseq_1.9.0                        
##  [8] edgeR_3.25.0                          
##  [9] Matrix_1.2-15                         
## [10] org.Hs.eg.db_3.7.0                    
## [11] EnsDb.Hsapiens.v86_2.99.0             
## [12] ensembldb_2.7.2                       
## [13] AnnotationFilter_1.7.0                
## [14] DropletUtils_1.3.1                    
## [15] pheatmap_1.0.10                       
## [16] cluster_2.0.7-1                       
## [17] dynamicTreeCut_1.63-1                 
## [18] limma_3.39.1                          
## [19] scran_1.11.4                          
## [20] scater_1.11.2                         
## [21] ggplot2_3.1.0                         
## [22] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0
## [23] GenomicFeatures_1.35.1                
## [24] org.Mm.eg.db_3.7.0                    
## [25] AnnotationDbi_1.45.0                  
## [26] SingleCellExperiment_1.5.0            
## [27] SummarizedExperiment_1.13.0           
## [28] DelayedArray_0.9.0                    
## [29] BiocParallel_1.17.1                   
## [30] matrixStats_0.54.0                    
## [31] Biobase_2.43.0                        
## [32] GenomicRanges_1.35.1                  
## [33] GenomeInfoDb_1.19.1                   
## [34] IRanges_2.17.1                        
## [35] S4Vectors_0.21.4                      
## [36] BiocGenerics_0.29.1                   
## [37] bindrcpp_0.2.2                        
## [38] BiocFileCache_1.7.0                   
## [39] dbplyr_1.2.2                          
## [40] knitr_1.20                            
## [41] BiocStyle_2.11.0                      
## 
## loaded via a namespace (and not attached):
##   [1] tidyselect_0.2.5         RSQLite_2.1.1           
##   [3] grid_3.6.0               trimcluster_0.1-2.1     
##   [5] Rtsne_0.15               munsell_0.5.0           
##   [7] destiny_2.13.0           statmod_1.4.30          
##   [9] sROC_0.1-2               withr_2.1.2             
##  [11] colorspace_1.3-2         highr_0.7               
##  [13] robustbase_0.93-3        vcd_1.4-4               
##  [15] VIM_4.7.0                TTR_0.23-4              
##  [17] labeling_0.3             GenomeInfoDbData_1.2.0  
##  [19] cvTools_0.3.2            bit64_0.9-7             
##  [21] rhdf5_2.27.1             rprojroot_1.3-2         
##  [23] xfun_0.4                 ggthemes_4.0.1          
##  [25] diptest_0.75-7           R6_2.3.0                
##  [27] ggbeeswarm_0.6.0         robCompositions_2.0.9   
##  [29] RcppEigen_0.3.3.4.0      locfit_1.5-9.1          
##  [31] mvoutlier_2.0.9          flexmix_2.3-14          
##  [33] bitops_1.0-6             reshape_0.8.8           
##  [35] assertthat_0.2.0         scales_1.0.0            
##  [37] nnet_7.3-12              beeswarm_0.2.3          
##  [39] gtable_0.2.0             rlang_0.3.0.1           
##  [41] scatterplot3d_0.3-41     splines_3.6.0           
##  [43] rtracklayer_1.43.0       lazyeval_0.2.1          
##  [45] BiocManager_1.30.4       yaml_2.2.0              
##  [47] reshape2_1.4.3           abind_1.4-5             
##  [49] backports_1.1.2          tools_3.6.0             
##  [51] bookdown_0.7             zCompositions_1.1.2     
##  [53] RColorBrewer_1.1-2       proxy_0.4-22            
##  [55] Rcpp_1.0.0               plyr_1.8.4              
##  [57] progress_1.2.0           zlibbioc_1.29.0         
##  [59] purrr_0.2.5              RCurl_1.95-4.11         
##  [61] prettyunits_1.0.2        viridis_0.5.1           
##  [63] cowplot_0.9.3            zoo_1.8-4               
##  [65] haven_1.1.2              magrittr_1.5            
##  [67] data.table_1.11.8        openxlsx_4.1.0          
##  [69] lmtest_0.9-36            truncnorm_1.0-8         
##  [71] mvtnorm_1.0-8            ProtGenerics_1.15.0     
##  [73] hms_0.4.2                evaluate_0.12           
##  [75] smoother_1.1             XML_3.98-1.16           
##  [77] rio_0.5.10               mclust_5.4.1            
##  [79] gridExtra_2.3            compiler_3.6.0          
##  [81] biomaRt_2.39.2           tibble_1.4.2            
##  [83] KernSmooth_2.23-15       crayon_1.3.4            
##  [85] htmltools_0.3.6          pcaPP_1.9-73            
##  [87] rrcov_1.4-4              DBI_1.0.0               
##  [89] MASS_7.3-51.1            fpc_2.1-11.1            
##  [91] rappdirs_0.3.1           boot_1.3-20             
##  [93] car_3.0-2                sgeostat_1.0-27         
##  [95] bindr_0.1.1              igraph_1.2.2            
##  [97] forcats_0.3.0            pkgconfig_2.0.2         
##  [99] GenomicAlignments_1.19.0 foreign_0.8-71          
## [101] laeken_0.4.6             sp_1.3-1                
## [103] vipor_0.4.5              XVector_0.23.0          
## [105] NADA_1.6-1               stringr_1.3.1           
## [107] digest_0.6.18            pls_2.7-0               
## [109] Biostrings_2.51.1        rmarkdown_1.10          
## [111] cellranger_1.1.0         curl_3.2                
## [113] kernlab_0.9-27           gtools_3.8.1            
## [115] Rsamtools_1.35.0         modeltools_0.2-22       
## [117] Rhdf5lib_1.5.0           carData_3.0-2           
## [119] BiocNeighbors_1.1.1      viridisLite_0.3.0       
## [121] pillar_1.3.0             lattice_0.20-38         
## [123] GGally_1.4.0             httr_1.3.1              
## [125] DEoptimR_1.0-8           survival_2.43-1         
## [127] xts_0.11-2               glue_1.3.0              
## [129] zip_1.0.0                prabclus_2.2-6          
## [131] bit_1.1-14               class_7.3-14            
## [133] stringi_1.2.4            HDF5Array_1.11.0        
## [135] blob_1.1.1               memoise_1.1.0           
## [137] dplyr_0.7.8              irlba_2.3.3             
## [139] e1071_1.7-0
```

# References

