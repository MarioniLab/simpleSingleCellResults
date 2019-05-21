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
date: "2019-05-20"
vignette: >
  %\VignetteIndexEntry{13. Further analysis strategies}
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
## 1       H2-Aa      H2-Ab1 0.537721745295378   1.999998000002e-06
## 2      H2-Ab1      H2-Eb1 0.427978052802737  1.7999982000018e-05
## 3       H2-D1       H2-K1 0.412280566558266  4.1999958000042e-05
## 4       H2-Aa      H2-Eb1 0.407417967726523  5.1999948000052e-05
## 5       H2-Q6       H2-Q7 0.339981196923693 0.000943999056000944
## 6       H2-K2       H2-Q6 0.308678698580522  0.00280599719400281
##                   FDR   limited
##             <numeric> <logical>
## 1 0.00086999913000087      TRUE
## 2 0.00391499608500391     FALSE
## 3 0.00565499434500566     FALSE
## 4 0.00565499434500566     FALSE
## 5  0.0821279178720821     FALSE
## 6   0.203434796565203     FALSE
```

Correction for multiple testing across many gene pairs is performed by controlling the FDR at 5%.


```r
sig.cor <- var.cor$FDR <= 0.05
summary(sig.cor)
```

```
##    Mode   FALSE    TRUE 
## logical     431       4
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
## 1         Fos         Jun 0.467634048980473 7.999992000008e-06
##                  FDR   limited
##            <numeric> <logical>
## 1 7.999992000008e-06     FALSE
```

Examination of the expression profiles in Figure \@ref(fig:fosjuncorplot) confirms the presence of a modest correlation between these two genes.


```r
library(scater)
plotExpression(sce.hsc, features="Fos", x="Jun")
```

<div class="figure">
<img src="misc_files/figure-html/fosjuncorplot-1.png" alt="Expression of _Fos_ plotted against the expression of _Jun_ for all cells in the HSC dataset." width="100%" />
<p class="caption">(\#fig:fosjuncorplot)Expression of _Fos_ plotted against the expression of _Jun_ for all cells in the HSC dataset.</p>
</div>

The use of `correlatePairs` is primarily intended to identify correlated gene pairs for validation studies.
Obviously, non-zero correlations do not provide evidence for a direct regulatory interaction, let alone specify causality.
To construct regulatory networks involving many genes, we suggest using dedicated packages such as *[WCGNA](https://CRAN.R-project.org/package=WCGNA)*.

__Comments from Aaron:__

- We suggest only computing correlations between a subset of genes of interest, known either _a priori_ or empirically defined, e.g., as HVGs.
Computing correlations across all genes will take too long; unnecessarily increase the severity of the multiple testing correction; 
and may prioritize strong but uninteresting correlations, e.g., between tightly co-regulated house-keeping genes.
- The `correlateGenes()` function can be used on the output of `correlatePairs()` to return gene-centric output. 
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
<img src="misc_files/figure-html/phaseplotth2-1.png" alt="Cell cycle phase scores from applying the pair-based classifier on the T~H~2 dataset, where each point represents a cell." width="100%" />
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
<img src="misc_files/figure-html/pcaplotth2-1.png" alt="PCA plots before (left) and after (right) removal of the cell cycle effect in the T~H~2 dataset. Each cell is represented by a point with colour and size determined by the G1 and G2/M scores, respectively." width="1152" />
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
<img src="misc_files/figure-html/diffusionth2-1.png" alt="A diffusion map for the T~H~2 dataset, where each cell is coloured by its expression of _Gata3_. A larger `sigma` is used compared to the default value to obtain a smoother plot." width="100%" />
<p class="caption">(\#fig:diffusionth2)A diffusion map for the T~H~2 dataset, where each cell is coloured by its expression of _Gata3_. A larger `sigma` is used compared to the default value to obtain a smoother plot.</p>
</div>

# Concluding remarks 

All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org).
The specific version numbers of the packages used are shown below, along with the version of the R installation.


```r
sessionInfo()
```

```
## R version 3.6.0 Patched (2019-05-02 r76458)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.2 LTS
## 
## Matrix products: default
## BLAS:   /home/luna/Software/R/R-3-6-branch-dev/lib/libRblas.so
## LAPACK: /home/luna/Software/R/R-3-6-branch-dev/lib/libRlapack.so
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
##  [1] limma_3.41.2                org.Mm.eg.db_3.8.2         
##  [3] AnnotationDbi_1.47.0        readxl_1.3.1               
##  [5] BiocFileCache_1.9.0         dbplyr_1.4.0               
##  [7] scater_1.13.3               ggplot2_3.1.1              
##  [9] scran_1.13.3                SingleCellExperiment_1.7.0 
## [11] SummarizedExperiment_1.15.1 DelayedArray_0.11.0        
## [13] BiocParallel_1.19.0         matrixStats_0.54.0         
## [15] Biobase_2.45.0              GenomicRanges_1.37.4       
## [17] GenomeInfoDb_1.21.1         IRanges_2.19.3             
## [19] S4Vectors_0.23.3            BiocGenerics_0.31.2        
## [21] knitr_1.23                  BiocStyle_2.13.0           
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0         colorspace_1.4-1        
##   [3] RcppEigen_0.3.3.5.0      class_7.3-15            
##   [5] rio_0.5.16               dynamicTreeCut_1.63-1   
##   [7] XVector_0.25.0           BiocNeighbors_1.3.1     
##   [9] proxy_0.4-23             bit64_0.9-7             
##  [11] ranger_0.11.2            codetools_0.2-16        
##  [13] robustbase_0.93-5        BiocManager_1.30.4      
##  [15] compiler_3.6.0           httr_1.4.0              
##  [17] dqrng_0.2.1              assertthat_0.2.1        
##  [19] Matrix_1.2-17            lazyeval_0.2.2          
##  [21] BiocSingular_1.1.1       htmltools_0.3.6         
##  [23] tools_3.6.0              rsvd_1.0.0              
##  [25] igraph_1.2.4.1           gtable_0.3.0            
##  [27] glue_1.3.1               GenomeInfoDbData_1.2.1  
##  [29] dplyr_0.8.1              ggthemes_4.2.0          
##  [31] rappdirs_0.3.1           Rcpp_1.0.1              
##  [33] carData_3.0-2            cellranger_1.1.0        
##  [35] DelayedMatrixStats_1.7.0 lmtest_0.9-37           
##  [37] laeken_0.5.0             xfun_0.7                
##  [39] stringr_1.4.0            ps_1.3.0                
##  [41] openxlsx_4.1.0           irlba_2.3.3             
##  [43] statmod_1.4.30           edgeR_3.27.3            
##  [45] DEoptimR_1.0-8           zoo_1.8-5               
##  [47] zlibbioc_1.31.0          MASS_7.3-51.4           
##  [49] scales_1.0.0             VIM_4.8.0               
##  [51] hms_0.4.2                yaml_2.2.0              
##  [53] curl_3.3                 memoise_1.1.0           
##  [55] gridExtra_2.3            stringi_1.4.3           
##  [57] RSQLite_2.1.1            highr_0.8               
##  [59] simpleSingleCell_1.9.3   e1071_1.7-1             
##  [61] destiny_2.15.0           TTR_0.23-4              
##  [63] boot_1.3-22              zip_2.0.2               
##  [65] rlang_0.3.4              pkgconfig_2.0.2         
##  [67] bitops_1.0-6             evaluate_0.13           
##  [69] lattice_0.20-38          purrr_0.3.2             
##  [71] labeling_0.3             cowplot_0.9.4           
##  [73] bit_1.1-14               processx_3.3.1          
##  [75] tidyselect_0.2.5         plyr_1.8.4              
##  [77] magrittr_1.5             bookdown_0.10           
##  [79] R6_2.4.0                 DBI_1.0.0               
##  [81] pillar_1.4.0             haven_2.1.0             
##  [83] foreign_0.8-71           withr_2.1.2             
##  [85] xts_0.11-2               scatterplot3d_0.3-41    
##  [87] abind_1.4-5              RCurl_1.95-4.12         
##  [89] sp_1.3-1                 nnet_7.3-12             
##  [91] tibble_2.1.1             crayon_1.3.4            
##  [93] car_3.0-2                rmarkdown_1.12          
##  [95] viridis_0.5.1            locfit_1.5-9.1          
##  [97] grid_3.6.0               data.table_1.12.2       
##  [99] blob_1.1.1               callr_3.2.0             
## [101] forcats_0.4.0            vcd_1.4-4               
## [103] digest_0.6.19            munsell_0.5.0           
## [105] beeswarm_0.2.3           viridisLite_0.3.0       
## [107] smoother_1.1             vipor_0.4.5
```

# References

