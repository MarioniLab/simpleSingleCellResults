---
title: Detecting differental expression from single-cell RNA-seq data
author: 
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
date: "2019-05-20"
vignette: >
  %\VignetteIndexEntry{10. Detecting differential expression}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---
    




# Overview

Here, we describe some of the more theoretical aspects of detecting differential expression (DE) from single-cell RNA sequencing (scRNA-seq) data.
This includes the basis of blocking on uninteresting factors of variation in `findMarkers()`;
     the use of Wilcoxon rank sum tests in `overlapExprs()`;
     incorporating other DE analysis results with `combineMarkers()`;
     and some caveats on the interpretation of DE $p$-values in scRNA-seq contexts.

# Blocking on uninteresting factors of variation

## Using the `block=` argument

Previous workflows ([here](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/reads.html#detecting-marker-genes-between-clusters) and [here](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/batch.html#using-the-corrected-values-in-downstream-analyses)) used the `block=` argument in `findMarkers()` to account for uninteresting factors of variation.
This will instruct `findMarkers()` to perform pairwise $t$-tests between clusters using only cells on the same level of the blocking factor. 
It will then combine $p$-values from different plates using Stouffer's Z method to obtain a single $p$-value per gene.


```r
library(SingleCellExperiment)
sce.pancreas <- readRDS("pancreas_data.rds") 

# Same code as in pancreas MNN correction workflow.
library(scran)
m.out <- findMarkers(sce.pancreas, sce.pancreas$Cluster, 
    block=sce.pancreas$batch, direction="up", assay.type="original") 
demo <- m.out[["1"]] 
as.data.frame(demo[demo$Top <= 5,1:3])
```

```
##          Top       p.value           FDR
## TMSB4X     1  0.000000e+00  0.000000e+00
## ANXA2      1 2.856151e-306 1.047279e-302
## KRT8       1 4.133937e-288 1.212649e-284
## TACSTD2    1 1.080810e-247 1.981529e-244
## CFTR       1 2.944056e-233 3.321574e-230
## SLC4A4     1 4.683639e-232 4.906781e-229
## KIAA1522   1 2.452923e-215 1.998724e-212
## IGFBP7     1 7.003013e-183 2.567830e-180
## VAMP8      1 2.517592e-173 6.838059e-171
## COL18A1    1 1.146431e-157 2.128444e-155
## CD24       2  0.000000e+00  0.000000e+00
## ANXA4      2 4.054275e-314 1.982135e-310
## TM4SF1     2 8.303170e-274 1.739751e-270
## BACE2      2 1.170391e-234 1.560557e-231
## SERPING1   2 2.126506e-220 2.079297e-217
## EPCAM      2 2.598941e-201 1.411803e-198
## DDR1       2 5.138729e-173 1.345888e-170
## ONECUT2    2 2.471746e-160 4.899067e-158
## ANXA3      2 3.836901e-139 5.024627e-137
## LSR        2 8.778561e-111 6.405729e-109
## S100A11    3 4.494620e-274 1.098710e-270
## CLDN7      3 1.302434e-207 7.347229e-205
## CDH1       3 9.071719e-182 3.167974e-179
## CTSH       3 1.513027e-174 4.351287e-172
## CLDN4      3 7.633508e-142 1.056233e-139
## CIB1       3  3.900792e-73  1.091850e-71
## LITAF      4 1.583840e-214 1.222641e-211
## KRT7       4 4.426852e-209 2.705360e-206
## PTPRF      4 8.732258e-180 2.784262e-177
## KRT19      4 2.534131e-172 6.520719e-170
## MYRF       4 5.739823e-170 1.426881e-167
## CLDN10     4 6.173173e-166 1.331499e-163
## RAB25      4 2.680846e-143 3.931997e-141
## UQCR10     4  4.585739e-71  1.192536e-69
## PMEPA1     5 2.501078e-214 1.834165e-211
## KRT18      5 4.318874e-174 1.195187e-171
## NFIB       5 4.083758e-169 9.819095e-167
## PROM1      5 7.321851e-113 5.564228e-111
## LMAN2      5  1.467657e-70  3.737175e-69
```

Intra-batch comparisons with `block=` are robust to difference in the log-fold changes or variance between batches.
However, we need to assume that each pair of clusters is present in at least one batch.
In scenarios where cells from two clusters never co-occur in the same batch, the comparison will be impossible and `NA`s will be reported in the output.

## Using the `design=` argument

Another approach is to define a design matrix containing the batch of origin as the sole factor.
`findMarkers()` will then fit a linear model to the log-expression values, similar to the use of *[limma](https://bioconductor.org/packages/3.10/limma)* for bulk RNA sequencing data [@ritchie2015limma].
This handles situations where multiple batches contain unique clusters, as comparisons can be implicitly performed via shared cell types in each batch.
There is also a slight increase in power when information is shared across clusters for variance estimation.


```r
# Setting up the design matrix (we remove intercept for full rank
# in the final design matrix with the cluster-specific terms).
design <- model.matrix(~sce.pancreas$batch)
design <- design[,-1,drop=FALSE]

m.alt <- findMarkers(sce.pancreas, sce.pancreas$Cluster, 
    design=design, direction="up", assay.type="original")
demo <- m.alt[["1"]]
as.data.frame(demo[demo$Top <= 5,1:3])
```

```
##          Top       p.value           FDR
## CFTR       1  0.000000e+00  0.000000e+00
## ANXA4      1  0.000000e+00  0.000000e+00
## SLC4A4     2  0.000000e+00  0.000000e+00
## SERPING1   2  0.000000e+00  0.000000e+00
## CD24       2  0.000000e+00  0.000000e+00
## KRT8       2  0.000000e+00  0.000000e+00
## B2M        2 6.785276e-166 5.156458e-164
## SPP1       3  0.000000e+00  0.000000e+00
## KRT19      3  0.000000e+00  0.000000e+00
## ANXA2      3  0.000000e+00  0.000000e+00
## RPL15      3 2.004173e-153 1.348404e-151
## TACSTD2    4  0.000000e+00  0.000000e+00
## EPCAM      4 1.572563e-267 2.957024e-265
## TMBIM6     4 1.740368e-149 1.105021e-147
## PMEPA1     5  0.000000e+00  0.000000e+00
## ALDH1A3    5  0.000000e+00  0.000000e+00
## CD59       5  0.000000e+00  0.000000e+00
```

The use of a linear model makes a few some strong assumptions, necessitating some caution when interpreting the results.
The batch effect across cell types is assumed to be homogeneous.
If this is not the case, the variance will be inflated and the log-fold change estimates will be distorted.
Variances are also assumed to be equal across groups, which is not true in general.
In particular, the presence of clusters in which a gene is silent will shrink the residual variance towards zero, 
preventing the model from penalizing genes with high variance in other clusters.
Thus, we generally recommend the use of `block=` where possible.

# Using the Wilcoxon rank sum test

The `overlapExprs()` function uses the Wilcoxon rank sum test to detect uneven mixing of the distributions of expression values between clusters.
The effect size is reported as the probability of randomly sampling one observation in one cluster that is greater than a random observation in another cluster.
This prioritizes genes where there is clear separation between the distributions of expression values of different clusters.
We demonstrate the use of `overlapExprs()` on the 416B data set from the [previous workflow](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/reads.html),
detecting DE genes between clusters while blocking on the plate of origin.


```r
sce.416b <- readRDS("416B_data.rds")
o.out <- overlapExprs(sce.416b, group=sce.416b$cluster, block=sce.416b$Plate)
head(o.out[["1"]]) # top DEGs for cluster 1 against the others.
```

```
## DataFrame with 6 rows and 7 columns
##                          Top              p.value                  FDR
##                    <integer>            <numeric>            <numeric>
## ENSMUSG00000107771         1 6.24804570972308e-22 1.48597271114344e-17
## Hbb-b1                     1 8.26717542905217e-21 6.74596321819895e-17
## ENSMUSG00000074026         1 4.93542635295853e-20 1.30421383280458e-16
## Oip5                       1 1.76459269663039e-15 1.83263354165767e-13
## ENSMUSG00000084220         2 5.70788175385249e-20 1.35750551751874e-16
## Eef1akmt1                  2 1.89928885228432e-17 5.23470093472323e-15
##                             overlap.2            overlap.3
##                             <numeric>            <numeric>
## ENSMUSG00000107771 0.0436177036561899    0.306451612903226
## Hbb-b1             0.0301475304682489   0.0930875576036866
## ENSMUSG00000074026 0.0814624759461193 0.000921658986175133
## Oip5                0.235086593970494    0.223963133640553
## ENSMUSG00000084220 0.0490699166132136  0.00368663594470042
## Eef1akmt1          0.0647851186658114   0.0829493087557603
##                             overlap.4          overlap.5
##                             <numeric>          <numeric>
## ENSMUSG00000107771  0.443514644351464   0.11605415860735
## Hbb-b1             0.0993723849372385  0.470019342359768
## ENSMUSG00000074026  0.178870292887029  0.131528046421663
## Oip5                0.381799163179916 0.0783365570599613
## ENSMUSG00000084220  0.326359832635983 0.0947775628626693
## Eef1akmt1          0.0899581589958159  0.421663442940039
```

Effect sizes close to zero indicate that the gene is downregulated, while effect sizes close to unity correspond to upregulation.
The top DE genes all exhibit strong separation between cluster 1 and the others (Figure \@ref(fig:viol-de-wilcox)).


```r
library(scater)
plotExpression(sce.416b, x="cluster", colour_by="Plate",
    features=head(rownames(o.out[[1]])))
```

<div class="figure">
<img src="de_files/figure-html/viol-de-wilcox-1.png" alt="Distribution of log-normalized expression values for the top 10 DE genes involving cluster 1 with the Wilcoxon rank sum test, stratified by cluster assignment and coloured by the plate of origin for each cell." width="100%" />
<p class="caption">(\#fig:viol-de-wilcox)Distribution of log-normalized expression values for the top 10 DE genes involving cluster 1 with the Wilcoxon rank sum test, stratified by cluster assignment and coloured by the plate of origin for each cell.</p>
</div>

Wilcoxon tests provide a stronger guarantee of cluster separation than the $t$-tests in `findMarkers()` as the latter can have large effect sizes driven by a minority of cells.
This promotes the identification of good marker genes that discriminate between clusters.
The downside is that they are slower, the effect size is more difficult to interpret and the test result is not entirely robust to differences in scaling biases across clusters.

# Using other DE analysis results

It is possible to perform marker gene detection based on results from other DE analysis methods.
For example, consider the `voom` approach in the *[limma](https://bioconductor.org/packages/3.10/limma)* package [@law2014voom].


```r
library(limma)
design <- model.matrix(~0 + cluster + Plate, data=colData(sce.416b))
colnames(design)
```

```
## [1] "cluster1"      "cluster2"      "cluster3"      "cluster4"     
## [5] "cluster5"      "Plate20160325"
```

```r
keep <- calcAverage(sce.416b) > 1 # filter to remove very low-abundance genes.
summary(keep)
```

```
##    Mode   FALSE    TRUE 
## logical   10630   13233
```

```r
y <- convertTo(sce.416b, subset.row=keep)
v <- voom(y, design)
fit <- lmFit(v, design)
```

We perform pairwise moderated $t$-tests between clusters while blocking on the plate of origin.
Here, we use the TREAT strategy [@mccarthy2009treat] to test for log-fold changes that are significantly greater than 0.5.


```r
clust.terms <- head(colnames(design), length(unique(sce.416b$cluster)))
all.results <- all.pairs <- list()
counter <- 1L

for (x in seq_along(clust.terms)) {
    for (y in seq_len(x-1L)) {
        con <- integer(ncol(design))
        con[x] <- 1
        con[y] <- -1
        fit2 <- contrasts.fit(fit, con)
        fit2 <- treat(fit2, robust=TRUE, lfc=0.5)

        res <- topTreat(fit2, n=Inf, sort.by="none")
        all.results[[counter]] <- res
        all.pairs[[counter]] <- c(clust.terms[x], clust.terms[y])
        counter <- counter+1L

        # Also filling the reverse comparison.
        res$logFC <- -res$logFC
        all.results[[counter]] <- res
        all.pairs[[counter]] <- c(clust.terms[y], clust.terms[x])
        counter <- counter+1L
    }
}
```

The results of this comparison are consolidated into a single marker list for each cluster with the `combineMarkers()` function.
This yields an ordering of genes that can be interpreted in the same manner as discussed [previously](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/reads.html#detecting-marker-genes-between-clusters) for `findMarkers()` output.


```r
all.pairs <- do.call(rbind, all.pairs)
combined <- combineMarkers(all.results, all.pairs, pval.field="P.Value")
as.data.frame(head(combined[["cluster1"]][,1:3]))
```

```
##        Top       p.value           FDR
## CYTB     1 9.466154e-162 1.246976e-157
## Pimreg   1  4.100252e-75  2.700631e-71
## Myh11    1  6.196199e-69  2.720751e-65
## Kif11    1  1.134306e-51  7.115337e-49
## Rrm2     2  1.046396e-64  2.756836e-61
## ND1      2  1.301918e-44  3.811149e-42
```

# Caveats with interpreting DE $p$-values

## Data dredging from clustering

It is worth noting that all of our DE strategies for detecting marker genes between clusters are statistically flawed to some extent.
The DE analysis is performed on the same data used to obtain the clusters, which represents "data dredging" (also known as fishing or data snooping).
The hypothesis of interest - that are there differences between clusters? - is formulated from the data, 
so we are more likely to get a positive result when we re-use the data set to test that hypothesis.

The practical effect of data dredging is best illustrated with a simple simulation.
We simulate i.i.d. normal values, perform k-means clustering and test for DE between clusters of cells with `findMarkers()`.
The resulting distribution of $p$-values is heavily skewed towards low values (Figure \@ref(fig:pval-dist)).
Thus, we can detect "significant" differences between clusters even in the absence of any real substructure in the data.
This effect arises from the fact that clustering, by definition, yields groups of cells that differ in their coordinates in expression space. 
Testing for DE genes between clusters will inevitably yield some significant results as that is how the clusters were defined in the first place.


```r
set.seed(0)
y <- matrix(rnorm(100000), ncol=200)
clusters <- kmeans(t(y), centers=2)$cluster
out <- findMarkers(y, clusters)
hist(out[[1]]$p.value, col="grey80", xlab="p-value")
```

<div class="figure">
<img src="de_files/figure-html/pval-dist-1.png" alt="Distribution of $p$-values from a DE analysis between two clusters in a simulation with no true subpopulation structure." width="100%" />
<p class="caption">(\#fig:pval-dist)Distribution of $p$-values from a DE analysis between two clusters in a simulation with no true subpopulation structure.</p>
</div>

By and large, this effect does not cause problems for marker gene detection as the DE statistics from `findMarkers()` and counterparts are primarily used for ranking.
It does become an issue when the $p$-values are used to define "significant differences" between clusters with respect to an error rate threshold.
Meaningful interpretation of error rates require consideration of the long-run behaviour, i.e., the rate of incorrect rejections if the experiment were repeated many times.
The concept of statistical significance for differences between clusters is not applicable if clusters are not stably reproducible across (hypothetical) replicate experiments.

To overcome this conceptual hurdle, we need to annotate our clusters based on a few marker genes.
This allows us to use the annotated clusters as proxies for the true (and presumably reproducible) biological subpopulations.
We might then be tempted to interpret the significant genes as being DE between subpopulations.
However, this would result in loss of error control when the clusters are not stable, due to overfitting of the cluster definitions for true null genes.
This effect is exacerbated as the clusters become more unstable, e.g., due to poor separation between the underlying populations.

## Considering sample effects

The naive application of DE analysis methods will treat counts from the same cluster of cells as replicate observations.
This is not the most relevant level of replication when cells are derived from the same biological sample (i.e., cell culture, animal or patient).
DE analyses that treat cells as replicates fail to properly model the sample-to-sample variability [@lun2017overcoming].
The latter is arguably the more important level of replication as different samples will necessarily be generated if the experiment is to be replicated.
Indeed, the use of cells as replicates only masks the fact that the sample size is actually one in an experiment involving a single biological sample.

The "most correct" strategy for accommodating two levels of replication is to use a (generalized) linear mixed model.
However, these are difficult to implement from both a theoretical and practical perspective - see [here](https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html) for an in-depth discussion.
A faster approach is to use a summation strategy [@lun2017overcoming], where all cells in each combination of sample and condition (or cluster) are summed together.
This yields a single pseudo-bulk count profile per combination, to which standard methods like *[edgeR](https://bioconductor.org/packages/3.10/edgeR)* or *[limma](https://bioconductor.org/packages/3.10/limma)* can be applied.

# Using pseudo-bulk counts 

We demonstrate this procedure on the [416B data set](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/reads.html) again.
We create a factor containing all combinations of the per-cell factors of interest.
Here, the factors of interest are the assigned cluster and the plate of origin for each cell^[Cluster is nested in oncogene induction status so that latter will not be used here.].
All cells from one plate with the same oncogene induction status were obtained from the same biological sample.


```r
combos <- with(colData(sce.416b), paste(cluster, Plate, sep="."))
(num <- table(combos))
```

```
## combos
## 1.20160113 1.20160325 2.20160113 2.20160325 3.20160113 3.20160325 4.20160113 
##         41         39         19         20         16         11         10 
## 4.20160325 5.20160113 5.20160325 
##         14          5          8
```

We sum the count profiles for all cells in each level of `combos`.
This yields a set of pseudo-bulk samples that are more amenable to standard DE analysis as the counts are higher and per-observation variance is lower.
It also ensures that the variance is modelled across samples rather than across cells.
Each sample is represented no more than once for each condition in the `summed` matrix, avoiding problems from unmodelled correlations between samples. 


```r
library(scater)
summed <- sumCountsAcrossCells(counts(sce.416b), combos)
head(summed)
```

```
##                    1.20160113 1.20160325 2.20160113 2.20160325 3.20160113
## ENSMUSG00000103377          0          0          0          0          0
## ENSMUSG00000103147          0          7          0          0          0
## ENSMUSG00000103161          0          0          0          2          0
## ENSMUSG00000102331          0          0          4          0          0
## ENSMUSG00000102948          1          1          1          3          0
## Rp1                        97         48          0          1          0
##                    3.20160325 4.20160113 4.20160325 5.20160113 5.20160325
## ENSMUSG00000103377          0          0          0          2          0
## ENSMUSG00000103147          0          0          0          0          0
## ENSMUSG00000103161          0          0          0          0          0
## ENSMUSG00000102331          2          0          0          0          0
## ENSMUSG00000102948          1          0          0          0          9
## Rp1                         0          0          0          0          2
```

We perform a standard *[edgeR](https://bioconductor.org/packages/3.10/edgeR)* analysis using the quasi-likelihood (QL) framework [@chen2016from] to identify differential expression between clusters.
We ignore spike-in transcripts and low-abundance genes, and we compute normalization factors^[Not the same as size factors!] using the `calcNormFactors()` function.


```r
library(edgeR)
y <- DGEList(summed)
y <- y[aveLogCPM(y) > 1 & !isSpike(sce.416b),]
y <- calcNormFactors(y)
```

We set up the design matrix so that the plate of origin is an additive effect.


```r
sum.terms <- strsplit(colnames(y), split="\\.")
sum.clust <- unlist(lapply(sum.terms, "[[", i=1))
sum.plate <- unlist(lapply(sum.terms, "[[", i=2))
design <- model.matrix(~0 + sum.clust + sum.plate)
```

We estimate the negative binomial and QL dispersions using `estimateDisp()` and `glmQLFit()`, respectively.
Larger dispersions represent greater variability between replicate plates, not cells.


```r
y <- estimateDisp(y, design)
summary(y$trended.dispersion)    
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.02974 0.11443 0.25894 0.38510 0.60977 1.06925
```

```r
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)    
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   2.737   2.737   2.737   2.788   2.859   2.920
```

We test for DE between the first two clusters using `glmQLFTest()`.
This yields a number of potentially interesting DE genes involved in cell division,
consistent with the differences in cell cycle activity between clusters.


```r
res <- glmQLFTest(fit, contrast=c(1, -1, 0, 0, 0, 0))
summary(decideTests(res))
```

```
##        1*sum.clust1 -1*sum.clust2
## Down                          707
## NotSig                      10439
## Up                            405
```

```r
(top <- topTags(res))
```

```
## Coefficient:  1*sum.clust1 -1*sum.clust2 
##            logFC   logCPM        F       PValue          FDR
## H2afx  -4.705118 7.660453 759.0851 2.483667e-08 0.0002868884
## Cenpf  -3.908109 7.787069 388.5249 2.930229e-07 0.0016923539
## Pmf1   -3.489725 6.825686 264.7124 9.031752e-07 0.0020504630
## Mki67  -3.071246 8.525429 263.9944 1.182572e-06 0.0020504630
## Cdk1   -3.878273 8.772554 261.4039 1.221613e-06 0.0020504630
## Dut    -3.222792 7.980844 259.1590 1.256806e-06 0.0020504630
## Ube2c  -6.562471 8.544501 243.4641 1.543614e-06 0.0020504630
## Mad2l1 -3.334600 7.228002 224.9613 1.773234e-06 0.0020504630
## Ctsg   -9.810366 4.476225 210.1129 1.968233e-06 0.0020504630
## Srm    -3.371245 8.179901 225.0141 1.999676e-06 0.0020504630
```

The DE analysis on pseudo-bulk counts does not explicitly penalize DE genes that are highly variable across cells within each sample.
One could argue that this is undesirable as DE genes with low expression variance across cells are more discriminative and should be more highly ranked.
In practice, this is less of an issue than might be expected (Figure \@ref(fig:violde)).
The effect size will naturally be smaller when cell-to-cell expression is highly variable, simply because the means of each group must be closer together when the distributions mix.
This means that the analysis still implicitly favours DE genes with low cell-to-cell variance.
Summation also avoids harshly penalizing genes for high technical variability, which might otherwise reduce the ranking of good low-abundance markers. 


```r
sub.sce <- sce.416b[,sce.416b$cluster %in% c("1", "2")]
plotExpression(sub.sce, x="cluster", colour_by="Plate",
    features=rownames(top)[1:10])
```

<div class="figure">
<img src="de_files/figure-html/violde-1.png" alt="Distribution of log-normalized expression values for the top 10 DE genes between clusters 1 and 2 from the summation strategy, stratified by cluster assignment and coloured by the plate of origin for each cell." width="100%" />
<p class="caption">(\#fig:violde)Distribution of log-normalized expression values for the top 10 DE genes between clusters 1 and 2 from the summation strategy, stratified by cluster assignment and coloured by the plate of origin for each cell.</p>
</div>

**Comments from Aaron:**

- Note that the data dredging problem mentioned previously is an orthogonal issue to variance modelling.
This will not be resolved by performing summation on data with empirically defined clusters.
Of course, the summation approach can also be easily used with _a priori_ defined groups for which data dredging is not a problem.

# References
