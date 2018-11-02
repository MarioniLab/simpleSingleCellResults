---
title: Detecting differental expression from single-cell RNA-seq data
author: 
- name: Aaron T. L. Lun
  affiliation: &CRUK Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
date: "2018-11-02"
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

Previous workflows ([here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/work-1-reads.html#detecting-marker-genes-between-clusters) and [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/work-5-mnn.html#using-the-corrected-values-in-downstream-analyses)) used the `block=` argument in `findMarkers()` to account for uninteresting factors of variation.
This will instruct `findMarkers()` to perform pairwise $t$-tests between clusters using only cells on the same level of the blocking factor. 
It will then combine $p$-values from different plates using Stouffer's Z method to obtain a single $p$-value per gene.


```r
library(SingleCellExperiment)
sce.pancreas <- readRDS("pancreas_data.rds") 

# Same code as in pancreas MNN correction workflow.
library(scran)
m.out <- findMarkers(sce.pancreas, sce.pancreas$Cluster, 
    block=sce.pancreas$Batch, direction="up") 
demo <- m.out[["1"]] 
as.data.frame(demo[demo$Top <= 5,1:3])
```

```
##           Top      p.value          FDR
## KCNQ1OT1    1 3.670707e-54 5.417229e-50
## PGM5P2      1 6.649785e-49 3.271251e-45
## TMEM212     1 4.539580e-35 1.116585e-31
## TLCD2       1 3.281998e-27 4.036310e-24
## SMG1        1 9.546806e-25 8.287751e-22
## TSIX        1 3.295468e-17 5.286360e-15
## UGDH-AS1    2 1.668432e-52 1.231136e-48
## LOC643406   2 1.427030e-32 3.008587e-29
## ODF2L       2 1.324449e-24 1.085901e-21
## NLRP12      2 2.134796e-23 1.369797e-20
## ARMC9       2 5.660849e-16 6.847771e-14
## MAB21L3     3 6.950271e-42 2.564302e-38
## CCL5        3 7.914665e-27 7.786975e-24
## FBXL20      3 2.919750e-22 1.267343e-19
## TUBA3FP     3 1.927714e-14 1.734707e-12
## TFDP2       4 3.875190e-28 5.719006e-25
## LRRC57      4 7.788701e-22 3.192935e-19
## GPR155      4 2.261856e-14 1.998830e-12
## FBLIM1      5 4.101208e-37 1.210513e-33
## PNPT1       5 2.596178e-20 7.819264e-18
## SLC16A7     5 4.626478e-14 3.969626e-12
```

Intra-batch comparisons with `block=` are robust to difference in the log-fold changes or variance between batches.
However, we need to assume that each pair of clusters is present in at least one batch.
In scenarios where cells from two clusters never co-occur in the same batch, the comparison will be impossible and `NA`s will be reported in the output.

## Using the `design=` argument

Another approach is to define a design matrix containing the batch of origin as the sole factor.
`findMarkers()` will then fit a linear model to the log-expression values, similar to the use of *[limma](https://bioconductor.org/packages/3.9/limma)* for bulk RNA sequencing data [@ritchie2015limma].
This handles situations where multiple batches contain unique clusters, as comparisons can be implicitly performed via shared cell types in each batch.
There is also a slight increase in power when information is shared across clusters for variance estimation.


```r
# Setting up the design matrix (we remove intercept for full rank
# in the final design matrix with the cluster-specific terms).
design <- model.matrix(~sce.pancreas$Batch)
design <- design[,-1,drop=FALSE]

m.alt <- findMarkers(sce.pancreas, sce.pancreas$Cluster, 
    design=design, direction="up")
demo <- m.alt[["1"]]
as.data.frame(demo[demo$Top <= 5,1:3])
```

```
##              Top       p.value           FDR
## CCL5           1  0.000000e+00  0.000000e+00
## PGM5P2         1  0.000000e+00 1.482197e-323
## LPAL2          2  0.000000e+00  0.000000e+00
## PCDH11Y        2  0.000000e+00  0.000000e+00
## TFDP2          2 2.721679e-314 8.033307e-311
## VSTM4          3 2.755068e-298 3.696299e-295
## TMEM212        4 1.341884e-308 3.300589e-305
## ARHGEF26-AS1   4 7.491362e-302 1.228417e-298
## LRRC57         4 4.976489e-264 2.225546e-261
## IBA57          5 3.062391e-308 6.456396e-305
## FBXL20         5 6.455013e-254 2.323490e-251
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
We demonstrate the use of `overlapExprs()` on the 416B data set from the [previous workflow](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/work-1-reads.html),
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
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3b-de_files/figure-html/viol-de-wilcox-1.png" alt="Distribution of log-normalized expression values for the top 10 DE genes involving cluster 1 with the Wilcoxon rank sum test, stratified by cluster assignment and coloured by the plate of origin for each cell." width="100%" />
<p class="caption">(\#fig:viol-de-wilcox)Distribution of log-normalized expression values for the top 10 DE genes involving cluster 1 with the Wilcoxon rank sum test, stratified by cluster assignment and coloured by the plate of origin for each cell.</p>
</div>

Wilcoxon tests provide a stronger guarantee of cluster separation than the $t$-tests in `findMarkers()` as the latter can have large effect sizes driven by a minority of cells.
This promotes the identification of good marker genes that discriminate between clusters.
The downside is that they are slower, the effect size is more difficult to interpret and the test result is not entirely robust to differences in scaling biases across clusters.

# Caveats with interpreting DE p-values

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
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3b-de_files/figure-html/pval-dist-1.png" alt="Distribution of p-values from a DE analysis between two clusters in a simulation with no true subpopulation structure." width="100%" />
<p class="caption">(\#fig:pval-dist)Distribution of p-values from a DE analysis between two clusters in a simulation with no true subpopulation structure.</p>
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
This yields a single pseudo-bulk count profile per combination, to which standard methods like *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* or *[limma](https://bioconductor.org/packages/3.9/limma)* can be applied.

We demonstrate this procedure on the [416B data set](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/work-1-reads.html) again.
We first create a factor containing all combinations of the per-cell factors of interest.
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

We average the count profiles for all cells in each level of `combos`^[After dividing by size factors in `sce.416b`, to avoid unnecessary variability due to technical scaling biases.].
This yields a set of pseudo-bulk samples that are more amenable to standard DE analysis as the counts are higher and per-observation variance is lower.
It also ensures that the variance is modelled across samples rather than across cells.
Each sample is represented no more than once for each condition in the `averaged` matrix, avoiding problems from unmodelled correlations between samples. 


```r
# Computing (non-transformed) normalized expression values. 
sce.416b <- normalize(sce.416b, return_log=FALSE)

# Averaging within each level of 'combos'.
library(DelayedMatrixStats) 
averaged <- colsum(normcounts(sce.416b), group=combos)
num <- num[match(colnames(averaged), names(num))]
averaged <- sweep(averaged, 2, num, "/")
head(averaged)
```

```
##                     1.20160113 1.20160325 2.20160113 2.20160325 3.20160113
## ENSMUSG00000103377 0.000000000 0.00000000 0.00000000 0.00000000          0
## ENSMUSG00000103147 0.000000000 0.25874285 0.00000000 0.00000000          0
## ENSMUSG00000103161 0.000000000 0.00000000 0.00000000 0.04663610          0
## ENSMUSG00000102331 0.000000000 0.00000000 0.16038423 0.00000000          0
## ENSMUSG00000102948 0.006832898 0.04588272 0.05571614 0.07970595          0
## Rp1                3.354680762 0.84325911 0.00000000 0.06534831          0
##                    3.20160325 4.20160113 4.20160325 5.20160113 5.20160325
## ENSMUSG00000103377 0.00000000          0          0  0.3567784  0.0000000
## ENSMUSG00000103147 0.00000000          0          0  0.0000000  0.0000000
## ENSMUSG00000103161 0.00000000          0          0  0.0000000  0.0000000
## ENSMUSG00000102331 0.14882442          0          0  0.0000000  0.0000000
## ENSMUSG00000102948 0.07015579          0          0  0.0000000  1.1876183
## Rp1                0.00000000          0          0  0.0000000  0.2775207
```

Finally, we perform a standard *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* analysis using the quasi-likelihood framework [@chen2016from].
We ignore spike-in transcripts and low-abundance genes;
set up the design matrix so that the plate of origin is an additive effect;
weight each pseudo-bulk count by the number of cells used to compute it; 
and test for DE between the first two clusters.


```r
library(edgeR)
y <- DGEList(averaged)
y <- y[aveLogCPM(y) > 1 & !isSpike(sce.416b),]
y <- calcNormFactors(y)

# Adding a matrix of sample weights.
w <- as.numeric(num)
y$weights <- makeCompressedMatrix(w, dim(y), byrow=TRUE)

# Building the design matrix.
ave.factors <- strsplit(colnames(y), split="\\.")
ave.clust <- unlist(lapply(ave.factors, "[[", i=1))
ave.plate <- unlist(lapply(ave.factors, "[[", i=2))
design <- model.matrix(~0 + ave.clust + ave.plate)

# Running through the remaining edgeR functions.
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, contrast=c(1, -1, 0, 0, 0, 0))
(top <- topTags(res))
```

```
## Coefficient:  1*ave.clust1 -1*ave.clust2 
##           logFC   logCPM         F       PValue          FDR
## H2afx -4.786267 7.301251 2075.9990 1.268793e-09 1.802320e-05
## Cenpf -4.019964 7.348619 1338.0335 5.519868e-09 3.920486e-05
## Ube2c -6.597884 8.051286  938.3938 2.751505e-08 1.302838e-04
## Plk4  -6.271306 5.292440  656.1603 5.954153e-08 1.358478e-04
## Lst1   3.311667 7.460206  612.8069 7.475172e-08 1.358478e-04
## Spc24 -5.172029 5.726840  594.7888 8.255381e-08 1.358478e-04
## Srm   -3.179802 7.995532  650.8267 8.981760e-08 1.358478e-04
## Cip2a -4.582759 5.752809  573.9298 9.296009e-08 1.358478e-04
## Mki67 -3.053341 8.010079  628.2231 1.006734e-07 1.358478e-04
## Cks1b -3.992811 7.592255  617.4474 1.064543e-07 1.358478e-04
```

```r
summary(decideTests(res))
```

```
##        1*ave.clust1 -1*ave.clust2
## Down                         1580
## NotSig                      11634
## Up                            991
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
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3b-de_files/figure-html/violde-1.png" alt="Distribution of log-normalized expression values for the top 10 DE genes between clusters 1 and 2 from the summation strategy, stratified by cluster assignment and coloured by the plate of origin for each cell." width="100%" />
<p class="caption">(\#fig:violde)Distribution of log-normalized expression values for the top 10 DE genes between clusters 1 and 2 from the summation strategy, stratified by cluster assignment and coloured by the plate of origin for each cell.</p>
</div>

**Comments from Aaron:**

- Note that the data dredging problem mentioned previously is an orthogonal issue to variance modelling.
This will not be resolved by performing summation on data with empirically defined clusters.
Of course, the summation approach can also be easily used with _a priori_ defined groups for which data dredging is not a problem.

# References
