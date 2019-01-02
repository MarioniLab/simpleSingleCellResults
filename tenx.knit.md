---
title: Analyzing single-cell RNA sequencing data from droplet-based protocols
author: 
- name: Aaron T. L. Lun
  affiliation: Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
date: "2018-12-26"
vignette: >
  %\VignetteIndexEntry{04. Droplet-based data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---



# Overview 

Droplet-based scRNA-seq protocols capture cells in droplets for massively multiplexed library prepation [@klein2015droplet;@macosko2015highly].
This greatly increases the throughput of scRNA-seq studies, allowing tens of thousands of individual cells to be profiled in a routine experiment.
However, it (again) involves some differences from the previous workflows to reflect some unique aspects of droplet-based data.

Here, we describe a brief analysis of the peripheral blood mononuclear cell (PBMC) dataset from 10X Genomics [@zheng2017massively].
The data are publicly available from the [10X Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k), 
from which we download the raw gene/barcode count matrices, i.e., before cell calling from the _CellRanger_ pipeline.


```r
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir="pbmc4k")
```

# Setting up the data

## Reading in a sparse matrix

We load in the raw count matrix using the `read10xCounts()` function from the *[DropletUtils](https://bioconductor.org/packages/3.9/DropletUtils)* package.
This will create a `SingleCellExperiment` object where each column corresponds to a cell barcode.


```r
library(DropletUtils)
fname <- "pbmc4k/raw_gene_bc_matrices/GRCh38"
sce <- read10xCounts(fname, col.names=TRUE)
sce
```

```
## class: SingleCellExperiment 
## dim: 33694 737280 
## metadata(0):
## assays(1): counts
## rownames(33694): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
##   ENSG00000268674
## rowData names(2): ID Symbol
## colnames(737280): AAACCTGAGAAACCAT-1 AAACCTGAGAAACCGC-1 ...
##   TTTGTCATCTTTAGTC-1 TTTGTCATCTTTCCTC-1
## colData names(2): Sample Barcode
## reducedDimNames(0):
## spikeNames(0):
```

Here, each count represents the number of unique molecular identifiers (UMIs) assigned to a gene for a cell barcode.
Note that the counts are loaded as a sparse matrix object - specifically, a `dgCMatrix` instance from the *[Matrix](https://CRAN.R-project.org/package=Matrix)* package.
This avoids allocating memory to hold zero counts, which is highly memory-efficient for low-coverage scRNA-seq data.


```r
class(counts(sce))
```

```
## [1] "dgCMatrix"
## attr(,"package")
## [1] "Matrix"
```

## Annotating the rows

We relabel the rows with the gene symbols for easier reading.
This is done using the `uniquifyFeatureNames()` function, which ensures uniqueness in the case of duplicated or missing symbols.


```r
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
head(rownames(sce))
```

```
## [1] "RP11-34P13.3"  "FAM138A"       "OR4F5"         "RP11-34P13.7" 
## [5] "RP11-34P13.8"  "RP11-34P13.14"
```

We also identify the chromosomal location for each gene.
The mitochondrial location is particularly useful for later quality control.


```r
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
    column="SEQNAME", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="MT")
```

```
##    Mode   FALSE    TRUE    NA's 
## logical   33537      13     144
```

# Calling cells from empty droplets

## Testing for deviations from ambient expression

An interesting aspect of droplet-based data is that we have no prior knowledge about which droplets (i.e., cell barcodes) actually contain cells, and which are empty.
Thus, we need to call cells from empty droplets based on the observed expression profiles.
This is not entirely straightforward as empty droplets can contain ambient (i.e., extracellular) RNA that can be captured and sequenced.
The distribution of total counts exhibits a sharp transition between barcodes with large and small total counts (Figure \@ref(fig:rankplot)),
probably corresponding to cell-containing and empty droplets respectively.


```r
bcrank <- barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
    xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
	col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
```

<div class="figure">
<img src="tenx_files/figure-html/rankplot-1.png" alt="Total UMI count for each barcode in the PBMC dataset, plotted against its rank (in decreasing order of total counts). The inferred locations of the inflection and knee points are also shown." width="100%" />
<p class="caption">(\#fig:rankplot)Total UMI count for each barcode in the PBMC dataset, plotted against its rank (in decreasing order of total counts). The inferred locations of the inflection and knee points are also shown.</p>
</div>

We use the `emptyDrops()` function to test whether the expression profile for each cell barcode is significantly different from the ambient RNA pool [@lun2018distinguishing].
Any significant deviation indicates that the barcode corresponds to a cell-containing droplet.
We call cells at a false discovery rate (FDR) of 0.1%, meaning that no more than 0.1% of our called barcodes should be empty droplets on average.


```r
set.seed(100)
e.out <- emptyDrops(counts(sce))
sum(e.out$FDR <= 0.001, na.rm=TRUE)
```

```
## [1] 4224
```

We then subset our `SingleCellExperiment` object to retain only the detected cells.


```r
# using which() to automatically remove NAs.
sce <- sce[,which(e.out$FDR <= 0.001)]
```

**Comments from Aaron:**

- `emptyDrops()` computes Monte Carlo _p_-values based on a Dirichlet-multinomial model of sampling molecules into droplets.
These _p_-values are stochastic so it is important to set the random seed to obtain reproducible results.
- The function assumes that cell barcodes with total UMI counts below a certain threshold (`lower=100` by default) correspond to empty droplets, 
and uses them to estimate the ambient expression profile.
By definition, these barcodes cannot be cell-containing droplets and are excluded from the hypothesis testing, hence the `NA`s in the output.
- Users wanting to use the cell calling algorithm from the _CellRanger_ pipeline can call `defaultDrops()` instead.
This tends to be quite conservative as it often discards genuine cells with low RNA content (and thus low total counts).
It also requires an estimate of the expected number of cells in the experiment.

## Examining cell-calling diagnostics

The number of Monte Carlo iterations (specified by the `niters` argument in `emptyDrops()`) determines the lower bound for the _p_values [@phipson2010permutation].
The `Limited` field in the output indicates whether or not the computed _p_-value for a particular barcode is bounded by the number of iterations.
If any non-significant barcodes are `TRUE` for `Limited`, we may need to increase the number of iterations.
A larger number of iterations will often result in a lower _p_-value for these barcodes, which may allow them to be detected after correcting for multiple testing.


```r
table(Sig=e.out$FDR <= 0.01, Limited=e.out$Limited)
```

```
##        Limited
## Sig     FALSE TRUE
##   FALSE   929    0
##   TRUE   1785 2575
```

As mentioned above, `emptyDrops()` assumes that barcodes with low total UMI counts are empty droplets.
Thus, the null hypothesis should be true for all of these barcodes. 
We can check whether the hypothesis test holds its size by examining the distribution of _p_-values for low-total barcodes.
Ideally, the distribution should be close to uniform.


```r
full.data <- read10xCounts(fname, col.names=TRUE)
set.seed(100)
limit <- 100   
all.out <- emptyDrops(counts(full.data), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
    xlab="P-value", main="", col="grey80") 
```

<div class="figure">
<img src="tenx_files/figure-html/ambientpvalhist-1.png" alt="Distribution of p-values for the assumed empty droplets." width="100%" />
<p class="caption">(\#fig:ambientpvalhist)Distribution of p-values for the assumed empty droplets.</p>
</div>

Large peaks near zero indicate that barcodes with total counts below `lower` are not all ambient in origin.
This can be resolved by decreasing `lower` further to exclude barcodes corresponding to droplets with very small cells.

# Quality control on the cells

The previous step only distinguishes cells from empty droplets, but makes no statement about the quality of the cells.
It is entirely possible for droplets to contain damaged or dying cells, which need to be removed prior to downstream analysis.
We compute some QC metrics using `calculateQCMetrics()` [@mccarthy2017scater] and examine their distributions in Figure \@ref(fig:qchist).


```r
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")))
par(mfrow=c(1,3))
hist(sce$log10_total_counts, breaks=20, col="grey80",
    xlab="Log-total UMI count")
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
    xlab="Log-total number of expressed features")
hist(sce$pct_counts_Mito, breaks=20, col="grey80",
	xlab="Proportion of reads in mitochondrial genes")
```

<div class="figure">
<img src="tenx_files/figure-html/qchist-1.png" alt="Histograms of QC metric distributions in the PBMC dataset." width="960" />
<p class="caption">(\#fig:qchist)Histograms of QC metric distributions in the PBMC dataset.</p>
</div>

Ideally, we would remove cells with low library sizes or total number of expressed features as described [previously](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/reads.html#quality-control-on-the-cells).
However, this would likely remove cell types with low RNA content, especially in a heterogeneous PBMC population with many different cell types.
Thus, we use a more relaxed strategy and only remove cells with large mitochondrial proportions, using it as a proxy for cell damage.
(Keep in mind that droplet-based datasets usually do not have spike-in RNA.)


```r
high.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")
sce <- sce[,!high.mito]
summary(high.mito)
```

```
##    Mode   FALSE    TRUE 
## logical    3917     307
```

**Comments from Aaron:**

- The above justification for using a more relaxed filter is largely retrospective.
In practice, we may not know _a priori_ the degree of population heterogeneity and whether it manifests in the QC metrics.
We recommend performing the analysis first with a stringent QC filter, and then relaxing it based on further diagnostics (see [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/qc.html#checking-for-discarded-cell-types) for an example).

# Examining gene expression

The average expression of each gene is much lower here compared to the previous datasets (Figure \@ref(fig:abhist)).
This is due to the reduced coverage per cell when thousands of cells are multiplexed together for sequencing.


```r
ave <- calcAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80")
```

<div class="figure">
<img src="tenx_files/figure-html/abhist-1.png" alt="Histogram of the log~10~-average counts for each gene in the PBMC dataset." width="100%" />
<p class="caption">(\#fig:abhist)Histogram of the log~10~-average counts for each gene in the PBMC dataset.</p>
</div>

The set of most highly expressed genes is dominated by ribosomal protein and mitochondrial genes (Figure \@ref(fig:highexpr)), as expected.


```r
plotHighestExprs(sce)
```

<div class="figure">
<img src="tenx_files/figure-html/highexpr-1.png" alt="Percentage of total counts assigned to the top 50 most highly-abundant features in the PBMC dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell." width="100%"  class="widefigure" />
<p class="caption">(\#fig:highexpr)Percentage of total counts assigned to the top 50 most highly-abundant features in the PBMC dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell.</p>
</div>

# Normalizing for cell-specific biases

We perform some pre-clustering as described [previously](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/umis.html#normalization-of-cell-specific-biases).
Recall that we need to set the seed when using `pc.approx=TRUE` due to the randomness of *[irlba](https://CRAN.R-project.org/package=irlba)*. 


```r
library(scran)
set.seed(1000)
clusters <- quickCluster(sce, use.ranks=FALSE, pc.approx=TRUE)
table(clusters)
```

```
## clusters
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
## 160 158 263 421 387 244 282 337 557 104 441 103 117 159 184
```

We apply the deconvolution method to compute size factors for all cells [@lun2016pooling].
Again, we use `min.mean=0.1` to account for the fact that UMI counts are lower.
The specification of `cluster=` also ensures that we do not pool cells that are very different.


```r
sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce))
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##  0.003323  0.708549  0.863859  1.000000  1.090025 11.550483
```

The size factors are well correlated against the library sizes (Figure \@ref(fig:sfplot)), indicating that capture efficiency and sequencing depth are the major biases.


```r
plot(sce$total_counts, sizeFactors(sce), log="xy")
```

<div class="figure">
<img src="tenx_files/figure-html/sfplot-1.png" alt="Size factors for all cells in the PBMC dataset, plotted against the library size." width="100%" />
<p class="caption">(\#fig:sfplot)Size factors for all cells in the PBMC dataset, plotted against the library size.</p>
</div>

Finally, we compute normalized log-expression values.
There is no need to call `computeSpikeFactors()` here, as there are no spike-in transcripts available.


```r
sce <- normalize(sce)
```

**Comments from Aaron:**

- Larger droplet-based datasets will often be generated in separate batches or runs.
In such cases, we can set `block=` in `quickCluster()` to cluster cells within each batch or run.
This reduces computational work considerably without compromising performance, provided that there are enough cells within each batch.
- Even in the absence of any known batch structure, we can improve speed by setting an arbitrary factor, e.g., using `block=cut(seq_len(ncol(sce)), 10)` to split the cells into ten "batches" of roughly equal size.
Recall that we are not interpreting the clusters themselves, so it is not a problem to have multiple redundant cluster labels.
Again, this assumes that each cluster is large enough to support deconvolution.
- On a similar note, both `quickCluster()` and `computeSumFactors()` can process blocks or clusters in parallel.
This is achieved using the *[BiocParallel](https://bioconductor.org/packages/3.9/BiocParallel)* framework, which provides support for a range of parallelization strategies.

# Modelling the mean-variance trend

The lack of spike-in transcripts complicates the modelling of the technical noise.
One option is to assume that most genes do not exhibit strong biological variation, and to fit a trend to the variances of endogenous genes
(see [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/var.html#when-spike-ins-are-unavailable) for details).
However, this assumption is generally unreasonable for a heterogeneous population.
Instead, we assume that the technical noise is Poisson and create a fitted trend on that basis using the `makeTechTrend()` function.


```r
new.trend <- makeTechTrend(x=sce)
```

We estimate the variances for all genes and compare the trend fits in Figure \@ref(fig:trendplot).
The Poisson-based trend serves as a lower bound for the variances of the endogenous genes.
This results in non-zero biological components for most genes, which is consistent with other UMI-based data sets 
(see the [corresponding analysis](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/umis.html#modelling-and-removing-technical_noise) of the @zeisel2015brain data set).


```r
fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
plot(fit$mean, fit$var, pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
curve(new.trend(x), col="red", add=TRUE)
```

<div class="figure">
<img src="tenx_files/figure-html/trendplot-1.png" alt="Variance of normalized log-expression values for each gene in the PBMC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances, while the red line represents the Poisson noise." width="100%" />
<p class="caption">(\#fig:trendplot)Variance of normalized log-expression values for each gene in the PBMC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances, while the red line represents the Poisson noise.</p>
</div>

We decompose the variance for each gene using the Poisson-based trend, and examine the genes with the highest biological components.


```r
fit0 <- fit
fit$trend <- new.trend
dec <- decomposeVar(fit=fit)
top.dec <- dec[order(dec$bio, decreasing=TRUE),] 
head(top.dec)
```

```
## DataFrame with 6 rows and 6 columns
##                     mean            total              bio              tech
##                <numeric>        <numeric>        <numeric>         <numeric>
## LYZ     1.97142558808358 5.06477413208539 4.44019387157634 0.624580260509048
## S100A9  1.94570307882118 4.54632685012533 3.91764011030351 0.628686739821822
## S100A8  1.71364262382412 4.42085855610524 3.76199505612535 0.658863499979884
## HLA-DRA 2.09617672718907 3.72725848932921 3.12423457486866 0.603023914460557
## CD74    2.90288417785368 3.31555642893024 2.87554945366814 0.440006975262104
## CST3    1.49017608516901 2.95604474928493 2.28351996396206 0.672524785322875
##           p.value       FDR
##         <numeric> <numeric>
## LYZ             0         0
## S100A9          0         0
## S100A8          0         0
## HLA-DRA         0         0
## CD74            0         0
## CST3            0         0
```

We can plot the genes with the largest biological components, to verify that they are indeed highly variable (Figure \@ref(fig:hvgplot)).


```r
plotExpression(sce, features=rownames(top.dec)[1:10])
```

<div class="figure">
<img src="tenx_files/figure-html/hvgplot-1.png" alt="Distributions of normalized log-expression values for the top 10 genes with the largest biological components in the PBMC dataset. Each point represents the log-expression value in a single cell." width="100%"  class="widefigure" />
<p class="caption">(\#fig:hvgplot)Distributions of normalized log-expression values for the top 10 genes with the largest biological components in the PBMC dataset. Each point represents the log-expression value in a single cell.</p>
</div>

**Comments from Aaron:**

- The Poisson-based trend from `makeTechTrend()` tends to yield large biological components for highly-expressed genes for which Poisson noise is low (in the log-expression space).
This often includes so-called "house-keeping" genes coding for essential cellular components such as ribosomal proteins.
These genes are often considered uninteresting for characterizing cellular heterogeneity, 
though this is debatable as they are often differentially expressed in a variety of conditions [@glare2002betaactin;@nazari2015gapdh;@guimaraes2016patterns].
Indeed, the fact that they have large biological components indicates that there is strong variation in their expression across cells, which warrants some further investigation.
Nonetheless, if they are deemed to be uninteresting, their impact can be reduced by fitting the mean-variance trend to the endogenous genes.

# Dimensionality reduction

We use the `denoisePCA()` function with the assumed Poisson technical trend to choose the number of dimensions to retain after PCA.
Recall that this involves a random initialization when `approximate=TRUE`, which motivates the call to `set.seed()` to obtain reproducible results.


```r
set.seed(1000)
sce <- denoisePCA(sce, technical=new.trend, approximate=TRUE)
ncol(reducedDim(sce, "PCA"))
```

```
## [1] 12
```


```r
plot(attr(reducedDim(sce), "percentVar"), xlab="PC",
	ylab="Proportion of variance explained")
abline(v=ncol(reducedDim(sce, "PCA")), lty=2, col="red")
```

<div class="figure">
<img src="tenx_files/figure-html/screeplot-1.png" alt="Variance explained by each principal component in the PBMC dataset. The red line represents the chosen number of PCs." width="100%" />
<p class="caption">(\#fig:screeplot)Variance explained by each principal component in the PBMC dataset. The red line represents the chosen number of PCs.</p>
</div>

Examination of the first few PCs already reveals some strong substructure in the data (Figure \@ref(fig:pcaplot-init)).


```r
plotPCA(sce, ncomponents=3, colour_by="log10_total_features_by_counts")
```

<div class="figure">
<img src="tenx_files/figure-html/pcaplot-init-1.png" alt="Pairwise PCA plots of the first three PCs in the PBMC dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured by the log-number of expressed features." width="864" />
<p class="caption">(\#fig:pcaplot-init)Pairwise PCA plots of the first three PCs in the PBMC dataset, constructed from normalized log-expression values of genes with positive biological components. Each point represents a cell, coloured by the log-number of expressed features.</p>
</div>

This is recapitulated with a _t_-SNE plot (Figure \@ref(fig:tsneplot-init)).
Again, note that we set `use_dimred=` to perform _t_-SNE on the denoised expression matrix.


```r
set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30)
plotTSNE(sce, colour_by="log10_total_features_by_counts")
```

<div class="figure">
<img src="tenx_files/figure-html/tsneplot-init-1.png" alt="_t_-SNE plots constructed from the denoised PCs of the PBMC dataset. Each point represents a cell and is coloured according to the log-number of expressed features." width="100%" />
<p class="caption">(\#fig:tsneplot-init)_t_-SNE plots constructed from the denoised PCs of the PBMC dataset. Each point represents a cell and is coloured according to the log-number of expressed features.</p>
</div>

# Clustering with graph-based methods

We build a shared nearest neighbour graph [@xu2015identification] and use the Walktrap algorithm to identify clusters.


```r
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)
table(sce$Cluster)
```

```
## 
##   1   2   3   4   5   6   7   8   9  10  11  12  13 
## 580 807 521 202 516 128 817  40  18  45 125  82  36
```

We look at the ratio of the observed and expected edge weights to confirm that the clusters are modular.
(We don't look at the modularity score itself, as that varies by orders of magnitudes across clusters and is difficult to interpret.)
Figure \@ref(fig:clustermod) indicates that most of the clusters are well seperated, with few strong off-diagonal entries. 


```r
cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)

library(pheatmap)
pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
    color=colorRampPalette(c("white", "blue"))(100))
```

<div class="figure">
<img src="tenx_files/figure-html/clustermod-1.png" alt="Heatmap of the log~10~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links." width="100%" />
<p class="caption">(\#fig:clustermod)Heatmap of the log~10~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links.</p>
</div>

We examine the cluster identities on a _t_-SNE plot (Figure \@ref(fig:tsneplot-cluster)) to confirm that different clusters are indeed separated.


```r
plotTSNE(sce, colour_by="Cluster")
```

<div class="figure">
<img src="tenx_files/figure-html/tsneplot-cluster-1.png" alt="_t_-SNE plots constructed from the denoised PCs of the PBMC dataset. Each point represents a cell and is coloured according to its cluster identity." width="100%" />
<p class="caption">(\#fig:tsneplot-cluster)_t_-SNE plots constructed from the denoised PCs of the PBMC dataset. Each point represents a cell and is coloured according to its cluster identity.</p>
</div>

# Marker gene detection



We detect marker genes for each cluster using `findMarkers()`.
Again, we only look at upregulated genes in each cluster, as these are more useful for positive identification of cell types in a heterogeneous population.


```r
markers <- findMarkers(sce, clusters=sce$Cluster, direction="up")
```

We examine the markers for cluster 10 in more detail.
The upregulation of genes such as _PF4_ and _PPBP_ suggests that this cluster contains platelets or their precursors.


```r
marker.set <- markers[["10"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity
```

```
## DataFrame with 10 rows and 8 columns
##              Top              p.value                  FDR          logFC.1
##        <integer>            <numeric>            <numeric>        <numeric>
## PF4            1 2.05211320403009e-33 6.91439022965899e-29 6.75405849561143
## TAGLN2         1  1.3570056813532e-22 1.52409831425049e-18 4.40171343309397
## TMSB4X         2 2.01222579811262e-27 3.38999680208033e-23 2.91149334103303
## SDPR           2 1.24366719866146e-21 1.04760306479248e-17 5.40452070062898
## NRGN           3 1.54715587347124e-20  1.0425974000148e-16 4.87415412422472
## GPX1           3  1.0125517343013e-19 4.87384544793541e-16 4.96339461602238
## PPBP           3 2.26706652282989e-20 1.27310899033717e-16 6.26803629158601
## CCL5           6 3.84932059196734e-18 1.44110008917497e-14 4.53622129028549
## GNG11          6 3.19834491053025e-18 1.34706291769258e-14 5.26980617527753
## ACTB           7 9.86149100619106e-18 3.32273077962602e-14 3.33683896701068
##                 logFC.2          logFC.3          logFC.4          logFC.5
##               <numeric>        <numeric>        <numeric>        <numeric>
## PF4    6.70037626955874 6.75801307630927 6.75174806118932  6.7594697388179
## TAGLN2 4.65781231393239 4.62270543168491 4.60193611046034  4.3142234992163
## TMSB4X 2.55041744363595 2.93418266758885 3.17202467596703  3.5211828266668
## SDPR   5.35624658055225 5.41279757131019 5.41279757131019  5.4063321460251
## NRGN   4.64537299048304 4.88147623823525 4.88171598638508 4.87776037924028
## GPX1   2.70569235242282 5.11793862354292 5.18504246282525 4.68311296235489
## PPBP    6.1819369239228 6.27756356394969 6.27623332449874 6.26929041978115
## CCL5   5.14246588490684 1.46434277184655 2.93036613128338 5.09880481356499
## GNG11  5.21613623290634 5.27103963288141 5.27667639902105 5.23782451264031
## ACTB   2.57052019669099 3.14346374297953 2.95979979581139 3.71564066950398
```



This is confirmed in Figure \@ref(fig:heatmap), where the transcriptional profile of cluster 10 is clearly distinct from the others.


```r
chosen <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
    zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
    colour_columns_by="Cluster", columns=order(sce$Cluster),
    show_colnames=FALSE)
```

<div class="figure">
<img src="tenx_files/figure-html/heatmap-1.png" alt="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 10 in the PBMC dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend." width="100%"  class="widefigure" />
<p class="caption">(\#fig:heatmap)Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 10 in the PBMC dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend.</p>
</div>

# Concluding remarks

Having completed the basic analysis, we save the `SingleCellExperiment` object with its associated data to file.
This avoids having to repeat all of the pre-processing steps described above prior to further analyses.


```r
saveRDS(sce, file="pbmc_data.rds")
```

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
##  [1] pheatmap_1.0.10             scran_1.11.10              
##  [3] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.7.3            
##  [5] AnnotationFilter_1.7.0      GenomicFeatures_1.35.4     
##  [7] AnnotationDbi_1.45.0        scater_1.11.5              
##  [9] ggplot2_3.1.0               DropletUtils_1.3.4         
## [11] SingleCellExperiment_1.5.1  SummarizedExperiment_1.13.0
## [13] DelayedArray_0.9.4          matrixStats_0.54.0         
## [15] Biobase_2.43.0              GenomicRanges_1.35.1       
## [17] GenomeInfoDb_1.19.1         IRanges_2.17.3             
## [19] S4Vectors_0.21.8            BiocGenerics_0.29.1        
## [21] BiocParallel_1.17.3         bindrcpp_0.2.2             
## [23] BiocFileCache_1.7.0         dbplyr_1.2.2               
## [25] knitr_1.21                  BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##  [1] Rtsne_0.15               ggbeeswarm_0.6.0        
##  [3] colorspace_1.3-2         dynamicTreeCut_1.63-1   
##  [5] XVector_0.23.0           BiocNeighbors_1.1.7     
##  [7] bit64_0.9-7              codetools_0.2-16        
##  [9] R.methodsS3_1.7.1        Rsamtools_1.35.0        
## [11] R.oo_1.22.0              HDF5Array_1.11.10       
## [13] BiocManager_1.30.4       compiler_3.6.0          
## [15] httr_1.4.0               assertthat_0.2.0        
## [17] Matrix_1.2-15            lazyeval_0.2.1          
## [19] limma_3.39.3             htmltools_0.3.6         
## [21] prettyunits_1.0.2        tools_3.6.0             
## [23] igraph_1.2.2             gtable_0.2.0            
## [25] glue_1.3.0               GenomeInfoDbData_1.2.0  
## [27] reshape2_1.4.3           dplyr_0.7.8             
## [29] rappdirs_0.3.1           Rcpp_1.0.0              
## [31] Biostrings_2.51.1        rtracklayer_1.43.1      
## [33] DelayedMatrixStats_1.5.0 xfun_0.4                
## [35] stringr_1.3.1            irlba_2.3.3             
## [37] statmod_1.4.30           XML_3.98-1.16           
## [39] edgeR_3.25.2             zlibbioc_1.29.0         
## [41] scales_1.0.0             hms_0.4.2               
## [43] ProtGenerics_1.15.0      rhdf5_2.27.4            
## [45] RColorBrewer_1.1-2       yaml_2.2.0              
## [47] curl_3.2                 memoise_1.1.0           
## [49] gridExtra_2.3            biomaRt_2.39.2          
## [51] stringi_1.2.4            RSQLite_2.1.1           
## [53] highr_0.7                simpleSingleCell_1.7.8  
## [55] rlang_0.3.0.1            pkgconfig_2.0.2         
## [57] bitops_1.0-6             evaluate_0.12           
## [59] lattice_0.20-38          purrr_0.2.5             
## [61] Rhdf5lib_1.5.1           bindr_0.1.1             
## [63] GenomicAlignments_1.19.0 labeling_0.3            
## [65] cowplot_0.9.3            bit_1.1-14              
## [67] tidyselect_0.2.5         plyr_1.8.4              
## [69] magrittr_1.5             bookdown_0.9            
## [71] R6_2.3.0                 DBI_1.0.0               
## [73] pillar_1.3.1             withr_2.1.2             
## [75] RCurl_1.95-4.11          tibble_1.4.2            
## [77] crayon_1.3.4             rmarkdown_1.11          
## [79] viridis_0.5.1            progress_1.2.0          
## [81] locfit_1.5-9.1           grid_3.6.0              
## [83] blob_1.1.1               digest_0.6.18           
## [85] R.utils_2.7.0            munsell_0.5.0           
## [87] beeswarm_0.2.3           viridisLite_0.3.0       
## [89] vipor_0.4.5
```

# References
