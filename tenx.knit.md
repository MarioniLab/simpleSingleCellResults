---
title: Analyzing single-cell RNA sequencing data from droplet-based protocols
author: 
- name: Aaron T. L. Lun
  affiliation: Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
date: "2019-02-08"
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
## [1] 4237
```

We then subset our `SingleCellExperiment` object to retain only the detected cells.


```r
# using which() to automatically remove NAs.
sce <- sce[,which(e.out$FDR <= 0.001)]
```

**Comments from Aaron:**

- `emptyDrops()` computes Monte Carlo $p$-values based on a Dirichlet-multinomial model of sampling molecules into droplets.
These $p$-values are stochastic so it is important to set the random seed to obtain reproducible results.
- The stability of the Monte Carlo $p$-values depends on the number of iterations used to compute them. 
This is controlled using the `niters=` argument in `emptyDrops()`, set to a default of 10000 for speed.
Larger values improve stability at the cost of time only, so users should set `niters=` to the largest value they are willing to tolerate.
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
##   FALSE   927    0
##   TRUE   1832 2530
```

As mentioned above, `emptyDrops()` assumes that barcodes with low total UMI counts are empty droplets.
Thus, the null hypothesis should be true for all of these barcodes. 
We can check whether the hypothesis test holds its size by examining the distribution of $p$-values for low-total barcodes.
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
<img src="tenx_files/figure-html/ambientpvalhist-1.png" alt="Distribution of $p$-values for the assumed empty droplets." width="100%" />
<p class="caption">(\#fig:ambientpvalhist)Distribution of $p$-values for the assumed empty droplets.</p>
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
## logical    3927     310
```

**Comments from Aaron:**

- The above justification for using a more relaxed filter is largely retrospective.
In practice, we may not know _a priori_ the degree of population heterogeneity and whether it manifests in the QC metrics.
We recommend performing an initial analysis with some QC, and then relaxing the filter (or making it more stringent) based on further diagnostics.
See [here](https://bioconductor.org/packages/3.9/simpleSingleCell/vignettes/qc.html#checking-for-discarded-cell-types) for an example of a potential diagnostic approach.
- Removal of cells with high mitochondrial content assumes that the damage to the cell membrane is modest enough that the mitochondria are retained.
It is possible for the cells to be so damaged that all cytoplasmic content is lost and only the stripped nucleus remains for sequencing.
This manifests as mitochondrial proportions of zero, usually accompanied by low library sizes/numbers of expressed genes and possibly low ribosomal protein gene expression.
If a cluster of such cells is observed in the data, they can be removed by using `isOutlier()` on the mitochondrial proportions with `type="lower"` and `log=TRUE`
(to improve resolution around zero).

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
##    1    2    3    4    5    6    7    8 
##  372  674  503  616  376 1099  160  127
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
##  0.000983  0.709306  0.870530  1.000000  1.094430 11.954870
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
fit$trend <- new.trend # overwrite trend.
dec <- decomposeVar(fit=fit) # use per-gene variance estimates in 'fit'.
top.dec <- dec[order(dec$bio, decreasing=TRUE),] 
head(top.dec)
```

```
## DataFrame with 6 rows and 6 columns
##                     mean            total              bio              tech
##                <numeric>        <numeric>        <numeric>         <numeric>
## LYZ     1.98439137942878 5.14036911600168 4.51227338997385 0.628095726027826
## S100A9  1.95891443219613 4.65018654869279  4.0181531762503 0.632033372442491
## S100A8  1.72794462749961  4.5172745961365 3.85603339530558 0.661241200830916
## HLA-DRA 2.10230972005318 3.75630516417143 3.14783809568087 0.608467068490557
## CD74    2.90443667945225 3.33713280593544  2.8834466304442 0.453686175491244
## CST3     1.5010372946582 3.02462480027789 2.34984258954562 0.674782210732269
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
## [1] 14
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
## 807 585 203 510 516  47 126 790  46  39 152  92  14
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

We examine the markers for cluster 9 in more detail.
The upregulation of genes such as _PF4_ and _PPBP_ suggests that this cluster contains platelets or their precursors.


```r
marker.set <- markers[["9"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity
```

```
## DataFrame with 10 rows and 8 columns
##                 Top              p.value                  FDR
##           <integer>            <numeric>            <numeric>
## PF4               1 1.70721625716368e-29 5.75229445688731e-25
## TAGLN2            1 2.99166656659263e-24 3.36004044315907e-20
## TMSB4X            2  4.4627578572629e-29 7.51840816213079e-25
## SDPR              2 1.22530295149563e-21 1.03213394119235e-17
## B2M               3 6.89735766236593e-21 3.87332615126263e-17
## NRGN              3 3.91954244119406e-21 2.64130126027185e-17
## GPX1              4 9.92135077150974e-20 4.17862491119062e-16
## PPBP              4 5.57911472831314e-20 2.68546702365403e-16
## HIST1H2AC         6  1.6059986319309e-18 5.41125179042799e-15
## ACTB              7    6.05048064744e-19 2.26516549927605e-15
##                    logFC.1          logFC.2          logFC.3
##                  <numeric>        <numeric>        <numeric>
## PF4       7.25679856477077  7.3045287421145  7.3042842273026
## TAGLN2    5.22185788080801 5.01280635916981 5.17647733126683
## TMSB4X    3.12406328360886 3.47124463540495 3.76944420404423
## SDPR      5.97127253659018 6.01775614734224 6.02432631992538
## B2M       2.18320635116997 1.66648176326608 1.65040834974051
## NRGN      5.21679045907656 5.44767004267579 5.45490583411146
## GPX1      3.24738954031262 5.55543584399255 5.74864823420056
## PPBP      6.81831106062476 6.90236145593299 6.91021363042703
## HIST1H2AC 5.89203212490953 5.92078959651487 5.94674603085075
## ACTB      3.12401577813728 3.91172509338832 3.52313142155693
##                    logFC.4          logFC.5
##                  <numeric>        <numeric>
## PF4       7.31034660245878 7.31189287548287
## TAGLN2    5.19557908201902 4.89096713081308
## TMSB4X    3.54096333125848 4.11450107743419
## SDPR      6.02432631992538 6.01800803163108
## B2M       1.49777214868441 2.38160588954119
## NRGN      5.45453680428722 5.45084540253421
## GPX1      5.68135911748578 5.23928445630973
## PPBP      6.91109566740711 6.90300646072055
## HIST1H2AC 5.92234095046611 5.84781580009597
## ACTB      3.72307886150859 4.28800048032883
```



This is confirmed in Figure \@ref(fig:heatmap), where the transcriptional profile of cluster 9 is clearly distinct from the others.


```r
chosen <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
    zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
    colour_columns_by="Cluster", columns=order(sce$Cluster),
    show_colnames=FALSE)
```

<div class="figure">
<img src="tenx_files/figure-html/heatmap-1.png" alt="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 9 in the PBMC dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend." width="100%"  class="widefigure" />
<p class="caption">(\#fig:heatmap)Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 9 in the PBMC dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend.</p>
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
## R Under development (unstable) (2019-01-14 r75992)
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
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] pheatmap_1.0.12             scran_1.11.20              
##  [3] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.7.8            
##  [5] AnnotationFilter_1.7.0      GenomicFeatures_1.35.6     
##  [7] AnnotationDbi_1.45.0        scater_1.11.11             
##  [9] ggplot2_3.1.0               DropletUtils_1.3.9         
## [11] SingleCellExperiment_1.5.2  SummarizedExperiment_1.13.0
## [13] DelayedArray_0.9.8          BiocParallel_1.17.9        
## [15] matrixStats_0.54.0          Biobase_2.43.1             
## [17] GenomicRanges_1.35.1        GenomeInfoDb_1.19.1        
## [19] IRanges_2.17.4              S4Vectors_0.21.10          
## [21] BiocGenerics_0.29.1         bindrcpp_0.2.2             
## [23] BiocFileCache_1.7.0         dbplyr_1.3.0               
## [25] knitr_1.21                  BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##  [1] Rtsne_0.15               ggbeeswarm_0.6.0        
##  [3] colorspace_1.4-0         dynamicTreeCut_1.63-1   
##  [5] XVector_0.23.0           BiocNeighbors_1.1.11    
##  [7] bit64_0.9-7              codetools_0.2-16        
##  [9] R.methodsS3_1.7.1        Rsamtools_1.35.2        
## [11] R.oo_1.22.0              HDF5Array_1.11.10       
## [13] BiocManager_1.30.4       compiler_3.6.0          
## [15] httr_1.4.0               assertthat_0.2.0        
## [17] Matrix_1.2-15            lazyeval_0.2.1          
## [19] limma_3.39.5             BiocSingular_0.99.0     
## [21] htmltools_0.3.6          prettyunits_1.0.2       
## [23] tools_3.6.0              rsvd_1.0.0              
## [25] igraph_1.2.2             gtable_0.2.0            
## [27] glue_1.3.0               GenomeInfoDbData_1.2.0  
## [29] reshape2_1.4.3           dplyr_0.7.8             
## [31] rappdirs_0.3.1           Rcpp_1.0.0              
## [33] Biostrings_2.51.2        rtracklayer_1.43.1      
## [35] DelayedMatrixStats_1.5.2 xfun_0.4                
## [37] stringr_1.3.1            ps_1.3.0                
## [39] irlba_2.3.3              statmod_1.4.30          
## [41] XML_3.98-1.16            edgeR_3.25.3            
## [43] zlibbioc_1.29.0          scales_1.0.0            
## [45] hms_0.4.2                ProtGenerics_1.15.0     
## [47] rhdf5_2.27.9             RColorBrewer_1.1-2      
## [49] yaml_2.2.0               curl_3.3                
## [51] memoise_1.1.0            gridExtra_2.3           
## [53] biomaRt_2.39.2           stringi_1.2.4           
## [55] RSQLite_2.1.1            highr_0.7               
## [57] simpleSingleCell_1.7.16  rlang_0.3.1             
## [59] pkgconfig_2.0.2          bitops_1.0-6            
## [61] evaluate_0.12            lattice_0.20-38         
## [63] purrr_0.3.0              Rhdf5lib_1.5.1          
## [65] bindr_0.1.1              GenomicAlignments_1.19.1
## [67] labeling_0.3             cowplot_0.9.4           
## [69] bit_1.1-14               tidyselect_0.2.5        
## [71] processx_3.2.1           plyr_1.8.4              
## [73] magrittr_1.5             bookdown_0.9            
## [75] R6_2.3.0                 DBI_1.0.0               
## [77] pillar_1.3.1             withr_2.1.2             
## [79] RCurl_1.95-4.11          tibble_2.0.1            
## [81] crayon_1.3.4             rmarkdown_1.11          
## [83] viridis_0.5.1            progress_1.2.0          
## [85] locfit_1.5-9.1           grid_3.6.0              
## [87] blob_1.1.1               callr_3.1.1             
## [89] digest_0.6.18            R.utils_2.7.0           
## [91] munsell_0.5.0            beeswarm_0.2.3          
## [93] viridisLite_0.3.0        vipor_0.4.5
```

# References
