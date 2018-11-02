---
title: Analyzing single-cell RNA-seq data containing UMI counts
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
date: "2018-11-02"
vignette: >
  %\VignetteIndexEntry{03. UMI count data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---



# Overview

In this workflow, we examine a heterogeneous dataset from a study of cell types in the mouse brain [@zeisel2015brain].
This contains approximately 3000 cells of varying types such as oligodendrocytes, microglia and neurons.
Individual cells were isolated using the Fluidigm C1 microfluidics system [@pollen2014lowcoverage] and library preparation was performed on each cell using a UMI-based protocol.
After sequencing, expression was quantified by counting the number of UMIs mapped to each gene.
Count data for all endogenous genes, mitochondrial genes and spike-in transcripts are available from http://linnarssonlab.org/cortex.


```r
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
base.url <- file.path("https://storage.googleapis.com",
    "linnarsson-lab-www-blobs/blobs/cortex")
mRNA.path <- bfcrpath(bfc, file.path(base.url, 
    "expression_mRNA_17-Aug-2014.txt"))
mito.path <- bfcrpath(bfc, file.path(base.url, 
    "expression_mito_17-Aug-2014.txt"))
spike.path <- bfcrpath(bfc, file.path(base.url, 
    "expression_spikes_17-Aug-2014.txt"))
```

# Setting up the data 

The count data are distributed across several files, so some work is necessary to consolidate them into a single matrix.
We define a simple utility function for loading data in from each file. 
(We stress that this function is only relevant to the current dataset, and should not be used for other datasets.
This kind of effort is generally not required if all of the counts are in a single file and separated from the metadata.)


```r
readFormat <- function(infile) { 
    # First column is empty.
    metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] 
    rownames(metadata) <- metadata[,1]
    metadata <- metadata[,-1]
    metadata <- as.data.frame(t(metadata))

    # First column after row names is some useless filler.
    counts <- read.delim(infile, stringsAsFactors=FALSE, 
        header=FALSE, row.names=1, skip=11)[,-1] 
    counts <- as.matrix(counts)
    return(list(metadata=metadata, counts=counts))
}
```

Using this function, we read in the counts for the endogenous genes, ERCC spike-in transcripts and mitochondrial genes.


```r
endo.data <- readFormat(mRNA.path)
spike.data <- readFormat(spike.path)
mito.data <- readFormat(mito.path)
```

We also need to rearrange the columns for the mitochondrial data, as the order is not consistent with the other files.


```r
m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)
mito.data$metadata <- mito.data$metadata[m,]
mito.data$counts <- mito.data$counts[,m]
```



In this particular dataset, some genes are represented by multiple rows corresponding to alternative genomic locations.
We sum the counts for all rows corresponding to a single gene for ease of interpretation.


```r
raw.names <- sub("_loc[0-9]+$", "", rownames(endo.data$counts))
new.counts <- rowsum(endo.data$counts, group=raw.names, reorder=FALSE)
endo.data$counts <- new.counts
```

The counts are then combined into a single matrix for constructing a `SingleCellExperiment` object.
For convenience, metadata for all cells are stored in the same object for later access.


```r
library(SingleCellExperiment)
all.counts <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)
sce <- SingleCellExperiment(list(counts=all.counts), colData=endo.data$metadata)
dim(sce)
```

```
## [1] 19896  3005
```

We add gene-based annotation identifying rows that correspond to each class of features.
We also determine the Ensembl identifier for each row.


```r
# Specifying the nature of each row.
nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))
is.spike <- rep(c(FALSE, FALSE, TRUE), nrows)
is.mito <- rep(c(FALSE, TRUE, FALSE), nrows)
isSpike(sce, "Spike") <- is.spike

# Adding Ensembl IDs.
library(org.Mm.eg.db)
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
rowData(sce)$ENSEMBL <- ensembl

sce
```

```
## class: SingleCellExperiment 
## dim: 19896 3005 
## metadata(0):
## assays(1): counts
## rownames(19896): Tspan12 Tshz1 ... ERCC-00170 ERCC-00171
## rowData names(1): ENSEMBL
## colnames(3005): V3 V4 ... V3006 V3007
## colData names(10): tissue group # ... level1class level2class
## reducedDimNames(0):
## spikeNames(1): Spike
```



# Quality control on the cells 

The original authors of the study have already removed low-quality cells prior to data publication.
Nonetheless, we compute some quality control metrics with *[scater](https://bioconductor.org/packages/3.9/scater)* [@mccarthy2017scater] to check whether the remaining cells are satisfactory.


```r
library(scater)
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 
```

We examine the distribution of the QC metrics across all cells (Figure \@ref(fig:libplotbrain)).
The library sizes here are at least one order of magnitude lower than observed in the 416B dataset.
This is consistent with the use of UMI counts rather than read counts, as each transcript molecule can only produce one UMI count but can yield many reads after fragmentation.
In addition, the spike-in proportions are more variable than observed in the 416B dataset.
This may reflect a greater variability in the total amount of endogenous RNA per cell when many cell types are present.


```r
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="", 
    breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Spike, xlab="ERCC proportion (%)",
    ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
    ylab="Number of cells", breaks=20, main="", col="grey80")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/libplotbrain-1.png" alt="Histograms of QC metrics including the library sizes, number of expressed genes and proportion of UMIs assigned to spike-in transcripts or mitochondrial genes for all cells in the brain dataset." width="100%"  class="widefigure" />
<p class="caption">(\#fig:libplotbrain)Histograms of QC metrics including the library sizes, number of expressed genes and proportion of UMIs assigned to spike-in transcripts or mitochondrial genes for all cells in the brain dataset.</p>
</div>

We remove small outliers for the library size and the number of expressed features, and large outliers for the spike-in proportions.
Again, the presence of spike-in transcripts means that we do not have to use the mitochondrial proportions.


```r
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce$pct_counts_Spike, nmads=3, type="higher")
```

Removal of low-quality cells is then performed by combining the filters for all of the metrics.
The majority of cells are retained, which suggests that the original quality control procedures were generally adequate.


```r
sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
    BySpike=sum(spike.drop), Remaining=ncol(sce))
```

```
##   ByLibSize ByFeature BySpike Remaining
## 1         8         3       8      2989
```

We could improve our cell filtering procedure further by setting `batch` in `isOutlier` to one or more known factors, e.g., mouse/plate of origin.
As previously mentioned, this would avoid inflation of the MAD and improve power to remove low-quality cells.
However, for simplicity, we will not do this as sufficient quality control has already been performed.



# Cell cycle classification

Application of `cyclone` [@scialdone2015computational] to the brain dataset suggests that most of the cells are in G1 phase (Figure \@ref(fig:phaseplotbrain)).
This requires the use of the Ensembl identifiers to match up with the pre-defined classifier.


```r
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
table(assignments$phase)
```

```
## 
##   G1  G2M    S 
## 2980    8    1
```

```r
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/phaseplotbrain-1.png" alt="Cell cycle phase scores from applying the pair-based classifier on the brain dataset, where each point represents a cell." width="100%" />
<p class="caption">(\#fig:phaseplotbrain)Cell cycle phase scores from applying the pair-based classifier on the brain dataset, where each point represents a cell.</p>
</div>

However, the intepretation of this result requires some caution due to differences between the training and test datasets.
The classifier was trained on C1 SMARTer data and accounts for the biases in that protocol. 
The brain dataset uses UMI counts, which has a different set of biases, e.g., 3'-end coverage only, no length bias, no amplification noise.
Furthermore, many neuronal cell types are expected to lie in the G0 resting phase, which is distinct from the other phases of the cell cycle [@coller2006new].
`cyclone` will generally assign such cells to the closest known phase in the training set, which would be G1.



# Examining gene-level metrics

Figure \@ref(fig:topgenebrain) shows the most highly expressed genes across the cell population in the brain dataset.
This is mostly occupied by spike-in transcripts, reflecting the use of spike-in concentrations that span the entire range of expression.
There are also a number of constitutively expressed genes, as expected.


```r
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotHighestExprs(sce, n=50) + fontsize
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/topgenebrain-1.png" alt="Percentage of total counts assigned to the top 50 most highly-abundant features in the brain dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature." width="100%"  class="widefigure" />
<p class="caption">(\#fig:topgenebrain)Percentage of total counts assigned to the top 50 most highly-abundant features in the brain dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature.</p>
</div>

Gene abundance is quantified by computing the average count across all cells (Figure \@ref(fig:abhistbrain)).
As previously mentioned, the UMI count is generally lower than the read count.


```r
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey",
    xlab=expression(Log[10]~"average count"))
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/abhistbrain-1.png" alt="Histogram of log-average counts for all genes in the brain dataset." width="100%" />
<p class="caption">(\#fig:abhistbrain)Histogram of log-average counts for all genes in the brain dataset.</p>
</div>

We save the average counts into the `SingleCellExperiment` object for later use.
We also remove genes that have average counts of zero, as this means that they are not expressed in any cell.


```r
rowData(sce)$ave.count <- ave.counts
to.keep <- ave.counts > 0
sce <- sce[to.keep,]
summary(to.keep)
```

```
##    Mode   FALSE    TRUE 
## logical       2   19894
```



# Normalization of cell-specific biases

For endogenous genes, normalization is performed using the `computeSumFactors` function as previously described.
Here, we cluster similar cells together and normalize the cells in each cluster using the deconvolution method [@lun2016pooling].
This improves normalization accuracy by reducing the number of DE genes between cells in the same cluster.
Scaling is then performed to ensure that size factors of cells in different clusters are comparable.


```r
set.seed(1000)
clusters <- quickCluster(sce, min.mean=0.1, method="igraph")
sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)
summary(sizeFactors(sce))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1255  0.4633  0.8207  1.0000  1.3366  4.7022
```

We use a average count threshold of 0.1 to define high-abundance genes to use during normalization.
This is lower than the default threshold of `min.mean=1` in `computeSumFactors`, reflecting the fact that UMI counts are generally smaller than read counts.



Compared to the 416B analysis, more scatter is observed around the trend between the total count and size factor for each cell (Figure \@ref(fig:normplotbrain)).
This is consistent with an increased amount of DE between cells of different types, which compromises the accuracy of library size normalization [@robinson2010scaling].
In contrast, the size factors are estimated based on median ratios and are more robust to the presence of DE between cells.


```r
plot(sizeFactors(sce), sce$total_counts/1e3, log="xy",
    ylab="Library size (thousands)", xlab="Size factor")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/normplotbrain-1.png" alt="Size factors from deconvolution, plotted against library sizes for all cells in the brain dataset. Axes are shown on a log-scale." width="100%" />
<p class="caption">(\#fig:normplotbrain)Size factors from deconvolution, plotted against library sizes for all cells in the brain dataset. Axes are shown on a log-scale.</p>
</div>

We also compute size factors specific to the spike-in set, as previously described.


```r
sce <- computeSpikeFactors(sce, type="Spike", general.use=FALSE)
```

Finally, normalized log-expression values are computed for each endogenous gene or spike-in transcript using the appropriate size factors.


```r
sce <- normalize(sce)
```



__Comments from Aaron:__

- Only a rough clustering is required to avoid pooling together very different cell types in `computeSumFactors()`.
This reduces the chance of violating the non-DE assumption that is made during any gene-based scaling normalization.
We stress that there is no need for precise clustering at this step, as we will not be interpreting these clusters at all.
Our normalization strategy is also robust to a moderate level of differential expression between cells in the same cluster, so careful definition of subclusters is not required.
We do ask that there are sufficient cells in each cluster for pooling, which can be guaranteed with the `min.size=` argument in `quickCluster()`.
- `quickCluster` uses distances based on Spearman's rank correlation for clustering.
This ensures that scaling biases in the counts do not affect clustering, but yields very coarse clusters and is not recommended for biological interpretation.
However, for the purposes of breaking up distinct cell types, this is more than sufficient.
- For large datasets, using `method="igraph"` in `quickCluster` will speed up clustering. 
This uses a graph-based clustering algorithm with a randomization step, hence the need for `set.seed()`.
See `?buildSNNGraph` for more details.

# Modelling and removing technical noise

We model the technical noise by fitting a mean-variance trend to the spike-in transcripts, as previously described.
In theory, we should block on the plate of origin for each cell.
However, only 20-40 cells are available on each plate, and the population is also highly heterogeneous.
This means that we cannot assume that the distribution of sampled cell types on each plate is the same.
Thus, to avoid regressing out potential biology, we will not block on any factors in this analysis.


```r
var.fit <- trendVar(sce, parametric=TRUE, loess.args=list(span=0.4))
var.out <- decomposeVar(sce, var.fit)
```

Figure \@ref(fig:hvgplotbrain) indicates that the trend is fitted accurately to the technical variances.
The technical and total variances are also much smaller than those in the 416B dataset.
This is due to the use of UMIs, which reduces the noise caused by variable PCR amplification [@islam2014quantitative].
Furthermore, the spike-in trend is consistently lower than the variances of the endogenous genes.
This reflects the heterogeneity in gene expression across cells of different types.


```r
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
points(var.out$mean[isSpike(sce)], var.out$total[isSpike(sce)], col="red", pch=16)
curve(var.fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/hvgplotbrain-1.png" alt="Variance of normalized log-expression values against the mean for each gene, calculated across all cells in the brain dataset after blocking on the sex effect. The blue line represents the mean-dependent trend in the technical variance of the spike-in transcripts (also highlighted as red points)." width="100%" />
<p class="caption">(\#fig:hvgplotbrain)Variance of normalized log-expression values against the mean for each gene, calculated across all cells in the brain dataset after blocking on the sex effect. The blue line represents the mean-dependent trend in the technical variance of the spike-in transcripts (also highlighted as red points).</p>
</div>

We check the distribution of expression values for the genes with the largest biological components to ensure that they are not driven by outliers (Figure \@ref(fig:hvgvioplotbrain)).
Some tweaking of the `plotExpression` parameters is necessary to visualize a large number of cells.


```r
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, rownames(var.out)[chosen.genes], 
    point_alpha=0.05, jitter_type="jitter") + fontsize
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/hvgvioplotbrain-1.png" alt="Violin plots of normalized log-expression values for the top 10 HVGs in the brain dataset. For each gene, each point represents the log-expression value for an individual cell." width="100%" />
<p class="caption">(\#fig:hvgvioplotbrain)Violin plots of normalized log-expression values for the top 10 HVGs in the brain dataset. For each gene, each point represents the log-expression value for an individual cell.</p>
</div>

Finally, we use PCA to denoise the expression values, yielding a set of coordinates for each cell where the technical noise has been removed.
Setting `approximate=TRUE` in `denoisePCA` will perform an approximate singular value decomposition (SVD), using methods from the *[irlba](https://CRAN.R-project.org/package=irlba)* package.
This is much faster than the exact algorithm on large datasets without much loss of accuracy.
The approximate algorithm involves a random initialization so we set the seed to guarantee reproducibility.


```r
set.seed(1000)
sce <- denoisePCA(sce, technical=var.fit$trend, approximate=TRUE)
ncol(reducedDim(sce, "PCA"))
```

```
## [1] 100
```

**Comments from Aaron:**

- The upper limit of PCs to retain in `denoisePCA()` is specified by the `max.rank=` argument.
This is set to 100 by default to ensure that the approximate SVD runs quickly.
A higher `max.rank` may be more appropriate for extremely heterogeneous populations, though the default setting is generally satisfactory for dimensionality reduction.



# Data exploration with dimensionality reduction

We perform dimensionality reduction on the denoised PCs to check if there is any substructure. 
Cells separate into clear clusters in the _t_-SNE plot [@van2008visualizing] in Figure \@ref(fig:tsneplotbrain), corresponding to distinct subpopulations.
This is consistent with the presence of multiple cell types in the diverse brain population.
We increase the perplexity to favour visualization of the overall structure at the expense of local scale.


```r
set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=50)
tsne1 <- plotTSNE(sce, colour_by="Neurod6") + fontsize
tsne2 <- plotTSNE(sce, colour_by="Mog") + fontsize
multiplot(tsne1, tsne2, cols=2)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/tsneplotbrain-1.png" alt="_t_-SNE plots constructed from the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its expression of _Neurod6_ (left) or _Mog_ (right)." width="1152" />
<p class="caption">(\#fig:tsneplotbrain)_t_-SNE plots constructed from the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its expression of _Neurod6_ (left) or _Mog_ (right).</p>
</div>

The PCA plot is less effective at separating cells into many different clusters (Figure \@ref(fig:pcaplotbrain)).
This is because the first two PCs are driven by strong differences between specific subpopulations, which reduces the resolution of more subtle differences between some of the other subpopulations.
Nonetheless, some substructure is still visible.


```r
pca1 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Neurod6") + fontsize
pca2 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Mog") + fontsize
multiplot(pca1, pca2, cols=2)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/pcaplotbrain-1.png" alt="PCA plots constructed from the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its expression of the _Neurod6_ (left) or _Mog_ (right)." width="1152" />
<p class="caption">(\#fig:pcaplotbrain)PCA plots constructed from the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its expression of the _Neurod6_ (left) or _Mog_ (right).</p>
</div>

For both methods, we colour each cell based on the expression of a particular gene.
This is a useful strategy for visualizing changes in expression across the lower-dimensional space.
It can also be used to characterise each cluster if the selected genes are known markers for particular cell types.
For example, _Mog_ can be used to identify clusters corresponding to oligodendrocytes.



# Clustering cells into putative subpopulations

## Using graph-based clustering

The reduced dimension coordinates are used to cluster cells into putative subpopulations.
We do so by constructing a shared-nearest-neighbour graph [@xu2015identification], in which cells are the nodes and edges are formed between cells that share nearest neighbours.
Clusters are then defined as highly connected communities of cells within this graph, using methods from the *[igraph](https://CRAN.R-project.org/package=igraph)* package.
This is more efficient than forming a pairwise distance matrix for hierarchical clustering of large numbers of cells.


```r
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
cluster.out <- igraph::cluster_walktrap(snn.gr)
my.clusters <- cluster.out$membership
table(my.clusters)
```

```
## my.clusters
##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
## 201 186 129 230 105 414 171 643 395  82  83 213  28  36  12  61
```

We visualize the cluster assignments for all cells on the _t_-SNE plot in Figure \@ref(fig:tsneclusterbrain).
Adjacent cells are generally assigned to the same cluster, indicating that the clustering procedure was applied correctly.


```r
sce$cluster <- factor(my.clusters)
plotTSNE(sce, colour_by="cluster") + fontsize
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/tsneclusterbrain-1.png" alt="_t_-SNE plot of the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its assigned cluster identity." width="100%" />
<p class="caption">(\#fig:tsneclusterbrain)_t_-SNE plot of the denoised PCs of the brain dataset. Each point represents a cell and is coloured according to its assigned cluster identity.</p>
</div>



An alternative approach is to use graph-based visualizations such as force-directed layouts (Figure \@ref(fig:fdlbrain)).
These are appealing as they directly represent the relationships used during clustering.
However, convergence tends to be slow for large graphs, so some tinkering with `niter=` may be required to ensure that the results are stable. 


```r
set.seed(2000)
reducedDim(sce, "force") <- igraph::layout_with_fr(snn.gr, niter=5000)
plotReducedDim(sce, colour_by="cluster", use_dimred="force")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/fdlbrain-1.png" alt="Force-directed layout for the shared nearest-neighbour graph of the brain dataset. Each point represents a cell and is coloured according to its assigned cluster identity." width="100%" />
<p class="caption">(\#fig:fdlbrain)Force-directed layout for the shared nearest-neighbour graph of the brain dataset. Each point represents a cell and is coloured according to its assigned cluster identity.</p>
</div>

Very heterogeneous datasets may yield a few large clusters on the first round of clustering.
It can be useful to repeat the variance modelling, denoising and clustering using only the cells within each of the initial clusters.
This can be achieved by subsetting `sce` according to a particular level of `my.clusters`, and re-applying the relevant functions on the subset.
Doing so may focus on a different set of genes that define heterogeneity _within_ an initial cluster, as opposed to those that define differences _between_ the initial clusters.
This would allow fine-scale structure within each cluster to be explored at greater resolution. 
For simplicity, though, we will only use the broad clusters corresponding to clear subpopulations in this workflow.

**Comments from Aaron:**

- Many different clustering methods are available in the *[igraph](https://CRAN.R-project.org/package=igraph)* package.
We find that the Walktrap algorithm is usually a good default choice [@yang2016comparative], though users are encouraged to experiment with different algorithms.
- Decreasing the number of neighbours `k` in `buildSNNGraph` will reduce the connectivity of the graph.
This will generally result in the formation of smaller clusters [@xu2015identification], which may be desirable if greater resolution is required.
- Notice that we do not run `library(igraph)`, but instead use `igraph::` to extract methods from the package. 
This is because *[igraph](https://CRAN.R-project.org/package=igraph)* contains a `normalize` method that will override its counterpart from *[scater](https://bioconductor.org/packages/3.9/scater)*, resulting in some unusual bugs.

## Evaluating graph-based clusters

The modularity score provides a global measure of clustering performance for community detection methods.
Briefly, it compares the number of within-cluster edges to the expected number under a null model of random edges.
A high modularity score (approaching the maximum of 1) indicates that the detected clusters are enriched for internal edges, with relatively few edges between clusters.


```r
igraph::modularity(cluster.out)
```

```
## [1] 0.7680668
```

We further investigate the clusters by examining the total weight of edges for each pair of clusters.
For each pair, the observed total weight is compared to what is expected under a null model, similar to the modularity calculation.
Most clusters contain more internal links than expected (Figure \@ref(fig:heatmodbrain)), while links between clusters are fewer than expected.
This indicates that we successfully clustered cells into highly-connected communities.


```r
mod.out <- clusterModularity(snn.gr, my.clusters, get.values=TRUE)
ratio <- mod.out$observed/mod.out$expected
lratio <- log10(ratio + 1)

library(pheatmap)
pheatmap(lratio, cluster_rows=FALSE, cluster_cols=FALSE, 
    color=colorRampPalette(c("white", "blue"))(100))
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/heatmodbrain-1.png" alt="Heatmap of the log~10~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links." width="100%" />
<p class="caption">(\#fig:heatmodbrain)Heatmap of the log~10~-ratio of the total weight between nodes in the same cluster or in different clusters, relative to the total weight expected under a null model of random links.</p>
</div>

To summarize the relationships between clusters, we use the ratio of the observed and expected total weights to build a graph across clusters.
The cluster-based graph can be visualized using a force-directed layout to identify "clusters of clusters" that are highly interconnected. 
This is similar to the "graph abstraction" strategy proposed by @wolf2017graph.


```r
cluster.gr <- igraph::graph_from_adjacency_matrix(ratio, 
    mode="undirected", weighted=TRUE, diag=FALSE)
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*10)  
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/graphbrain-1.png" alt="Force-directed layout showing the relationships between clusters based on the ratio of observed to expected total weights between nodes in different clusters. The thickness of the edge between a pair of clusters is proportional to the corresponding ratio." width="100%" />
<p class="caption">(\#fig:graphbrain)Force-directed layout showing the relationships between clusters based on the ratio of observed to expected total weights between nodes in different clusters. The thickness of the edge between a pair of clusters is proportional to the corresponding ratio.</p>
</div>

**Comments from Aaron:**

- We do not use the silhouette coefficient to assess clustering for large datasets.
This is because `cluster::silhouette` requires the construction of a distance matrix, which may not be feasible when many cells are involved.
- Technically, the modularity score is obtained by subtracting the observed from expected total weights. 
We use the ratio instead as this is guaranteed to be positive and does not exhibit differences in magnitude due to differences in the number of cells in each cluster.

# Detecting marker genes between subpopulations

We use the `findMarkers` function with `direction="up"` to identify upregulated marker genes for each cluster.
As previously mentioned, we focus on upregulated genes as these can quickly provide positive identification of cell type in a heterogeneous population.
We examine the table for cluster 1, in which log-fold changes are reported between cluster 1 and every other cluster.
The same output is provided for each cluster in order to identify genes that discriminate between clusters.




```r
markers <- findMarkers(sce, my.clusters, direction="up")
marker.set <- markers[["1"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity
```

```
## DataFrame with 10 rows and 8 columns
##               Top               p.value                   FDR
##         <integer>             <numeric>             <numeric>
## Snap25          1  1.5867064696874e-268 3.14754962391882e-264
## Mllt11          1 2.31447867657062e-196 9.18246270142618e-193
## Atp1a3          1 1.88672817822383e-177 5.34671812448942e-174
## Gad1            1 4.97164699225956e-176 9.86225613854515e-173
## Celf4           1 5.46552707358057e-160 6.02331447547864e-157
## Gad2            1 2.49777876244356e-155  2.4774218655296e-152
## Vstm2a          1 1.17315765416478e-111  4.4753708433975e-109
## Synpr           1  2.97884331198009e-70  3.17695240751335e-68
## Slc32a1         1  1.30132881331111e-66  1.18960643638953e-64
## Ndrg4           2 2.36919621691788e-249 2.34988726775007e-245
##                    logFC.2           logFC.3          logFC.4
##                  <numeric>         <numeric>        <numeric>
## Snap25  -0.354294095149351 0.340117299504504 3.65521123984669
## Mllt11   0.365170419908248 0.924996225036646 2.98325911354461
## Atp1a3   0.462504066481883  1.12674016284263 3.23363921881298
## Gad1      4.18705718512123  3.57360914317586 3.98791017467217
## Celf4   -0.187435032806592 0.704927320176858 2.70010040039546
## Gad2      3.81079423774551  3.28063593061477 3.63281114318632
## Vstm2a    2.14747846790715  2.35624004184546 2.89937207148961
## Synpr      2.9288511065565  2.48276527902932 3.08673521450369
## Slc32a1   1.72021427539416  1.57696733397569 1.70655532822886
## Ndrg4    0.251601283679209  1.05192547512691  3.6272374874884
##                     logFC.5            logFC.6
##                   <numeric>          <numeric>
## Snap25      1.2457243604206   0.79590537709112
## Mllt11     1.49884172626987  0.622323985512458
## Atp1a3  -0.0305212989553296 -0.115321746617353
## Gad1       3.91590008743585   4.11434816863818
## Celf4     0.430907554882714  0.195423912602373
## Gad2       3.45750199364546   3.76682908641341
## Vstm2a     2.68316613092912   2.79126151922907
## Synpr      3.02164798455159   3.08438332201148
## Slc32a1    1.61148118807856   1.69842283940163
## Ndrg4     0.965091775807663  0.806032578350594
```



We save the list of candidate marker genes for further examination, using compression to reduce the file size.
The `overlapExprs` function may also be useful for summarizing differences between clusters, as previously mentioned.


```r
gzout <- gzfile("brain_marker_1.tsv.gz", open="wb")
write.table(marker.set, file=gzout, sep="\t", quote=FALSE, col.names=NA)
close(gzout)
```

Figure \@ref(fig:heatmapmarkerbrain) indicates that most of the top markers are strongly DE in cells of cluster 1 compared to some or all of the other clusters.
We can use these markers to identify cells from cluster 1 in validation studies with an independent population of cells.
A quick look at the markers suggest that cluster 1 represents interneurons based on expression of *Gad1* and *Slc6a1* [@zeng2012largescale],
differing from closely related cells in cluster 10 by virtue of high *Synpr* expression.


```r
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=top.markers, columns=order(my.clusters),
    colour_columns_by="cluster", cluster_cols=FALSE, 
    center=TRUE, symmetric=TRUE, zlim=c(-5, 5))
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/work-2-umis_files/figure-html/heatmapmarkerbrain-1.png" alt="Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 1 in the brain dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend." width="100%"  class="widefigure" />
<p class="caption">(\#fig:heatmapmarkerbrain)Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 1 in the brain dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend.</p>
</div>

# Concluding remarks

Having completed the basic analysis, we save the `SingleCellExperiment` object with its associated data to file.
This is especially important here as the brain dataset is quite large.
If further analyses are to be performed, it would be inconvenient to have to repeat all of the pre-processing steps described above.


```r
saveRDS(file="brain_data.rds", sce)
```



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
##  [1] pheatmap_1.0.10             scran_1.11.1               
##  [3] scater_1.11.1               ggplot2_3.1.0              
##  [5] org.Mm.eg.db_3.7.0          AnnotationDbi_1.45.0       
##  [7] SingleCellExperiment_1.5.0  SummarizedExperiment_1.13.0
##  [9] DelayedArray_0.9.0          BiocParallel_1.17.0        
## [11] matrixStats_0.54.0          Biobase_2.43.0             
## [13] GenomicRanges_1.35.0        GenomeInfoDb_1.19.0        
## [15] IRanges_2.17.0              S4Vectors_0.21.0           
## [17] BiocGenerics_0.29.0         bindrcpp_0.2.2             
## [19] BiocFileCache_1.7.0         dbplyr_1.2.2               
## [21] knitr_1.20                  BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] RColorBrewer_1.1-2       httr_1.3.1              
##  [5] rprojroot_1.3-2          dynamicTreeCut_1.63-1   
##  [7] tools_3.6.0              backports_1.1.2         
##  [9] R6_2.3.0                 irlba_2.3.2             
## [11] HDF5Array_1.11.0         vipor_0.4.5             
## [13] DBI_1.0.0                lazyeval_0.2.1          
## [15] colorspace_1.3-2         withr_2.1.2             
## [17] tidyselect_0.2.5         gridExtra_2.3           
## [19] bit_1.1-14               curl_3.2                
## [21] compiler_3.6.0           BiocNeighbors_1.1.0     
## [23] labeling_0.3             bookdown_0.7            
## [25] scales_1.0.0             rappdirs_0.3.1          
## [27] stringr_1.3.1            digest_0.6.18           
## [29] rmarkdown_1.10           XVector_0.23.0          
## [31] pkgconfig_2.0.2          htmltools_0.3.6         
## [33] limma_3.39.0             highr_0.7               
## [35] rlang_0.3.0.1            RSQLite_2.1.1           
## [37] DelayedMatrixStats_1.5.0 bindr_0.1.1             
## [39] dplyr_0.7.7              RCurl_1.95-4.11         
## [41] magrittr_1.5             GenomeInfoDbData_1.2.0  
## [43] Matrix_1.2-15            Rcpp_0.12.19            
## [45] ggbeeswarm_0.6.0         munsell_0.5.0           
## [47] Rhdf5lib_1.5.0           viridis_0.5.1           
## [49] stringi_1.2.4            yaml_2.2.0              
## [51] edgeR_3.25.0             zlibbioc_1.29.0         
## [53] rhdf5_2.27.0             Rtsne_0.13              
## [55] plyr_1.8.4               grid_3.6.0              
## [57] blob_1.1.1               crayon_1.3.4            
## [59] lattice_0.20-35          cowplot_0.9.3           
## [61] locfit_1.5-9.1           pillar_1.3.0            
## [63] igraph_1.2.2             reshape2_1.4.3          
## [65] glue_1.3.0               evaluate_0.12           
## [67] BiocManager_1.30.3       gtable_0.2.0            
## [69] purrr_0.2.5              assertthat_0.2.0        
## [71] xfun_0.4                 viridisLite_0.3.0       
## [73] tibble_1.4.2             beeswarm_0.2.3          
## [75] memoise_1.1.0            statmod_1.4.30
```

# References

