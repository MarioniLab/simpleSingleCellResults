---
title: Characterizing gene expression variance in single-cell RNA-seq data
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
date: "2018-11-10"
vignette: >
  %\VignetteIndexEntry{09. Advanced variance modelling}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---
    


# Overview

Here, we describe some more detailed aspects of modelling variability in gene expression across the cell population.
This includes the detection of significant highly variable genes (HVGs), 
as well as advanced modelling of the mean-variance trend in the presence of confounding factors.

# Detecting highly variable genes

## Overview

HVGs are defined as genes with biological components that are significantly greater than zero.
These genes are interesting as they drive differences in the expression profiles between cells, and should be prioritized for further investigation.
Formal detection of HVGs allows us to avoid genes that are highly variable due to technical factors such as sampling noise during RNA capture and library preparation.
This adds another level of statistical rigour to our previous analyses, in which we only modelled the technical component.

## Setting up the data 

### Loading the dataset

To demonstrate, we use data from haematopoietic stem cells (HSCs) [@wilson2015combined], generated using the Smart-seq2 protocol [@picelli2014fulllength] with ERCC spike-ins.
Counts were obtained from NCBI GEO as a supplementary file using the accession number [GSE61533](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61533).


```r
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask=FALSE)
wilson.fname <- bfcrpath(bfc, file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
    "GSE61nnn/GSE61533/suppl/GSE61533_HTSEQ_count_results.xls.gz"))

library(R.utils)
wilson.name2 <- "GSE61533_HTSEQ_count_results.xls"
gunzip(wilson.fname, destname=wilson.name2, remove=FALSE, overwrite=TRUE)
```

Our first task is to load the count matrix into memory.
In this case, some work is required to retrieve the data from the Gzip-compressed Excel format.


```r
library(gdata)
all.counts <- read.xls(wilson.name2, sheet=1, header=TRUE)
rownames(all.counts) <- all.counts$ID
all.counts <- as.matrix(all.counts[,-1])
```

We store the results in a `SingleCellExperiment` object and identify the rows corresponding to the spike-ins based on the row names.


```r
library(SingleCellExperiment)
sce.hsc <- SingleCellExperiment(list(counts=all.counts))
dim(sce.hsc)
```

```
## [1] 38498    96
```

```r
is.spike <- grepl("^ERCC", rownames(sce.hsc))
isSpike(sce.hsc, "ERCC") <- is.spike
summary(is.spike)
```

```
##    Mode   FALSE    TRUE 
## logical   38406      92
```

### Quality control and normalization

For each cell, we calculate quality control metrics using the `calculateQCMetrics` function from *[scater](https://bioconductor.org/packages/3.9/scater)* [@mccarthy2017scater] as previously described.
We filter out HSCs that are outliers for any metric, under the assumption that these represent low-quality libraries. 


```r
library(scater)
sce.hsc <- calculateQCMetrics(sce.hsc)
libsize.drop <- isOutlier(sce.hsc$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce.hsc$total_features_by_counts, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce.hsc$pct_counts_ERCC, nmads=3, type="higher")
sce.hsc <- sce.hsc[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    BySpike=sum(spike.drop), Remaining=ncol(sce.hsc))
```

```
##   ByLibSize ByFeature BySpike Remaining
## 1         2         2       3        92
```

We remove genes that are not expressed in any cell to reduce computational work in downstream steps. 


```r
to.keep <- nexprs(sce.hsc, byrow=TRUE) > 0
sce.hsc <- sce.hsc[to.keep,]
summary(to.keep)
```

```
##    Mode   FALSE    TRUE 
## logical   17229   21269
```

We apply the deconvolution method [@lun2016pooling] to compute size factors for the endogenous genes [@lun2016pooling].
Separate size factors for the spike-in transcripts are also calculated, as previously discussed.
We then calculate log-transformed normalized expression values for further use.


```r
library(scran)
sce.hsc <- simpleSumFactors(sce.hsc)
summary(sizeFactors(sce.hsc))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4238  0.8177  0.9363  1.0000  1.1804  2.0151
```

```r
sce.hsc <- computeSpikeFactors(sce.hsc, type="ERCC", general.use=FALSE)
summary(sizeFactors(sce.hsc, "ERCC"))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2562  0.6198  0.8623  1.0000  1.2122  3.0289
```

```r
sce.hsc <- normalize(sce.hsc)
```

## Testing for significantly positive biological components

We fit a mean-variance trend to the spike-in transcripts to quantify the technical component of the variance, as previously described.
The biological component for each gene is defined as the difference between its total variance and the fitted value of the trend (Figure \@ref(fig:hvgplothsc)).


```r
var.fit <- trendVar(sce.hsc, parametric=TRUE, loess.args=list(span=0.3))
var.out <- decomposeVar(sce.hsc, var.fit)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce.hsc)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3-var_files/figure-html/hvgplothsc-1.png" alt="Variance of normalized log-expression values for each gene in the HSC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)." width="100%" />
<p class="caption">(\#fig:hvgplothsc)Variance of normalized log-expression values for each gene in the HSC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red).</p>
</div>

We define HVGs as those genes that have a biological component that is significantly greater than zero.
We use a false discovery rate (FDR) of 5% after correcting for multiple testing with the Benjamini-Hochberg method.


```r
hvg.out <- var.out[which(var.out$FDR <= 0.05),]
nrow(hvg.out)
```

```
## [1] 504
```

We rank the results to focus on genes with larger biological components.
This highlights an interesting aspect of the underlying hypothesis test, which is based on the ratio of the total variance to the expected technical variance.
Ranking based on _p_-value tends to prioritize HVGs that are more likely to be true positives but, at the same time, less likely to be interesting.
This is because the ratio can be very large for HVGs that have very low total variance and do not contribute much to the cell-cell heterogeneity.


```r
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
write.table(file="hsc_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
head(hvg.out)
```

```
## DataFrame with 6 rows and 6 columns
##                      mean            total              bio             tech
##                 <numeric>        <numeric>        <numeric>        <numeric>
## Fos      6.46196799882072 19.5610328073659 12.8045771643956 6.75645564297032
## Dusp1    6.82477646028322 15.5991684453524 10.0847547346845 5.51441371066791
## Rgs1     5.30721777600443 20.2441159375576 9.92804397134609 10.3160719662115
## Ppp1r15a 6.66598045859203 14.5028371912445 8.45276980013916 6.05006739110531
## Ly6a     8.39740181201201 10.0467504816539 8.03818961156809 2.00856087008578
## Egr1     6.71502634373587 13.8389703181346 7.95641515628638 5.88255516184822
##                       p.value                  FDR
##                     <numeric>            <numeric>
## Fos      1.10398969403887e-18 7.56093753689073e-16
## Dusp1    8.15129503928221e-18 4.80269775217264e-15
## Rgs1     1.18496143793071e-07 1.35612234428254e-05
## Ppp1r15a 1.91754808050171e-12 5.49322855056976e-10
## Ly6a     7.87714321198765e-50 2.38688692384957e-46
## Egr1     6.40390775621137e-12 1.65650350508536e-09
```

We check the distribution of expression values for the genes with the largest biological components.
We see that the variance estimate is not driven by one or two outlier cells (Figure \@ref(fig:hvgvioplothsc)).


```r
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotExpression(sce.hsc, features=rownames(hvg.out)[1:10]) + fontsize
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3-var_files/figure-html/hvgvioplothsc-1.png" alt="Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the HSC dataset. Each point represents the log-expression value in a single cell." width="100%" />
<p class="caption">(\#fig:hvgvioplothsc)Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the HSC dataset. Each point represents the log-expression value in a single cell.</p>
</div>

## Further comments

There are many other strategies for defining HVGs, based on a variety of metrics:

- the coefficient of variation, using the `technicalCV2()` function [@brennecke2013accounting] or the `DM()` function [@kim2015characterizing] in *[scran](https://bioconductor.org/packages/3.9/scran)*.
- the dispersion parameter in the negative binomial distribution, using the `estimateDisp()` function in *[edgeR](https://bioconductor.org/packages/3.9/edgeR)* [@mccarthy2012differential].
- a proportion of total variability, using methods in the *[BASiCS](https://bioconductor.org/packages/3.9/BASiCS)* package [@vallejos2015basics].

Here, we used the variance of the log-expression values because the log-transformation protects against genes with strong expression in only one or two cells.
This reduces the risk that the set of top HVGs is not dominated by genes with (mostly uninteresting) outlier expression patterns.

We also save the HSC dataset to file for later use, using `saveRDS()` as previously described.


```r
saveRDS(sce.hsc, file="hsc_data.rds")
```

# Fine-tuning the mean-variance trend fit

## Details of trend fitting parameters

Fitting of the mean-variance trend with `trendVar()` is complicated by the small number of spike-in transcripts and the uneven distribution of their abundances.
For low numbers of cells, these issues are exacerbated by the low precision of the variance estimates.
Thus, some tuning of trend parameters such as `span` may be required to achieve a suitable fit - see `?trendVar` for more details.
Setting `parametric=TRUE` is especially useful for modelling the expected wave-like shape of the mean-variance relationship.
(This is not the default setting as it is not robust for arbitrary trend shapes.)

The `trendVar` function will also automatically filter out low-abundance genes prior to trend fitting.
This ensures that low-abundance genes do not interfere with the fit due to discreteness, which biases the estimate of variability of the variances around the trend;
or due to the frequency of low-abundance genes, which reduces the sensitivity of span-based smoothing algorithms at higher abundances.
The internal choice of filtering strategy involves a number of considerations:

- Filtering uses the average of log-expression values rather than the (library size-adjusted) average count.
The mean log-expression is independent of the variance estimate in a linear modelling framework [@bourgon2010independent], 
which ensures that the filter does not introduce spurious trends in the variances at the filter boundary.
- The filter threshold is specified with the `min.mean` argument in `trendVar`.
We use the default threshold of 0.1 (`min.mean`) based on the appearance of discrete patterns in the variance estimates for simulated Poisson-distributed counts.
Lower thresholds of 0.001-0.01 may be more suitable for very sparse data, e.g., from droplet-based protocols.
- The filter used in `trendVar` is _not_ applied in `decomposeVar` by default.
Retention of all genes ensures that weak biological signal from rare subpopulations is not discarded.
To apply the filter in `decomposeVar`, users should set `subset.row=rowMeans(logcounts(sce)) > 0.1` in the function call.

**Comments from Aaron:**

- On occasion, users may observe a warning from `trendVar()` about the lack of centering in the size factors.
Recall that the trend is fitted to the mean and variances of the spike-in transcripts, and the technical component for each _endogenous_ gene is estimated by interpolation.
This assumes that an endogenous gene is comparable to a spike-in transcript of the same abundance.
In particular, we assume that variation is primarily driven by the magnitude of the counts, based on the well-known mean-variance relationships in count models.
Thus, we need to ensure that similarities in the average counts are preserved in the normalized expression values.
This is achieved by centering the gene- and spike-in-based size factors in `normalize()`, such that features with similar average counts will also have similar normalized abundances.
However, if the `SingleCellExperiment` object was manipulated (e.g., subsetted) _after_ `normalize()` and _before_ `trendVar()`, centering may not be preserved - hence the warning. 

## When spike-ins are unavailable

In some datasets, spike-in RNA may not have been added in appropriate quantities (or indeed at all).
It may also be inappropriate to assume Poisson technical noise with `makeTechTrend()`, especially for read count data where amplification noise is non-negligible.
In such cases, an alternative approach is to fit the trend to the variance estimates of the endogenous genes.
This is done using the `use.spikes=FALSE` setting in `trendVar`, as shown below for the HSC dataset.


```r
var.fit.nospike <- trendVar(sce.hsc, parametric=TRUE, 
    use.spikes=FALSE, loess.args=list(span=0.2))
var.out.nospike <- decomposeVar(sce.hsc, var.fit.nospike)
```

The simplest interpretation of the results assumes that the majority of genes are not variably expressed.
This means that the technical component dominates the total variance for most genes, such that the fitted trend can be treated as an estimate of the technical component.
In Figure \@ref(fig:hvgplot416b2), the trend passes through or close to most of the spike-in variances, indicating that our assumption is valid.


```r
plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
    xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3-var_files/figure-html/hvgplot416b2-1.png" alt="Variance of normalized log-expression values for each gene in the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the endogenous genes (black), with spike-in transcripts shown in red." width="100%" />
<p class="caption">(\#fig:hvgplot416b2)Variance of normalized log-expression values for each gene in the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the endogenous genes (black), with spike-in transcripts shown in red.</p>
</div>

If our assumption does not hold, the output of decomposeVar is more difficult to interpret. 
The fitted value of the trend can no longer be generally interpreted as the technical component, as it contains some biological variation as well.
Instead, recall that the biological component reported by `decomposeVar` represents the residual for each gene over the majority of genes with the same abundance.
One could assume that the variabilities of most genes are driven by constitutive "house-keeping" processes, which are biological in origin but generally uninteresting.
Any gene with an increase in its variance is _relatively_ highly variable and can be prioritized for further study.

# Blocking on uninteresting factors of variation

## Using the `block=` argument

Our previous analysis of the 416B dataset specified `block=` in `trendVar()` to ensure that systematic differences between plates do not inflate the variance.
This involves estimating the mean and variance of the log-expression _separately_ in each plate,
followed by fitting a single trend to the plate-specific means and variances of all spike-in transcripts.
In doing so, we implicitly assume that the trend is the same between plates, which is reasonable for this dataset (Figure \@ref(fig:trendplotblock-416b)).


```r
# Loading the saved object.
sce.416B <- readRDS("416B_data.rds") 

# Repeating the trendVar() call.
var.fit <- trendVar(sce.416B, parametric=TRUE, block=sce.416B$Plate,
    loess.args=list(span=0.3))

matplot(var.fit$means, var.fit$vars, col=c("darkorange", "forestgreen"),
    xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit$trend(x), add=TRUE, col="red")
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3-var_files/figure-html/trendplotblock-416b-1.png" alt="Plate-specific variance estimates for all spike-in transcripts in the 416B dataset, plotted against the plate-specific means. Each point represents a spike-in transcript, numbered by the plate from which the values were estimated. The red line denotes the fitted mean-variance trend." width="100%" />
<p class="caption">(\#fig:trendplotblock-416b)Plate-specific variance estimates for all spike-in transcripts in the 416B dataset, plotted against the plate-specific means. Each point represents a spike-in transcript, numbered by the plate from which the values were estimated. The red line denotes the fitted mean-variance trend.</p>
</div>

The use of `block=` also assumes that the average size factor within each plate is close to unity for both endogenous genes and spike-in transcripts.
This means that scaling normalization will preserve the magnitude of the counts, allowing genes of the same average abundance to be compared within and across plates.
Here, the distributions of size factors exhibit only modest deviations from unity in the averages for the endogenous genes and spike-in transcripts (Figure \@ref(fig:sizefacplot-416b)), indicating that our assumption is again reasonable.


```r
tmp.416B <- sce.416B
tmp.416B$log_size_factor <- log(sizeFactors(sce.416B))
tmp.416B$log_size_factor_ERCC <- log(sizeFactors(sce.416B, "ERCC"))
p1 <- plotColData(tmp.416B, x="Plate", y="log_size_factor")
p2 <- plotColData(tmp.416B, x="Plate", y="log_size_factor_ERCC")
multiplot(p1, p2, cols=2)
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3-var_files/figure-html/sizefacplot-416b-1.png" alt="Plate-specific distribution of the size factors for endogenous genes (left) and spike-in transcripts (right)." width="960" />
<p class="caption">(\#fig:sizefacplot-416b)Plate-specific distribution of the size factors for endogenous genes (left) and spike-in transcripts (right).</p>
</div>

The most obvious violation of the above assumption occurs when more spike-in RNA is added in a particular batch.
Spike-in size factors would then be systematically larger than unity in that batch, meaning that average abundances after normalization would not be comparable between spike-in transcripts and endogenous genes.
This would compromise the accuracy of estimation of the technical component by interpolation from the trend.
More generally, batches that are processed separately may not exhibit to the same technical variation, such that the use of a single trend would be inappropriate.

## Fitting batch-specific trends

For datasets containing multiple batches, an alternative strategy is to perform trend fitting and variance decomposition separately for each batch.
This accommodates differences in the mean-variance trends between batches, especially if a different amount of spike-in RNA was added to the cells in each batch.
We demonstrate this approach by treating each plate in the 416B dataset as a different batch, using the `multiBlockVar()` function.
This yields plate-specific estimates of the biological and technical components for each gene.


```r
sce.416B.2 <- multiBlockNorm(sce.416B, sce.416B$Plate)
comb.out <- multiBlockVar(sce.416B.2, block=sce.416B.2$Plate,
    trend.args=list(parametric=TRUE, loess.args=list(span=0.4)))
```

Statistics are combined across multiple batches using the `combineVar()` function within `multiBlockVar()`.
This function computes a weighted average across batches for the means and variances, and applies Fisher's method for combining the _p_-values.
These results can be used in downstream functions such as `denoisePCA`, or for detecting highly variable genes (see below).


```r
head(comb.out[,1:6])
```

```
## DataFrame with 6 rows and 6 columns
##                                   mean               total
##                              <numeric>           <numeric>
## ENSMUSG00000103377  0.0077852662424069  0.0110910243874505
## ENSMUSG00000103147  0.0342704714401078  0.0710324324735765
## ENSMUSG00000103161 0.00533642721877204 0.00521168730727283
## ENSMUSG00000102331  0.0186505112011929  0.0333123382277785
## ENSMUSG00000102948  0.0588935602840672  0.0879186539054153
## Rp1                 0.0968879289012201   0.452061385087204
##                                    bio               tech           p.value
##                              <numeric>          <numeric>         <numeric>
## ENSMUSG00000103377 -0.0272027560050436 0.0382937803924941                 1
## ENSMUSG00000103147 -0.0757714743788675  0.146803906852444 0.999999999968904
## ENSMUSG00000103161 -0.0176478876202435 0.0228595749275163                 1
## ENSMUSG00000102331 -0.0537822238874146  0.087094562115193 0.999999999999963
## ENSMUSG00000102948  -0.169169236387022  0.257087890292437                 1
## Rp1                0.00876187354127483  0.443299511545929 0.175559732314035
##                                  FDR
##                            <numeric>
## ENSMUSG00000103377                 1
## ENSMUSG00000103147                 1
## ENSMUSG00000103161                 1
## ENSMUSG00000102331                 1
## ENSMUSG00000102948                 1
## Rp1                0.569934085943856
```

We visualize the quality of the batch-specific trend fits by extracting the relevant statistics from `comb.out` (Figure \@ref(fig:hvgplotbatch416b)).


```r
par(mfrow=c(1,2))
is.spike <- isSpike(sce.416B.2)
for (plate in levels(sce.416B.2$Plate)) {
    cur.out <- comb.out$per.block[[plate]]
    plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
        ylab="Variance of log-expression", main=plate)
    curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
    points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
}
```

<div class="figure">
<img src="/home/cri.camres.org/lun01/AaronDocs/Research/simpleSingleCell/results/xtra-3-var_files/figure-html/hvgplotbatch416b-1.png" alt="Variance of normalized log-expression values for each gene in each plate of the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)." width="960" />
<p class="caption">(\#fig:hvgplotbatch416b)Variance of normalized log-expression values for each gene in each plate of the 416B dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red).</p>
</div>

By fitting separate trends, we avoid the need to assume that a single trend is present across batches.
However, this also reduces the precision of each trend fit, as less information is available within each batch.
We recommend using `block=` as the default unless there is clear evidence for differences in the trends between batches.

**Comments from Aaron:**

- We run `multiBlockNorm()` to adjust the size factors within each level of the blocking factor.
Specifically, the spike-in size factors across cells in a given batch is scaled so that the mean is equal to that of the gene-based size factors for the same set of cells.
Log-normalized expression values are then recalculated using these centred size factors.
This procedure ensures that the average abundances of the spike-in transcripts are comparable to the endogenous genes,
avoiding problems due to differences in the quantity of spike-in RNA between batches.
Otherwise, if the globally-centred size factors (from `normalize()`) were used, 
there would be a systematic difference in the scaling of spike-in transcripts compared to endogenous genes in batches with more or less spike-in RNA than the dataset average.
The fitted trend would then be shifted along the x-axis and fail to accurately capture the technical component for each gene.

## Using the `design=` argument

For completeness, it is worth mentioning the `design=` argument in `trendVar()`.
This will estimate the residual variance from a linear model fitted to the log-normalized expression values for each gene.
The linear model can include blocking factors for known unwanted factors of variation, ensuring that they do not inflate the variance estimate.
The technical component for each gene is obtained at the average abundance across all cells.


```r
lfit <- trendVar(sce.416B, design=model.matrix(~sce.416B$Plate))
```

We do not recommend using this approach for categorical blocking factors in one-way layouts.
This is because it does not consider the mean of each blocking level, resulting in an inaccurate estimate of the technical component in the presence of a strong blocking effect. 
However, it is the only choice for dealing with real covariates or multiple blocking factors in an additive model.

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
##  [1] gdata_2.18.0                R.utils_2.7.0              
##  [3] R.oo_1.22.0                 R.methodsS3_1.7.1          
##  [5] scran_1.11.1                scater_1.11.2              
##  [7] ggplot2_3.1.0               SingleCellExperiment_1.5.0 
##  [9] SummarizedExperiment_1.13.0 DelayedArray_0.9.0         
## [11] BiocParallel_1.17.1         matrixStats_0.54.0         
## [13] Biobase_2.43.0              GenomicRanges_1.35.0       
## [15] GenomeInfoDb_1.19.0         IRanges_2.17.1             
## [17] S4Vectors_0.21.1            BiocGenerics_0.29.1        
## [19] bindrcpp_0.2.2              BiocFileCache_1.7.0        
## [21] dbplyr_1.2.2                knitr_1.20                 
## [23] BiocStyle_2.11.0           
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7             
##  [3] httr_1.3.1               rprojroot_1.3-2         
##  [5] dynamicTreeCut_1.63-1    tools_3.6.0             
##  [7] backports_1.1.2          R6_2.3.0                
##  [9] HDF5Array_1.11.0         vipor_0.4.5             
## [11] DBI_1.0.0                lazyeval_0.2.1          
## [13] colorspace_1.3-2         withr_2.1.2             
## [15] tidyselect_0.2.5         gridExtra_2.3           
## [17] bit_1.1-14               curl_3.2                
## [19] compiler_3.6.0           BiocNeighbors_1.1.1     
## [21] labeling_0.3             bookdown_0.7            
## [23] scales_1.0.0             rappdirs_0.3.1          
## [25] stringr_1.3.1            digest_0.6.18           
## [27] rmarkdown_1.10           XVector_0.23.0          
## [29] pkgconfig_2.0.2          htmltools_0.3.6         
## [31] limma_3.39.1             highr_0.7               
## [33] rlang_0.3.0.1            RSQLite_2.1.1           
## [35] DelayedMatrixStats_1.5.0 bindr_0.1.1             
## [37] gtools_3.8.1             dplyr_0.7.8             
## [39] RCurl_1.95-4.11          magrittr_1.5            
## [41] GenomeInfoDbData_1.2.0   Matrix_1.2-15           
## [43] Rcpp_1.0.0               ggbeeswarm_0.6.0        
## [45] munsell_0.5.0            Rhdf5lib_1.5.0          
## [47] viridis_0.5.1            stringi_1.2.4           
## [49] yaml_2.2.0               edgeR_3.25.0            
## [51] zlibbioc_1.29.0          rhdf5_2.27.0            
## [53] plyr_1.8.4               grid_3.6.0              
## [55] blob_1.1.1               crayon_1.3.4            
## [57] lattice_0.20-38          cowplot_0.9.3           
## [59] locfit_1.5-9.1           pillar_1.3.0            
## [61] igraph_1.2.2             reshape2_1.4.3          
## [63] glue_1.3.0               evaluate_0.12           
## [65] BiocManager_1.30.3       gtable_0.2.0            
## [67] purrr_0.2.5              assertthat_0.2.0        
## [69] xfun_0.4                 viridisLite_0.3.0       
## [71] tibble_1.4.2             beeswarm_0.2.3          
## [73] memoise_1.1.0            statmod_1.4.30
```

# References

