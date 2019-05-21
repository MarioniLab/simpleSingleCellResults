---
title: Multi-step correction of batch effects in single-cell RNA-seq data
author: 
- name: Aaron T. L. Lun
  affiliation: Cancer Research UK Cambridge Institute, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, United Kingdom
- name: Michael D. Morgan
  affiliation: Wellcome Trust Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridge CB10 1SA, United Kingdom
date: "2019-05-20"
vignette: >
  %\VignetteIndexEntry{11. Advanced batch correction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}    
output: 
  BiocStyle::html_document:
    titlecaps: false
    toc_float: true
bibliography: ref.bib
---





# Introduction

In the [previous batch correction workflow](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/batch.html), we merged two scRNA-seq data sets involving human pancreas cell populations.
This workflow extends the previous example to describe how to perform a more complicated merge operation involving different levels of batch effects.
Here, we will use pancreas data sets generated using Smart-based technologies [@segerstolpe2016singlecell;@lawlor2017singlecell],
and merge them with the previous data sets generated using CEL-seq-based methods [@grun2016denovo;@muraro2016singlecell].
Our overall strategy is to use a hierarchical merge to remove batch effects within each technology,
followed by the removal of batch effects between technologies.

# Loading in the data

## SMARTer, GSE86469

### Reading in data

Here, we use data from the @lawlor2017singlecell study of pancreatic islet cells from healthy and type II diabetic donors.
This was generated using the SMARTer protocol on the Fluidigm C1 system.
We download and cache the count matrix using the *[BiocFileCache](https://bioconductor.org/packages/3.10/BiocFileCache)* package.


```r
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)    
count.tab <- bfcrpath(bfc, file.path(
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86469/suppl",
    "GSE86469_GEO.islet.single.cell.processed.data.RSEM.raw.expected.counts.csv.gz"
))
```

We read in the count matrix in sparse format using *[scater](https://bioconductor.org/packages/3.10/scater)* [@mccarthy2017scater].


```r
library(scater)
mat <- readSparseCounts(count.tab, sep=",", quote='"', row.names=1)
dim(mat)
```

```
## [1] 26616   638
```

We load in the metadata from NCBI GEO using the *[GEOquery](https://bioconductor.org/packages/3.10/GEOquery)* package [@davis2007geoquery].


```r
library(GEOquery)
metadata <- pData(getGEO("GSE86469")[[1]])
metadata <- metadata[,c("title", "cell type:ch1", "islet unos id:ch1")]

rownames(metadata) <- metadata$title
metadata <- metadata[,-1]
colnames(metadata) <- c("CellType", "Donor")

stopifnot(identical(colnames(mat), rownames(metadata)))
head(metadata)
```

```
##                 CellType   Donor
## 10th_C10_S104 None/Other ACIW009
## 10th_C11_S96        Beta ACIW009
## 10th_C13_S61        Beta ACIW009
## 10th_C14_S53        Beta ACIW009
## 10th_C16_S105 None/Other ACIW009
## 10th_C17_S97        Beta ACIW009
```

Finally, we create a `SingleCellExperiment` object.


```r
library(SingleCellExperiment)
sce.gse86469 <- SingleCellExperiment(list(counts=mat), colData=metadata)
isSpike(sce.gse86469, "ERCC") <- grep("^ERCC-", rownames(sce.gse86469))
sce.gse86469
```

```
## class: SingleCellExperiment 
## dim: 26616 638 
## metadata(0):
## assays(1): counts
## rownames(26616): ENSG00000229483 ENSG00000232849 ... ENSG00000251576
##   ENSG00000082898
## rowData names(0):
## colnames(638): 10th_C10_S104 10th_C11_S96 ... 9th-C96_S81 9th-C9_S13
## colData names(2): CellType Donor
## reducedDimNames(0):
## spikeNames(1): ERCC
```

### Quality control and normalization

We remove low-quality cells based on outliers for various quality control metrics,
such as the total library size and the number of expressed genes.
This is similar to what was described [previously](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/reads.html#3_quality_control_on_the_cells).
Note that this data does not contain any counts for spike-in transcripts, 
so the spike-in percentage is not used here.


```r
sce.gse86469 <- calculateQCMetrics(sce.gse86469, compact=TRUE)
qc.mat <- cbind(
    NFeatures=isOutlier(sce.gse86469$scater_qc$all$total_features_by_counts, 
        log=TRUE, type="lower", nmads=3),
    LibSize=isOutlier(sce.gse86469$scater_qc$all$total_counts, 
        log=TRUE, type="lower", nmads=3)
)
colSums(qc.mat)
```

```
## NFeatures   LibSize 
##         0         3
```

```r
discard <- rowMeans(qc.mat) > 0
sce.gse86469 <- sce.gse86469[,!discard]
summary(discard)
```

```
##    Mode   FALSE    TRUE 
## logical     635       3
```

We compute size factors with the deconvolution method from the *[scran](https://bioconductor.org/packages/3.10/scran)* package [@lun2016pooling].
Pre-clustering is performed using `quickCluster()` to avoid pooling together very different cells.
Note the use of `IrlbaParam()` from *[BiocSingular](https://bioconductor.org/packages/3.10/BiocSingular)* to speed up the PCA calculations.


```r
library(scran)
library(BiocSingular)
clusters <- quickCluster(sce.gse86469, BSPARAM=IrlbaParam())
table(clusters)
```

```
## clusters
##   1   2   3 
## 242 270 123
```

```r
sce.gse86469 <- computeSumFactors(sce.gse86469, clusters=clusters)
summary(sizeFactors(sce.gse86469))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2912  0.7694  0.9596  1.0000  1.2007  2.9185
```

There is no need for spike-in normalization, as there are actually no spike-in counts. 
We thus proceed directly to calculation of the log-transformed normalized expression values for endogenous genes.


```r
# Ignore warnings due to no spike-in size factors.
suppressWarnings(sce.gse86469 <- normalize(sce.gse86469))
```

### Modelling variability
    
Given that no spike-ins are available, we fit a mean-dependent trend to the variances of the endogenous genes to model technical noise.
This requires the assumptions that have been stated [elsewhere](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/var.html#when-spike-ins-are-unavailable).


```r
fit <- trendVar(sce.gse86469, use.spikes=FALSE, 
    block=sce.gse86469$Donor, loess.args=list(span=0.05))
dec.gse86469 <- decomposeVar(sce.gse86469, fit)
dec.gse86469[order(dec.gse86469$bio, decreasing=TRUE),]
```

```
## DataFrame with 26616 rows and 6 columns
##                             mean             total               bio
##                        <numeric>         <numeric>         <numeric>
## ENSG00000115263 10.1536955344348  43.9282794110398  40.7827047057491
## ENSG00000254647 8.11583744343173  38.6406987550627  33.8414096768278
## ENSG00000169903  5.7360754103119  29.0826673417615  22.6937447312962
## ENSG00000118271 12.0056875945825  20.6222933766225  18.6428162898466
## ENSG00000197249 4.12395232568971  26.0384575981412   18.574407947614
## ...                          ...               ...               ...
## ENSG00000214105 3.57156880461368  1.47686769949045 -6.26100443295033
## ENSG00000144218 3.74525076548908  1.38625572839096   -6.295430365288
## ENSG00000261177 4.20932168304822  1.24423250040548 -6.35184676330201
## ENSG00000167995  4.4658056980619 0.975647920550666 -6.49349680189003
## ENSG00000141150 4.04259310145711  1.12290389575419 -6.54082049078115
##                             tech               p.value                   FDR
##                        <numeric>             <numeric>             <numeric>
## ENSG00000115263 3.14557470529066                     0                     0
## ENSG00000254647 4.79928907823488                     0                     0
## ENSG00000169903 6.38892261046531 2.41407654662433e-308 1.28506122729906e-304
## ENSG00000118271 1.97947708677594                     0                     0
## ENSG00000197249  7.4640496505272 5.76405676392468e-174  1.1801241140663e-170
## ...                          ...                   ...                   ...
## ENSG00000214105 7.73787213244079                     1                     1
## ENSG00000144218 7.68168609367896                     1                     1
## ENSG00000261177  7.5960792637075                     1                     1
## ENSG00000167995 7.46914472244069                     1                     1
## ENSG00000141150 7.66372438653534                     1                     1
```

Figure \@ref(fig:var-gse86469) shows the strong mean-variance relationship that is typical of read count data.


```r
plot(fit$mean, fit$var, xlab="Mean log-expression",
    ylab="Variance of log-expression")
curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
```

<div class="figure">
<img src="multibatch_files/figure-html/var-gse86469-1.png" alt="Variance of normalized log-expression values for each gene in the GSE86469 data set, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of all genes." width="100%" />
<p class="caption">(\#fig:var-gse86469)Variance of normalized log-expression values for each gene in the GSE86469 data set, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of all genes.</p>
</div>

## Smart-seq2, E-MTAB-5061

### Reading in the data

Here, we use data from the @@segerstolpe2016singlecell study.
The good news is that the authors have provided a count table in the ArrayExpress entry for this project.
We download it using *[BiocFileCache](https://bioconductor.org/packages/3.10/BiocFileCache)* to cache the results:


```r
bfc <- BiocFileCache(ask=FALSE)    
emat <- bfcrpath(bfc, file.path("https://www.ebi.ac.uk/arrayexpress",
    "experiments/E-MTAB-5061/files/E-MTAB-5061.processed.1.zip"))
count.file <- "pancreas_refseq_rpkms_counts_3514sc.txt"
```

The bad news is that the count table is needlessly complex:

- The first 2 columns contain the gene symbol and NCBI GenBank transcript identifiers for each row.
- The next `X` columns are the RPKMs, for `X` cells.
- The remaining `X` columns are the counts.

This requires some additional work to extract the useful data.
The first line contains the names of the cells, 
so we can use this to determine the number and indices of the columns with per-cell counts.


```r
col.names <- read.table(unz(emat, count.file), header=FALSE, sep="\t", 
    stringsAsFactors=FALSE, comment.char="", nrows = 1)[,-1]
ncells <- length(col.names)

what <- vector("list", ncells*2 + 2)
what[[1]] <- "character"
what[seq_len(ncells) + ncells + 2] <- "integer"
```

We then read in the gene symbols and the counts.
We use the gene symbols as the GenBank IDs have been rather clumsily formatted.


```r
emtab.df <- read.table(unz(emat, count.file), header=FALSE, sep="\t", 
    stringsAsFactors=FALSE, colClasses=what, skip=1)
row.names <- emtab.df[,1]
emtab.df <- emtab.df[,-1]
colnames(emtab.df) <- col.names
dim(emtab.df)
```

```
## [1] 26271  3514
```

Some work is required to translate the symbols into Ensembl gene identifiers via the *[org.Hs.eg.db](https://bioconductor.org/packages/3.10/org.Hs.eg.db)* package.


```r
library(org.Hs.eg.db)
ens.id <- mapIds(org.Hs.eg.db, keys=row.names, 
    keytype="SYMBOL", column="ENSEMBL")
ens.id <- ifelse(is.na(ens.id), row.names, ens.id)

keep <- !duplicated(ens.id)
emtab.df <- emtab.df[keep,]
rownames(emtab.df) <- ens.id[keep]
head(rownames(emtab.df))
```

```
## [1] "ENSG00000118473" "ENSG00000142920" "ENSG00000169504" "ENSG00000186094"
## [5] "ENSG00000157191" "ENSG00000162426"
```

We read in the metadata and extract the appropriate columns.


```r
meta.fname <- bfcrpath(bfc, file.path("https://www.ebi.ac.uk/arrayexpress",
    "files/E-MTAB-5061/E-MTAB-5061.sdrf.txt"))
emtab.sdrf <- read.delim(meta.fname, stringsAsFactors=FALSE)
stopifnot(identical(sort(emtab.sdrf$Source.Name), sort(colnames(emtab.df))))    

emtab.sdrf <- emtab.sdrf[match(colnames(emtab.df), emtab.sdrf$Source.Name),]
emtab.meta <- emtab.sdrf[, c("Assay.Name", 
    "Characteristics.cell.type.", "Characteristics.individual.",
    "Characteristics.single.cell.well.quality.")]
colnames(emtab.meta) <- c("Sample", "CellType", "Donor", "Quality")
emtab.meta$Study <- "E-MTAB-5061"
head(emtab.meta)
```

```
##            Sample       CellType     Donor          Quality       Study
## 406 HP1502401_N13 not applicable HP1502401 low quality cell E-MTAB-5061
## 172 HP1502401_D14 not applicable HP1502401 low quality cell E-MTAB-5061
## 219 HP1502401_F14 not applicable HP1502401 low quality cell E-MTAB-5061
## 312 HP1502401_J13 not applicable HP1502401 low quality cell E-MTAB-5061
## 124 HP1502401_B13 not applicable HP1502401 low quality cell E-MTAB-5061
## 265 HP1502401_H13     gamma cell HP1502401               OK E-MTAB-5061
```

Some editing of the cell type labels is necessary for consistency with GSE86469.


```r
emtab.meta$CellType <- gsub(" cell", "", emtab.meta$CellType)
emtab.meta$CellType <- paste0(
    toupper(substr(emtab.meta$CellType, 1, 1)),
    substring(emtab.meta$CellType, 2))
table(emtab.meta$CellType)
```

```
## 
##                 Acinar                  Alpha                   Beta 
##                    185                    886                    270 
##          Co-expression                  Delta                 Ductal 
##                     39                    114                    386 
##            Endothelial                Epsilon                  Gamma 
##                     16                      7                    197 
##                   Mast           MHC class II         Not applicable 
##                      7                      5                   1305 
##                    PSC           Unclassified Unclassified endocrine 
##                     54                      2                     41
```

Finally, we create a `SingleCellExperiment` object.


```r
sce.emtab <- SingleCellExperiment(list(counts=as.matrix(emtab.df)), 
    colData=emtab.meta)
isSpike(sce.emtab, "ERCC") <- grep("^ERCC_", rownames(sce.emtab))
sce.emtab
```

```
## class: SingleCellExperiment 
## dim: 25471 3514 
## metadata(0):
## assays(1): counts
## rownames(25471): ENSG00000118473 ENSG00000142920 ...
##   ERCC_0.01430512:mix1_0.02861023:mix2 eGFP
## rowData names(0):
## colnames(3514): HP1502401_N13 HP1502401_D14 ... HP1526901T2D_O11
##   HP1526901T2D_A8
## colData names(5): Sample CellType Donor Quality Study
## reducedDimNames(0):
## spikeNames(1): ERCC
```

### Quality control and normalization

We first remove the low quality cells that were marked by the authors.


```r
low.qual <- sce.emtab$Quality == "low quality cell"
sce.emtab <- sce.emtab[,!low.qual]
summary(low.qual)
```

```
##    Mode   FALSE    TRUE 
## logical    2337    1177
```

We also remove low quality cells based on our own quality control metrics.
It is debatable whether these two separate rounds of quality control are necessary, 
but we do this for consistency with respect to the preprocessing performed across all data sets.


```r
sce.emtab <- calculateQCMetrics(sce.emtab, compact=TRUE)
qc.mat <- cbind(
    NFeatures=isOutlier(sce.emtab$scater_qc$all$total_features_by_counts, 
        log=TRUE, type="lower", nmads=3),
    LibSize=isOutlier(sce.emtab$scater_qc$all$total_counts, 
        log=TRUE, type="lower", nmads=3),
    SpikePct=isOutlier(sce.emtab$scater_qc$feature_control_ERCC$pct_counts, 
        type="higher", nmads=3)
)
colSums(qc.mat)
```

```
## NFeatures   LibSize  SpikePct 
##       175       106       341
```

```r
discard <- rowMeans(qc.mat) > 0
sce.emtab <- sce.emtab[,!discard]
summary(discard)
```

```
##    Mode   FALSE    TRUE 
## logical    1951     386
```

We compute size factors using the pre-clustering and deconvolution approach.


```r
clusters <- quickCluster(sce.emtab, BSPARAM=IrlbaParam())
table(clusters)
```

```
## clusters
##   1   2   3   4   5   6   7   8 
## 181 175 253 202 344 272 414 110
```

```r
sce.emtab <- computeSumFactors(sce.emtab, clusters=clusters)
summary(sizeFactors(sce.emtab))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##  0.02297  0.39389  0.69609  1.00000  1.30830 11.70117
```

We also compute separate size factors for the spike-in counts.
Note that some cells have no spike-in counts and will not be useful for downstream steps that rely on spike-ins.


```r
sce.emtab <- computeSpikeFactors(sce.emtab, general.use=FALSE)
summary(sizeFactors(sce.emtab, "ERCC"))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.2279  0.4872  1.0000  1.2211 12.3395
```

### Modelling variability

Variance modelling requires some care in this data set.
The mean-variance trend is highly variable across cell populations from different donors,
presumably because each donor was processed in a separate batch.
Thus, we have to block on `Donor` using the `multiBlockVar()` function as discussed [here](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/var.html#fitting-batch-specific-trends).
We also have to remove cells with no spike-ins, as they are not useful for modelling technical noise;
and remove cells from donor `AZ`, which has very low spike-in concentrations.


```r
for.hvg <- sce.emtab[,sizeFactors(sce.emtab, "ERCC") > 0
    & sce.emtab$Donor!="AZ"]
for.hvg <- multiBlockNorm(for.hvg, for.hvg$Donor) 
dec.emtab <- multiBlockVar(for.hvg, for.hvg$Donor)
head(dec.emtab[order(dec.emtab$bio,decreasing=TRUE),-7])
```

```
## DataFrame with 6 rows and 6 columns
##                             mean            total              bio
##                        <numeric>        <numeric>        <numeric>
## ENSG00000115263 9.88190442484548 25.8969975556933 25.5852459231553
## ENSG00000118271 10.3307233103184 20.2742000137247 20.1086878178325
## ENSG00000089199 8.78634125225357 18.6524078339323 18.2663604764259
## ENSG00000166922 7.84098626894756 16.4774012002994 15.8457158517458
## ENSG00000197249 6.17277037945504 17.4619843294184 15.5437275596294
## ENSG00000100604 6.93096513654717 16.6103476582257 15.4730842154729
##                              tech   p.value       FDR
##                         <numeric> <numeric> <numeric>
## ENSG00000115263 0.311751632538025         0         0
## ENSG00000118271  0.16551219589229         0         0
## ENSG00000089199 0.386047357506408         0         0
## ENSG00000166922 0.631685348553599         0         0
## ENSG00000197249  1.91825676978892         0         0
## ENSG00000100604  1.13726344275278         0         0
```

Figure \@ref(fig:var-emtab) demonstrates the diversity of mean-variance relationships across different donors.


```r
all.donors <- unique(for.hvg$Donor)
par(mfrow=c(ceiling(length(all.donors)/3),3))
is.spike <- isSpike(for.hvg)
for (plate in all.donors) {
    cur.out <- dec.emtab$per.block[[plate]]
    plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, main=plate, 
        xlab="Mean log-expression", ylab="Variance of log-expression")
    curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
    points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
}
```

<div class="figure">
<img src="multibatch_files/figure-html/var-emtab-1.png" alt="Variance of normalized log-expression values for each gene in each donor of the EMTAB-5061 data set, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of spike-in transcripts within each donor." width="100%"  class="widefigure" />
<p class="caption">(\#fig:var-emtab)Variance of normalized log-expression values for each gene in each donor of the EMTAB-5061 data set, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of spike-in transcripts within each donor.</p>
</div>

# Feature selection across batches

Recall that our aim is to merge data from all four pancreas data sets [@segerstolpe2016singlecell;@lawlor2017singlecell;grun2016denovo;@muraro2016singlecell].
To do so, we first load in the CEL-seq data sets that we processed [previously](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/batch.html).


```r
sce.gse81076 <- readRDS("gse81076_sce.rds")
dec.gse81076 <- readRDS("gse81076_dec.rds") 
sce.gse85241 <- readRDS("gse85241_sce.rds") 
dec.gse85241 <- readRDS("gse85241_dec.rds") 
```

We define the universe of genes that are common across all batches.
This is made straightforward by the presence of common Ensembl identifiers.


```r
universe <- Reduce(intersect, list(rownames(dec.gse81076), 
    rownames(dec.gse85241), rownames(dec.gse86469), 
    rownames(dec.emtab)))
universe <- universe[!grepl("^ERCC-", universe)] # removing spike-ins.
length(universe)
```

```
## [1] 15478
```

We adjust the size factors with `multiBatchNorm()` to make them more comparable across batches.
This mitigates differences in scale and variance in the log-expression values between batches, especially between technologies.


```r
library(batchelor)
nout <- multiBatchNorm(sce.gse81076[universe,], sce.gse85241[universe,],
    sce.gse86469[universe,], sce.emtab[universe,])
sce.gse81076 <- nout[[1]]
sce.gse85241 <- nout[[2]]
sce.gse86469 <- nout[[3]]
sce.emtab <- nout[[4]]
```

We keep all genes with positive average biological components across all batches.
This is a relaxed approach to feature selection that ensures that interesting features in one or more batches are retained.


```r
mean.bio <- (dec.gse81076[universe,"bio"] + dec.gse85241[universe,"bio"]
    + dec.gse86469[universe,"bio"] + dec.emtab[universe,"bio"])/4
chosen <- universe[mean.bio > 0]
length(chosen)
```

```
## [1] 8697
```

We then subset all of the `SingleCellExperiment` objects so that only these features of interest are retained.


```r
sce.gse81076 <- sce.gse81076[chosen,]
sce.gse85241 <- sce.gse85241[chosen,]
sce.gse86469 <- sce.gse86469[chosen,]
sce.emtab <- sce.emtab[chosen,]
```

# Multi-batch principal components analysis

We use the `multiBatchPCA()` function to perform a PCA across _all_ batches to be merged.
This ensures that all cells are placed onto the same coordinate space, which would obviously not be possible if a PCA was performed for each batch separately.
Specifically, `multiBatchPCA()` performs a modified PCA to ensure that each supplied matrix contributes equally to the definition of the PC space.
This avoids problems with imbalances in the number of cells across batches, meaning that smaller batches (possibly with unique cell types) are not ignored.


```r
set.seed(1000)
pcs <- multiBatchPCA(
    gse81076=sce.gse81076,
    gse85241=sce.gse85241,
    gse86469=sce.gse86469,
    emtab=sce.emtab,
    BSPARAM=IrlbaParam(deferred=TRUE)
)
names(pcs)
```

```
## [1] "gse81076" "gse85241" "gse86469" "emtab"
```

Typical applications of `fastMNN()` will automatically call the `multiBatchPCA()` function on gene expression inputs.
However, this is not appropriate here as we will be performing a hierarchical merge.
Each call to `fastMNN()` will only involve a subset of batches,
and it would be difficult to try to merge results from two separate PCAs involving different subsets of data.
We need to run `multiBatchPCA()` manually on all batches to ensure that they are on the same coordinate system during merging.

**Comments from Aaron:**

- The `IrlbaParam(deferred=TRUE)` setting instructs `multiBatchPCA()` to perform a fast approximate PCA with methods from the *[irlba](https://CRAN.R-project.org/package=irlba)* package.
This involves some randomization and thus requires the seed to be set to obtain reproducible results.
- Here, we have applied `multiBatchPCA()` to the batch-level inputs for convenience.
It is also possible to supply donor-level matrices to equalize contributions across donors,
but this requires a bit more data manipulation that we will omit for the sake of simplicity.
- For full consistency with the `fastMNN()` defaults, we would call `cosineNorm()` on each log-expression matrix prior to running `multiBatchPCA()`. 
However, this is not technically necessary as all batches should be on the same scale already (see `?cosineNorm` for a discussion of this).

# Hierarchical merging

## Merging the Smart-based data

The @segerstolpe2016singlecell study contains strong donor effects that interfere with correction.
Several cell types exhibit strong per-donor effects such that the multiple clusters cannot be fully merged with the corresponding single cluster in GSE86469 (Figure \@ref(fig:smart-raw)).


```r
direct.smart <- fastMNN(pcs$gse86469, pcs$emtab, k=20, pc.input=TRUE)

set.seed(2000)
tsne.out <- Rtsne::Rtsne(direct.smart$corrected, perplexity=30, pca=FALSE)
df <- data.frame(x=tsne.out$Y[,1], y=tsne.out$Y[,2], 
    batch=rep(c("GSE86469", "E-MTAB-5061"), 
        c(ncol(sce.gse86469), ncol(sce.emtab))),
    donor=c(rep("unknown", ncol(sce.gse86469)), sce.emtab$Donor))

multiplot(
    ggplot(df) + geom_point(aes(x=x, y=y, color=donor)) +
        xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("By donor"),
    ggplot(df) + geom_point(aes(x=x, y=y, color=batch)) +
        xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("By batch"),
    cols=2)
```

<div class="figure">
<img src="multibatch_files/figure-html/smart-raw-1.png" alt="t-SNE plots of the merged Smart-based data sets. Each point is a cell coloured by the batch of origin (left) or the donor of origin in E-MTAB-5061 (right)." width="100%"  class="widefigure" />
<p class="caption">(\#fig:smart-raw)t-SNE plots of the merged Smart-based data sets. Each point is a cell coloured by the batch of origin (left) or the donor of origin in E-MTAB-5061 (right).</p>
</div>

To overcome this, the first step of our hierarchical merge is to remove differences between donors _within_ E-MTAB-5061.
This simply involves calling `fastMNN()` with the `batch=` argument set to the `Donor` variable. 
Note the use of `pc.input=TRUE` to avoid performing a PCA on a matrix of low-dimensional coordinates.


```r
fixed.emtab <- fastMNN(pcs$emtab, batch=sce.emtab$Donor, k=20, pc.input=TRUE)
```

It is then straightforward to merge these corrected expression values with the data from GSE86469.
This removes differences between studies and represents the second level of the merge hierarchy.


```r
mnn.smart <- fastMNN(pcs$gse86469, fixed.emtab, k=20, pc.input=TRUE)
```

This strategy eliminates the donor-based structure in the merged data (Figure \@ref(fig:smart-fixed)),
collapsing the previously distinct per-donor clusters into a single entity.


```r
set.seed(2000)
tsne.out <- Rtsne::Rtsne(mnn.smart$corrected, perplexity=30, pca=FALSE)
df <- data.frame(x=tsne.out$Y[,1], y=tsne.out$Y[,2], 
    batch=rep(c("GSE86469", "E-MTAB-5061"), 
        c(ncol(sce.gse86469), ncol(sce.emtab))),
    donor=c(rep(NA, ncol(sce.gse86469)), sce.emtab$Donor))

multiplot(
    ggplot(df) + geom_point(aes(x=x, y=y, color=donor)) +
        xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("By donor"),
    ggplot(df) + geom_point(aes(x=x, y=y, color=batch)) +
        xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("By batch"),
    cols=2)
```

<div class="figure">
<img src="multibatch_files/figure-html/smart-fixed-1.png" alt="t-SNE plots of the merged Smart-based data sets. Each point is a cell coloured by the batch of origin (left) or the donor of origin in E-MTAB-5061 (right)." width="100%"  class="widefigure" />
<p class="caption">(\#fig:smart-fixed)t-SNE plots of the merged Smart-based data sets. Each point is a cell coloured by the batch of origin (left) or the donor of origin in E-MTAB-5061 (right).</p>
</div>

## Merging the CEL-seq data

We directly merge together the two CEL-seq(2)-based data sets,
equivalent to our approach in the [previous workflow](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/batch.html).
We note that each of these data sets also contains some donor-based structure,
but we will ignore it as it does not seem to interfere with the merge.


```r
mnn.umi <- fastMNN(pcs$gse81076, pcs$gse85241, k=20, pc.input=TRUE)
```

Figure \@ref(fig:umi-merge) demonstrates that these two data sets are successfully merged.
The quality of this merge is probably due to the fact that both data sets were generated by the same provider,
combined with the reduction in technical variability that is offered by UMI-based protocols.


```r
set.seed(2000)
tsne.out <- Rtsne::Rtsne(mnn.umi$corrected, perplexity=30, pca=FALSE)
df <- data.frame(x=tsne.out$Y[,1], y=tsne.out$Y[,2], 
    batch=rep(c("GSE81076", "GSE85241"), 
        c(ncol(sce.gse81076), ncol(sce.gse85241))))

ggplot(df) + geom_point(aes(x=x, y=y, color=batch)) +
    xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("By batch")
```

<div class="figure">
<img src="multibatch_files/figure-html/umi-merge-1.png" alt="t-SNE plot of the merged CEL-seq data sets. Each point is a cell coloured by the batch of origin." width="100%" />
<p class="caption">(\#fig:umi-merge)t-SNE plot of the merged CEL-seq data sets. Each point is a cell coloured by the batch of origin.</p>
</div>

Finally, we merge the merged data sets across technologies by calling `fastMNN()` on the output of the per-technology merges.
This represents the final level of the hierarchical merge.


```r
mnn.final <- fastMNN(mnn.umi, mnn.smart, k=20, pc.input=TRUE)
```

Cells from multiple batches group together in distinct clusters in Figure \@ref(fig:overall-merge).
Each of the large clusters corresponds to a single cell type as annotated in the individual studies,
which suggests that the merge is largely satisfactory.


```r
set.seed(3000)
tsne.out <- Rtsne::Rtsne(mnn.final$corrected, perplexity=30, pca=FALSE)
df <- data.frame(x=tsne.out$Y[,1], y=tsne.out$Y[,2], 
    batch=rep(c("GSE81076", "GSE85241", "GSE86469", "E-MTAB-5061"), 
        c(ncol(sce.gse81076), ncol(sce.gse85241),
            ncol(sce.gse86469), ncol(sce.emtab))),
    type=c(rep("unknown", ncol(sce.gse81076)+ncol(sce.gse85241)), 
        sce.gse86469$CellType, sce.emtab$CellType),
    stringsAsFactors=FALSE)

# Restricting colors to certain approved cell types for visibility.
approved <- c(Acinar="#ffff00", Alpha="#ff0000", Beta="#c400ff", 
    Delta="#ff7800", Ductal="#00f5ff", Gamma="#0000ff",
    Other="#000000", unknown="grey80")
df$type[df$type=="Gamma/PP"] <- "Gamma"
df$type[!df$type %in% names(approved)] <- "Other"

multiplot(
    ggplot(df) + geom_point(aes(x=x, y=y, color=batch)) +
        xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("By batch"),
    ggplot(df) + geom_point(aes(x=x, y=y, color=type)) +
        scale_color_manual(values=approved) +
        xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("By cell type"),
    cols=2)
```

<div class="figure">
<img src="multibatch_files/figure-html/overall-merge-1.png" alt="t-SNE plot of all merged pancreas data sets. Each point is a cell coloured by the batch of origin (left) or the annotated cell type (right)." width="100%"  class="widefigure" />
<p class="caption">(\#fig:overall-merge)t-SNE plot of all merged pancreas data sets. Each point is a cell coloured by the batch of origin (left) or the annotated cell type (right).</p>
</div>

We can verify this by clustering on the corrected low-dimensional values using a graph-based method [@xu2015identification].
Each cluster contains contributions from all batches and is cleanly separated in Figure \@ref(fig:overall-cluster).


```r
g <- buildSNNGraph(mnn.final$corrected, d=NA, transposed=TRUE)
clusters <- igraph::cluster_walktrap(g)
df$cluster <- factor(clusters$membership)
table(df$cluster, df$batch) # Good mixing between batches.
```

```
##     
##      E-MTAB-5061 GSE81076 GSE85241 GSE86469
##   1          185      356      285       36
##   2          254      183      479      264
##   3          377      353      252       26
##   4          818      217      851      241
##   5           89       63      197       21
##   6           53       28      109       19
##   7          156       23      128       18
##   8            0       52        0        0
##   9            5       10       24        2
##   10          14        7       21        8
```

```r
ggplot(df) + geom_point(aes(x=x, y=y, color=cluster)) +
    xlab("t-SNE 1") + ylab("t-SNE 2") + ggtitle("By cluster")
```

<div class="figure">
<img src="multibatch_files/figure-html/overall-cluster-1.png" alt="t-SNE plot of all merged pancreas data sets. Each point is a cell coloured by the assigned cluster." width="100%" />
<p class="caption">(\#fig:overall-cluster)t-SNE plot of all merged pancreas data sets. Each point is a cell coloured by the assigned cluster.</p>
</div>

# Variance-based diagnostics 

As mentioned [previously](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/batch.html#with-diagnostics), the proportion of lost variance within each batch serves as a useful diagnostic for whether biological structure was inadvertently removed during correction.
Here, the proportion lost at each merge step is low, indicating that the variation within each batch is mostly preserved.


```r
summary(metadata(fixed.emtab)$merge.info$lost.var)
```

```
##        V1                 V2                 V3                 V4          
##  Min.   :0.001531   Min.   :0.001743   Min.   :0.001617   Min.   :0.001468  
##  1st Qu.:0.002145   1st Qu.:0.001874   1st Qu.:0.002220   1st Qu.:0.003676  
##  Median :0.004696   Median :0.006004   Median :0.004906   Median :0.006827  
##  Mean   :0.008038   Mean   :0.007645   Mean   :0.010360   Mean   :0.009625  
##  3rd Qu.:0.006884   3rd Qu.:0.007182   3rd Qu.:0.009930   3rd Qu.:0.007237  
##  Max.   :0.031898   Max.   :0.025134   Max.   :0.047892   Max.   :0.039907  
##        V5                 V6                 V7                 V8          
##  Min.   :0.001204   Min.   :0.001844   Min.   :0.001591   Min.   :0.004154  
##  1st Qu.:0.002948   1st Qu.:0.005558   1st Qu.:0.005968   1st Qu.:0.005660  
##  Median :0.004659   Median :0.006014   Median :0.006952   Median :0.009350  
##  Mean   :0.007533   Mean   :0.010315   Mean   :0.009142   Mean   :0.011774  
##  3rd Qu.:0.008815   3rd Qu.:0.007839   3rd Qu.:0.009812   3rd Qu.:0.012523  
##  Max.   :0.027470   Max.   :0.043382   Max.   :0.025755   Max.   :0.031818  
##        V9                V10          
##  Min.   :0.001657   Min.   :0.001925  
##  1st Qu.:0.004080   1st Qu.:0.003543  
##  Median :0.007327   Median :0.009597  
##  Mean   :0.008000   Mean   :0.008938  
##  3rd Qu.:0.009857   3rd Qu.:0.010073  
##  Max.   :0.022084   Max.   :0.025488
```

```r
metadata(mnn.smart)$merge.info$lost.var
```

```
##             [,1]        [,2]
## [1,] 0.006079988 0.004685284
```

```r
metadata(mnn.umi)$merge.info$lost.var
```

```
##            [,1]       [,2]
## [1,] 0.01397788 0.01112676
```

```r
metadata(mnn.final)$merge.info$lost.var
```

```
##             [,1]        [,2]
## [1,] 0.001241117 0.002401186
```

In a hierarchical merge, an additional metric is the proportion of variation lost during re-orthogonalization.
For example, consider `mnn.final` that results from merging `mnn.smart` and `mnn.umi`.
The orthogonalization step that was performed in `mnn.umi` is applied to the corrected values in `mnn.smart` and vice versa,
to ensure both inputs are processed in a consistent manner before attempting to merge them.
(In particular, it aims to avoid differences in variance that could interfere with a successful merge,
and/or compromise the utility of the proportion of variance as a diagnostic.)

Each re-orthogonalization step discards a proportion of variance that is recorded in the metadata of the `fastMNN()` output.
The intepretation of these propotions are the same as those from `merge.info$lost.var`, i.e.,
large proportions suggest that biological structure was in appropriate removed.
We can inspect the proportions at each of the "higher level" merge steps to ensure that this is not the case.


```r
summary(metadata(mnn.smart)$pre.orthog$lost.var)
```

```
##        V1                 V2           
##  Min.   :0.001631   Min.   :0.0006439  
##  1st Qu.:0.002285   1st Qu.:0.0015109  
##  Median :0.003331   Median :0.0025192  
##  Mean   :0.006580   Mean   :0.0025151  
##  3rd Qu.:0.004899   3rd Qu.:0.0029994  
##  Max.   :0.022044   Max.   :0.0059944
```

```r
summary(metadata(mnn.final)$pre.orthog$lost.var)
```

```
##        V1                 V2           
##  Min.   :0.000000   Min.   :0.0005297  
##  1st Qu.:0.003042   1st Qu.:0.0007589  
##  Median :0.003519   Median :0.0015502  
##  Mean   :0.006398   Mean   :0.0021724  
##  3rd Qu.:0.007437   3rd Qu.:0.0026115  
##  Max.   :0.027724   Max.   :0.0072834
```

# Obtaining corrected expression values

If `fastMNN()` was run on low-dimensional inputs, only the low-dimensional output will be reported.
Nonetheless, users can obtain per-gene corrected values by manually computing the cross-product using the PCA rotation vectors.
For example, the code below obtains corrected expression values for _GCG_ from our hierarchical merge.


```r
rotations <- metadata(pcs)$rotation
cor.exp <- tcrossprod(mnn.final$corrected,
    rotations["ENSG00000115263",,drop=FALSE])
summary(cor.exp)
```

```
##  ENSG00000115263  
##  Min.   :-5.1092  
##  1st Qu.:-1.8695  
##  Median :-0.9161  
##  Mean   : 0.1680  
##  3rd Qu.: 2.9965  
##  Max.   : 6.2226
```

Explicit calculation of all per-gene corrected values is probably ill-advised as this would involve the construction of a dense matrix.
This may be prohibitively memory-consuming for large data sets that are otherwise representable as sparse matrices.
Rather, corrected values can be computed for specific genes as they are needed, e.g., using the `LowRankMatrix` class.


```r
lrm <- LowRankMatrix(rotations, mnn.final$corrected)
lrm
```

```
## <8697 x 6224> LowRankMatrix object of type "double":
##                          D2ex_1          D2ex_2          D2ex_3 ...
## ENSG00000148584   -0.1986679546   -0.0328872638   -0.0443295645   .
## ENSG00000175899   -0.0603128831    0.0055706534   -0.1520464816   .
## ENSG00000094914   -0.0008825524   -0.2574246388   -0.1147313175   .
## ENSG00000114771    0.6259440629    0.4149066801    0.4024239683   .
## ENSG00000103591    0.0894946971   -0.0738051983   -0.0497041029   .
##             ...               .               .               .   .
## ENSG00000122952      0.15210901      0.03300412      0.03116294   .
## ENSG00000198455     -0.02226360     -0.06031804     -0.08798175   .
## ENSG00000070476      0.27552491      0.56472269      0.48913136   .
## ENSG00000159840      0.46303791      0.01769438     -0.14231625   .
## ENSG00000036549     -0.20340786     -0.14558450     -0.12167209   .
##                 HP1526901T2D_A8
## ENSG00000148584     0.165198080
## ENSG00000175899    -0.004638561
## ENSG00000094914    -0.178480012
## ENSG00000114771    -0.061191434
## ENSG00000103591    -0.101374717
##             ...               .
## ENSG00000122952     -0.09349929
## ENSG00000198455      0.02606419
## ENSG00000070476      0.01556047
## ENSG00000159840      0.07999075
## ENSG00000036549     -0.37610220
```

Users are referred to the [previous workflow](https://bioconductor.org/packages/3.10/simpleSingleCell/vignettes/batch.html) for some caveats on the direct use of the corrected expression values.

# References
