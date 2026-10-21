# Canine RNA-seq Differential-Expression Benchmark

Benchmarking DESeq2, edgeR, and limma-voom using simulated and canine lymphoma RNA-seq data.

## Project overview

Different statistical methods can produce different results when analyzing the same RNA-seq count data. This project investigates how three widely used methods - DESeq2, edgeR, and limma-voom - compare when identifying differentially expressed genes.

The methods will first be tested on simulated RNA-seq data where the truly differentially expressed genes are known. This makes it possible to measure how accurately each method recovers the known signal and how often it reports genes that are not truly differentially expressed. The methods will then be applied to a real canine lymphoma dataset to compare their results in a biologically relevant setting.

## Background

Canine lymphoma is one of the most common cancers diagnosed in dogs. It is a disease of lymphocytes, which are cells involved in the immune system. Because lymphoma affects immune-cell biology, studying gene-expression patterns may help researchers identify genes and biological processes associated with the disease. This makes canine lymphoma a useful biological context for evaluating RNA-seq differential-expression methods.

RNA sequencing, or RNA-seq, provides a way to study these gene-expression patterns. RNA molecules reflect which genes are active in a cell or tissue. During an RNA-seq experiment, RNA is converted into sequencing fragments, and those fragments are assigned to genes. The number of fragments assigned to each gene produces a count matrix, with genes arranged as rows and samples arranged as columns.

For example, a simplified count matrix might look like this:

| Gene | sample_1 | sample_2 | sample_3 |
|---|---:|---:|---:|
| gene_1 | 120 | 98 | 210 |
| gene_2 | 45 | 51 | 49 |
| gene_3 | 300 | 275 | 410 |

Each row represents a gene. Each sample column contains the number of sequencing fragments assigned to that gene.

Once a count matrix has been created, researchers can use it to ask whether gene activity differs between biological groups. This process is called differential-expression analysis. In this project, the main comparison is between blood samples from dogs with lymphoma and samples from dogs in the non-lymphoma comparison group.

A higher average count in one group does not automatically mean that a gene is truly differentially expressed. The difference could be influenced by natural variation among the biological samples or by differences in how deeply each sample was sequenced. The analysis must also account for the fact that thousands of genes are tested at the same time. For these reasons, differential-expression methods use statistical models to separate meaningful group differences from background variation.

DESeq2, edgeR, and limma-voom all address this problem, but they model count data and variation in different ways. Applying the three methods to the same datasets makes it possible to examine how those statistical choices affect the genes identified as differentially expressed. The real canine dataset provides biological context, while simulated data make it possible to measure accuracy because the truly differentially expressed genes are known in advance.

## What is benchmarking?

Benchmarking means comparing different methods under the same conditions. In this project, DESeq2, edgeR, and limma-voom will analyze the same RNA-seq count data and receive the same information about the sample groups.

Their results will be compared by examining which genes they identify and how consistently they rank those genes. The simulated datasets will also make it possible to evaluate whether each method recovers the genes that were intentionally made different. Later in the project, runtime and memory use will also be examined.

The purpose of benchmarking is not to declare that one method is always the best. Instead, it is to understand how different statistical approaches behave when they analyze the same data.

## Why use simulated data?

The main challenge in benchmarking differential-expression methods is that the correct answer is usually unknown for a real biological dataset. Even if several methods identify the same genes, that agreement does not prove that the genes are truly differentially expressed. Without known truth, it is difficult to measure accuracy directly.

Simulated data provides a controlled way to address this problem. A researcher can create count data in which a specific set of genes is intentionally made more or less expressed between groups. Because those genes are known in advance, the results from DESeq2, edgeR, and limma-voom can be compared with the simulated truth. This makes it possible to calculate measures such as precision, recall, and false discovery rate.

## Tools

The project uses three differential-expression workflows and one simulation package. DESeq2, edgeR, and limma-voom are the methods being compared. compcodeR will support the simulation part of the project by creating datasets with known differential-expression truth.

### DESeq2

DESeq2 is an R/Bioconductor package designed for differential-expression analysis of unnormalized RNA-seq count data. It uses a negative-binomial model to account for variability in sequencing counts. Its results include estimated fold changes, p-values, and adjusted p-values for each gene.

### edgeR

edgeR is an R/Bioconductor package for analyzing sequencing count data. It also models count variability using a negative-binomial framework, but it uses a different statistical workflow from DESeq2. This project uses edgeR's quasi-likelihood approach for the preliminary analysis.

### limma-voom

limma-voom is a workflow that uses the limma package together with the voom transformation. Voom transforms the count data while estimating the precision of each observation. The limma framework then fits linear models and applies empirical Bayes moderation to make the results more stable, especially when the number of samples is limited.

### compcodeR

compcodeR is a simulation package for generating RNA-seq count data. It is not one of the three differential-expression methods being benchmarked. Instead, it will be used to create simulated datasets in which the truly differentially expressed genes are known. These known genes will serve as the reference truth when method accuracy is evaluated.

## Installation

The required software can be installed from within R. BiocManager is used to install packages from Bioconductor, which is the software repository that contains DESeq2, edgeR, limma, and compcodeR.

Install BiocManager first:

```r
install.packages("BiocManager")
```
Then install the required Bioconductor packages:

```r
BiocManager::install(c("DESeq2", "edgeR", "limma", "compcodeR"))
```
The current scripts were tested using:

- R 4.5.2
- DESeq2 1.50.2
- edgeR 4.8.2
- limma 3.66.0
- compcodeR 1.46.0

All four packages loaded successfully during the package diagnostic test.

## Data

This project uses two types of RNA-seq count data. Simulated data will be used to evaluate method accuracy when the true differentially expressed genes are known. Real canine lymphoma data will be used to compare the methods in a biological setting where the complete truth is not known.

### Real canine lymphoma dataset

The real-data analysis will use GEO dataset GSE196186. This dataset contains whole-blood RNA-seq count data from seven dogs with lymphoma and seven dogs in a non-lymphoma comparison group.

The non-lymphoma samples should not be interpreted as seven healthy controls because the comparison group includes dogs with other medical conditions. This is an important limitation when interpreting biological differences between the groups.

The count matrix can be downloaded from the [GEO GSE196186 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196186). The associated study is Sheet et al. (2022), available through [PubMed](https://pubmed.ncbi.nlm.nih.gov/35835403/).

The raw count matrix is not currently stored in this repository. After downloading it from GEO, place it at:

```text
data/GSE196186_Raw_gene_counts_matrix.tsv
```

## How to run the scripts

Open `canine-rnaseq-benchmark.Rproj` in RStudio. Opening the project helps R use the correct project folder when reading files. The current scripts are stored in the main repository folder.

In RStudio, run each script from the project directory using the `source()` command. The `source()` function reads the script and runs its code.

### `01_inspect_data.R`

Before running this script, download the GSE196186 count matrix from GEO. Place the file at:
```text
data/GSE196186_Raw_gene_counts_matrix.tsv
```
Run the script with:

```r
source("01_inspect_data.R")
```
This script reads the count matrix and selects the 14 samples used in the study. It then creates a table describing each sample’s group and checks whether the data have the expected dimensions, missing values, negative counts, non-integer counts, or mismatched sample names. The results are printed in the R console.

### `02_package_diagnostics.R`

Run the package diagnostic script with:

```r
source("02_package_diagnostics.R")
```
This script prints the R version used for the project. It then loads DESeq2, edgeR, limma, and compcodeR and prints the installed version of each package. A successful run confirms that the required packages are installed and available in the current R environment.

### `03_test_de_methods.R`

Run the preliminary differential-expression test with:

```r
source("03_test_de_methods.R")
```
This script creates a small synthetic RNA-seq count matrix with 100 genes and 8 samples. Four samples are assigned to the non-lymphoma group, and four samples are assigned to the lymphoma group. The script intentionally increases the counts for the first 10 genes in the lymphoma samples. These genes are therefore the known differentially expressed genes in this preliminary test.

The script then gives the same count matrix to DESeq2, edgeR, and limma-voom. Each method produces a result table and prints its 10 highest-ranked genes. The first 10 genes should generally rank highly, although other genes may also appear because the dataset contains random variation.

#### Expected results for 03_test_de_methods.R

The three methods use different names for similar types of information.

DESeq2 reports log2FoldChange, pvalue, and padj. The fold-change column describes the estimated difference between the groups. The adjusted p-value accounts for the fact that many genes are tested at the same time.

edgeR reports logFC, logCPM, PValue, and FDR. The logFC column describes the estimated fold change, while FDR is the multiple-testing-adjusted significance value.

limma-voom reports logFC, P.Value, and adj.P.Val. These columns describe the estimated fold change, the unadjusted p-value, and the adjusted p-value.

The current script is a preliminary test of the analysis workflow. The full benchmark will later save complete result tables and compare the methods using consistent significance thresholds and accuracy measures.

## Current progress

- Created and organized the GitHub repository.
- Selected DESeq2, edgeR, limma-voom, compcodeR, and GSE196186.
- Created and tested scripts for inspecting data and checking packages.
- Recorded the R and package versions.
- Created and tested a small script that runs all three methods on synthetic data.

## Planned work

- Use compcodeR to create simulated RNA-seq datasets with known DE genes.
- Test the methods with different sample sizes, including 7 and 15 samples per group.
- Run DESeq2, edgeR, and limma-voom on the same simulated datasets.
- Compare how well each method finds the known DE genes.
- Measure runtime and memory using datasets scaled by factors of 10, 100, 1000, and 10000, with 3 replicates at each size until the tools begin to struggle or crash.
- Apply the methods to the canine lymphoma dataset and compare the results.

## References

- [DESeq2 Bioconductor vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [edgeR User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)
- [limma User's Guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
- [compcodeR Bioconductor package](https://bioconductor.org/packages/release/bioc/html/compcodeR.html)
- [GEO GSE196186 dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196186)
- [Sheet et al. (2022) canine lymphoma study](https://pubmed.ncbi.nlm.nih.gov/35835403/)
