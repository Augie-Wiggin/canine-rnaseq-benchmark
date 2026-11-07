# Canine RNA-seq Differential-Expression Benchmark

Benchmarking DESeq2, edgeR, and limma-voom using simulated and canine lymphoma RNA-seq data.

## Project overview

This project compares three R/Bioconductor methods for identifying differentially expressed genes from RNA-seq count data. It investigates how DESeq2, edgeR, and limma-voom compare when analyzing the same datasets.

The methods will first be tested using simulated RNA-seq data where the truly differentially expressed genes are known. This makes it possible to evaluate how well each method identifies the expected genes. The methods will then be applied to a real canine lymphoma dataset so their results can be compared in a biological setting.

## Background

Canine lymphoma is one of the most common cancers diagnosed in dogs. It affects lymphocytes, which are cells involved in the immune system. Studying gene-expression patterns in lymphoma can help researchers examine which genes may be associated with the disease.

RNA sequencing, or RNA-seq, measures gene activity by sequencing RNA molecules and counting the fragments assigned to each gene. These counts are organized into a matrix, with genes as rows and biological samples as columns.

For example, a simplified count matrix might look like this:

| Gene | sample_1 | sample_2 | sample_3 |
|---|---:|---:|---:|
| gene_1 | 120 | 98 | 210 |
| gene_2 | 45 | 51 | 49 |
| gene_3 | 300 | 275 | 410 |

Each row represents a gene. Each sample column contains the number of sequencing fragments assigned to that gene.

Differential-expression analysis uses this count matrix to test whether gene expression differs between biological groups. In this project, the main comparison is between blood samples from dogs with lymphoma and dogs in the non-lymphoma comparison group.

A higher average count does not automatically mean that a gene is differentially expressed. Statistical methods must account for variation between samples and adjust the results because many genes are tested at the same time.

## Benchmarking approach

Benchmarking means comparing different methods under the same conditions. In this project, DESeq2, edgeR, and limma-voom will analyze the same RNA-seq count data and use the same sample-group information. Their results will be compared based on the genes they identify and how closely they agree. Runtime and memory use will also be examined later in the project. The goal is not to claim that one method is always best, but to understand how the methods behave when analyzing the same data.

The correct answer is usually unknown for a real biological dataset, which makes accuracy difficult to measure directly. Simulated data provides a controlled solution because the truly differentially expressed genes are known in advance. The results from each method can then be compared with that known truth using measures such as precision, recall, and false discovery rate.

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
This script reads the count matrix and selects the 14 samples used in the study. It then creates a table describing each sample’s group and checks whether the data has the expected dimensions, missing values, negative counts, non-integer counts, or mismatched sample names. The results are printed in the R console.

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

#### Expected results for `03_test_de_methods.R`

The main output columns are summarized below:

| Method | Important columns | What they mean |
|---|---|---|
| DESeq2 | `log2FoldChange`, `pvalue`, `padj` | The fold-change column shows the direction and size of the difference. `pvalue` is the unadjusted result, and `padj` is the adjusted result. |
| edgeR | `logFC`, `logCPM`, `PValue`, `FDR` | `logFC` shows the estimated fold change. `logCPM` describes average expression. `PValue` is unadjusted, and `FDR` is adjusted. |
| limma-voom | `logFC`, `AveExpr`, `P.Value`, `adj.P.Val` | `logFC` shows the estimated fold change. `AveExpr` describes average expression after voom. `P.Value` is unadjusted, and `adj.P.Val` is adjusted. |

In the current scripts, all three methods use the Benjamini-Hochberg adjustment by default. The adjusted values account for testing many genes at the same time. Genes with an adjusted value below 0.05 will generally be considered statistically significant.

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
