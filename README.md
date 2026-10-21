# Canine RNA-seq Differential-Expression Benchmark
Benchmarking DESeq2, edgeR, and limma-voom using simulated and canine lymphoma RNA-seq data.

## Project goal
This project benchmarks DESeq2, edgeR, and limma-voom for identifying differentially expressed genes from RNA-seq count data. The methods will first be evaluated using simulated data where the truly differentially expressed genes are known, and then applied to a real canine lymphoma dataset. The project connects statistical method evaluation with biologically relevant questions about canine lymphoma.

## Tools
This project compares three R/Bioconductor methods for RNA-seq differential-expression analysis:

- DESeq2
- edgeR
- limma-voom

The project will also use compcodeR to generate simulated RNA-seq count data with known differential-expression truth.

## Data

The real-data analysis will use GSE196186, a canine whole-blood RNA-seq count dataset containing seven lymphoma samples and seven non-lymphoma samples. The count matrix is available from the [GEO GSE196186 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196186). The associated study is Sheet et al. (2022), available through [PubMed](https://pubmed.ncbi.nlm.nih.gov/35835403/). It is important to note that the non-lymphoma samples should not be interpreted as seven healthy controls because the group includes dogs with other medical conditions. Simulated datasets will be generated with known differential-expression truth for evaluating method accuracy.

## Installation

This project uses R and Bioconductor packages.

Install BiocManager and the required packages in R:

```r
install.packages("BiocManager")

BiocManager::install(
  c("DESeq2", "edgeR", "limma", "compcodeR")
)
```
Tested environment:

- R 4.5.2
- DESeq2 1.50.2
- edgeR 4.8.2
- limma 3.66.0
- compcodeR 1.46.0

All four packages loaded successfully during the diagnostic test.

Package installation and version information will be tested and documented as development continues.

## How to run the current script

Open the `canine-rnaseq-benchmark.Rproj` project in RStudio.

Before running the script, download the count matrix from GEO and place it at `data/GSE196186_Raw_gene_counts_matrix.tsv`.

```text
data/GSE196186_Raw_gene_counts_matrix.tsv
```

Run the data-inspection script:

```r
source("01_inspect_data.R")
```

The script reads the GSE196186 count matrix, selects the 14 study samples, creates the sample-information table, and checks:

- the dimensions of the count table,
- missing values,
- nonnegative count values,
- integer count values,
- sample names and ordering,
- the number of samples in each group.

The script currently prints the sample-information table and diagnostic results to the R console.

## Current progress

- Created the GitHub repository.
- Chose DESeq2, edgeR, and limma-voom for comparison.
- Selected compcodeR for simulation.
- Downloaded the GSE196186 count matrix.
- Created and tested a basic R data-inspection script.
- Recorded R and package versions and added a package-diagnostics script.

## Planned work

- Install and test the required R packages.
- Generate simulated RNA-seq data with known differentially expressed genes.
- Run the three methods on the simulated data.
- Compare accuracy, runtime, and memory use.
- Apply the methods to the canine lymphoma dataset.
- Compare and interpret the results.
