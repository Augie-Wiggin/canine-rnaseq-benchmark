# This script tests DESeq2, edgeR, and limma-voom on a small synthetic 
# RNA-seq count dataset. The analysis workflows were adapted from 
# official Bioconductor documentation.

# Create a small synthetic RNA-seq count dataset.
set.seed(123)

n_genes <- 100
n_samples <- 8

# Generate random RNA-seq-like counts. 
# Each cell represents the number of reads assigned to one gene in one sample.

counts <- matrix(
  rnbinom(n_genes * n_samples, mu = 100, size = 10),
  nrow = n_genes,
  ncol = n_samples)

# Give the genes and samples simple names.
rownames(counts) <- paste0("gene_", 1:n_genes)
colnames(counts) <- paste0("sample_", 1:n_samples)

# Check the first few rows and the matrix dimensions.
head(counts)
dim(counts)

# Create information about which group each sample belongs to.
# The first four samples are non-lymphoma and the last four samples are lymphoma.

sample_info <- data.frame(
  sample = colnames(counts),
  group = factor(c(rep("non_lymphoma", 4), rep("lymphoma", 4)),
  levels = c("non_lymphoma", "lymphoma")))

# Use sample names as row names.
rownames(sample_info) <- sample_info$sample

# Display the sample information.
sample_info

# Make the first 10 genes more highly expressed in the lymphoma samples.
# These will be the true differentially expressed genes.

counts[1:10, 5:8] <- counts[1:10, 5:8] * 3

# Check the average counts for the first gene in each group.
mean(counts[1, 1:4])
mean(counts[1, 5:8])

# Load DESeq2.
library(DESeq2)

# Create a DESeq2 dataset from the counts and sample information.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ group)

# Run the DESeq2 analysis.
dds <- DESeq(dds)

# Compare lymphoma samples with non-lymphoma samples.
deseq2_results <- results(
  dds, contrast = c("group", "lymphoma", "non_lymphoma"))

# Sort results by adjusted p-value (padj).
deseq2_results <- deseq2_results[order(deseq2_results$padj),]

# Display the 10 most significant genes.
head(deseq2_results, 10)

# Run edgeR.
library(edgeR)

# Create an edgeR object from the count matrix.
dge <- DGEList(
  counts = counts,
  group = sample_info$group)

# Normalize library sizes between samples.
dge <- normLibSizes(dge)

# Create a design matrix describing the two groups.
design <- model.matrix(~ group, data = sample_info)

# Fit the edgeR model.
fit <- glmQLFit(dge, design)

# Test lymphoma samples against non-lymphoma samples.
# Coefficient 2 represents the group comparison.
edger_test <- glmQLFTest(fit, coef = 2)

# Display the 10 most significant genes.
edger_results <- topTags(edger_test, n = 10)$table

edger_results

# Run limma-voom.
library(limma)

# Transform the count data using voom.
voom_data <- voom(
  dge,
  design,
  plot = FALSE)

# Fit a linear model to the transformed data.
limma_fit <- lmFit(voom_data, design)

# Apply empirical Bayes moderation.
limma_fit <- eBayes(limma_fit)

# Display the 10 most significant genes.
limma_results <- topTable(
  limma_fit,
  coef = 2,
  number = 10,
  sort.by = "P")

limma_results
