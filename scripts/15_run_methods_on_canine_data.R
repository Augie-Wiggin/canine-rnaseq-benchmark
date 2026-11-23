# Run DESeq2, edgeR, and limma-voom on the canine RNA-seq data.
# Compare 7 lymphoma samples with 7 non-lymphoma samples.
# Save the full results from each method for Script 16.

# Load the packages for the three methods.
library(DESeq2)
library(edgeR)
library(limma)

# Run Script 01 to read and check the canine data.
source("scripts/01_inspect_data.R")

# Make group a factor, with non-lymphoma as the reference group.
sample_info$group <- factor(
  sample_info$group,
  levels = c("non_lymphoma", "lymphoma"))

# Make sure the sample information matches the count columns.
rownames(sample_info) <- colnames(counts)

# Describe the comparison between the two groups.
design <- model.matrix(~ group, data = sample_info)

# Check the number of samples in each group.
print(table(sample_info$group))

# Display the design matrix.
print(design)

# -------------------------
# DESeq2 analysis
# -------------------------

# Create a DESeq2 dataset using the canine counts and sample information.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ group)

# Fit the DESeq2 model to the data.
dds <- DESeq(dds)

# Compare lymphoma with non-lymphoma.
# Set DESeq2's filtering to match the adjusted p-value threshold
# used in Script 16.
deseq2_results <- results(
  dds,
  contrast = c("group", "lymphoma", "non_lymphoma"),
  alpha = 0.05)

# Sort results by adjusted p-value and display the top 10 genes.
deseq2_results <- deseq2_results[order(deseq2_results$padj), ]
print(head(deseq2_results, 10))

# Create the canine-results folder if it does not exist.
if (!dir.exists("results/canine_data")) {
  dir.create("results/canine_data", recursive = TRUE)}

# Save the complete DESeq2 results table.
write.csv(
  as.data.frame(deseq2_results),
  "results/canine_data/deseq2_canine_results.csv")

# -------------------------
# edgeR analysis
# -------------------------

# Create an edgeR object from the count matrix.
dge <- DGEList(
  counts = counts,
  group = sample_info$group)

# Normalize library sizes and estimate dispersion.
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)

# Fit the model and test lymphoma against non-lymphoma.
fit <- glmQLFit(dge, design)
edger_test <- glmQLFTest(fit, coef = 2)

# Get the full results table and display the top 10 genes.
edger_results <- topTags(edger_test, n = Inf)$table
print(head(edger_results, 10))

# Save the complete edgeR results table.
write.csv(
  edger_results,
  "results/canine_data/edger_canine_results.csv")

# -------------------------
# limma-voom analysis
# -------------------------

# Make an edgeR object from the canine counts.
dge <- DGEList(
  counts = counts,
  group = sample_info$group)

# Normalize the library sizes.
dge <- calcNormFactors(dge)

# Use voom to prepare the counts for limma.
voom_data <- voom(
  dge,
  design,
  plot = FALSE)

# Fit the model and test the group difference.
limma_fit <- lmFit(voom_data, design)
limma_fit <- eBayes(limma_fit)

# Get the results for lymphoma compared with non-lymphoma.
limma_results <- topTable(
  limma_fit,
  coef = 2,
  number = Inf,
  sort.by = "P")

# Display the top 10 genes.
print(head(limma_results, 10))

write.csv(
  limma_results,
  "results/canine_data/limma_voom_canine_results.csv")
