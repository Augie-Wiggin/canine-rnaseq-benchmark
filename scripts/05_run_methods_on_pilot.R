# Run DESeq2, edgeR, and limma-voom
# on the compcodeR pilot dataset.

# Load the required packages.
library(DESeq2)
library(edgeR)
library(limma)

# Read the simulated dataset.
pilot_data <- readRDS("data/compcodeR_pilot_n7_repl1.rds")

# Extract the count matrix and sample information.
counts <- pilot_data@count.matrix
sample_info <- pilot_data@sample.annotations

# Check the simulated data.
dim(counts)
table(sample_info$condition)
head(counts)

# Give the two compcodeR conditions readable names.
sample_info$group <- factor(
  sample_info$condition,
  levels = c(1, 2),
  labels = c("condition_1", "condition_2"))

# Match the sample-information rows to the count-matrix columns.
rownames(sample_info) <- colnames(counts)

# Create a design matrix for the two conditions.
design <- model.matrix(~ group, data = sample_info)

# -------------------------
# DESeq2
# -------------------------

# Create a DESeq2 dataset.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ group)

# Run the DESeq2 analysis.
dds <- DESeq(dds)

# Compare condition 2 with condition 1.
deseq2_results <- results( dds,
  contrast = c("group", "condition_2", "condition_1"))

# Sort the results by adjusted p-value.
deseq2_results <- deseq2_results[order(deseq2_results$padj),]

# Display the 10 highest-ranked genes.
head(deseq2_results, 10)

# -------------------------
# edgeR
# -------------------------

# Create an edgeR object.
dge <- DGEList(
  counts = counts,
  group = sample_info$group)

# Normalize library sizes and estimate dispersion.
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)

# Fit the edgeR model.
fit <- glmQLFit(dge, design)

# Test condition 2 against condition 1.
edger_test <- glmQLFTest( fit, coef = 2)

# Display the 10 highest-ranked genes.
edger_results <- topTags(edger_test, n = Inf)$table

# Display the 10 highest-ranked genes.
head(edger_results, 10)

# -------------------------
# limma-voom
# -------------------------

# Transform the count data with voom.
voom_data <- voom(
  dge,
  design,
  plot = FALSE)

# Fit and moderate the linear model.
limma_fit <- lmFit(voom_data, design)

limma_fit <- eBayes(limma_fit)

# Save the complete limma-voom result table.
limma_results <- topTable(
  limma_fit,
  coef = 2,
  number = Inf,
  sort.by = "P")

# Display the 10 highest-ranked genes.
head(limma_results, 10)

# Create the pilot-results folder if needed.
if (!dir.exists("results/pilot")) {
  dir.create("results/pilot", recursive = TRUE)
}

# Save the complete result tables as CSV files.

write.csv(as.data.frame(deseq2_results),
  "results/pilot/deseq2_pilot_results.csv")

write.csv(edger_results,
  "results/pilot/edger_pilot_results.csv")

write.csv(limma_results,
  "results/pilot/limma_voom_pilot_results.csv")