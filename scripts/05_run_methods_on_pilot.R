# Run DESeq2, edgeR, and limma-voom on the same compcodeR pilot dataset.
# The simulation has known differentially expressed genes.
# This script saves the full results for Script 6 to compare with that truth.

# Load the packages used for the three analyses.
library(DESeq2)
library(edgeR)
library(limma)

# Read the simulated dataset saved by Script 4.
pilot_data <- readRDS("data/compcodeR_pilot_dataset.rds")

# Extract the count matrix and sample information from the compcodeR object.
# The @ symbol accesses a named part, or "slot," of this object.
counts <- pilot_data@count.matrix
sample_info <- pilot_data@sample.annotations

# Check the data dimensions, group labels, and first few counts.
dim(counts)
table(sample_info$condition)
head(counts)

# Replace compcodeR's numeric group labels with readable names.
# Condition 2 will be compared with Condition 1.
sample_info$group <- factor(
  sample_info$condition,
  levels = c(1, 2),
  labels = c("condition_1", "condition_2"))

# Match each sample-information row to its count-matrix column.
rownames(sample_info) <- colnames(counts)

# Describe the two groups for edgeR and limma-voom.
design <- model.matrix(~ group, data = sample_info)

# -------------------------
# DESeq2
# -------------------------

# Create a DESeq2 object from the counts and sample information.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ group)

# Run the DESeq2 analysis.
dds <- DESeq(dds)

# Test whether expression differs between condition 2 and condition 1.
deseq2_results <- results( dds,
  contrast = c("group", "condition_2", "condition_1"))

# Sort the results by adjusted p-value. 
# Put the smallest adjusted p-values at the top.
deseq2_results <- deseq2_results[order(deseq2_results$padj),]

# Show the first 10 rows as a preview.
# The complete results table is saved near the end of this script.
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

# Extract the complete edgeR results table.
# $table selects the table stored inside the topTags result.
edger_results <- topTags(edger_test, n = Inf)$table

# Display the 10 highest-ranked genes.
head(edger_results, 10)

# -------------------------
# limma-voom
# -------------------------

# Use the normalized edgeR object as input to voom.
# voom prepares count data for limma's linear-model analysis.
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
  dir.create("results/pilot", recursive = TRUE)}

# Save the complete result tables for Script 6.
# Row names are kept so the gene IDs are included in the CSV files.

write.csv(as.data.frame(deseq2_results),
  "results/pilot/deseq2_pilot_results.csv")

write.csv(edger_results,
  "results/pilot/edger_pilot_results.csv")

write.csv(limma_results,
  "results/pilot/limma_voom_pilot_results.csv")