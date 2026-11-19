# Run the full 20-dataset accuracy benchmark.
# For each dataset, this script simulates data, runs all three methods,
# and compares their results with compcodeR's known DE genes.
# The benchmark includes 10 datasets with 7 samples per group and
# 10 datasets with 15 samples per group.

# Load compcodeR for simulation.
library(compcodeR)

# Load the packages for the three differential-expression methods.
library(DESeq2)
library(edgeR)
library(limma)

# Create folders for the benchmark data and results if needed.
if (!dir.exists("data/simulation")) {
  dir.create("data/simulation", recursive = TRUE)}

if (!dir.exists("results/simulation")) {
  dir.create("results/simulation", recursive = TRUE)}

# Settings for the full 20-dataset accuracy benchmark.
benchmark_settings <- data.frame(
  samples_per_group = c(rep(7, 10), rep(15, 10)),
  replicate = rep(1:10, times = 2)
)

# Start an empty table to collect metrics from each dataset.
benchmark_metrics <- data.frame()

# Repeat the workflow for each row in the settings table.
for (i in 1:nrow(benchmark_settings)) {
  samples_per_group <- benchmark_settings$samples_per_group[i]
  replicate_id <- benchmark_settings$replicate[i]

# --------------------------------
# Simulate one benchmark dataset
# --------------------------------

# Create a name from the sample size and replicate number.
dataset_name <- paste0(
  "simulation_n", samples_per_group,
  "_repl", replicate_id)

# Set a different, reproducible seed for each dataset.
set.seed(samples_per_group * 1000 + replicate_id)

# Generate a dataset with 1,000 genes and 100 known DE genes.
simulation_data <- generateSyntheticData(
  dataset = dataset_name,
  n.vars = 1000,
  samples.per.cond = samples_per_group,
  n.diffexp = 100,
  repl.id = replicate_id,
  seqdepth = 1e7,
  fraction.upregulated = 0.5,
  between.group.diffdisp = FALSE,
  filter.threshold.total = 1,
  filter.threshold.mediancpm = 0,
  fraction.non.overdispersed = 0,
  output.file = paste0(
    "data/simulation/", dataset_name, ".rds"))

# Check whether compcodeR considers the simulated data valid.
check_compData(simulation_data)

# Check the number of genes and samples.
dim(simulation_data@count.matrix)

# Check the number of samples in each group.
table(simulation_data@sample.annotations$condition)

# Identify the genes labeled as truly differentially expressed.
true_de <- simulation_data@variable.annotations$differential.expression == 1

# Count the known DE genes; this should be 100.
sum(true_de)

# Extract the count matrix and sample information for the method analyses.
counts <- simulation_data@count.matrix
sample_info <- simulation_data@sample.annotations

# --------------------------------
# Prepare sample information
# --------------------------------

# Give the two compcodeR conditions readable names.
sample_info$group <- factor(
  sample_info$condition,
  levels = c(1, 2),
  labels = c("condition_1", "condition_2"))

# Match the sample-information rows to the count-matrix columns.
rownames(sample_info) <- colnames(counts)

# Create a design matrix describing the two groups.
design <- model.matrix(~ group, data = sample_info)

# --------------------------------
# DESeq2
# --------------------------------

# Create a DESeq2 dataset.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ group)

# Run the DESeq2 analysis.
dds <- DESeq(dds)

# Compare condition_2 with condition_1.
deseq2_results <- results(
  dds,
  contrast = c("group", "condition_2", "condition_1"))

# Display the 10 highest-ranked genes by adjusted p-value.
deseq2_results <- deseq2_results[order(deseq2_results$padj), ]
head(deseq2_results, 10)

# Save the complete DESeq2 results.
write.csv(
  as.data.frame(deseq2_results),
  paste0("results/simulation/", dataset_name, "_DESeq2.csv"))

# --------------------------------
# edgeR
# --------------------------------

# Create an edgeR object.
dge <- DGEList(
  counts = counts,
  group = sample_info$group )

# Normalize library sizes and estimate dispersion.
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)

# Fit the edgeR model and test condition_2 against condition_1.
edger_fit <- glmQLFit(dge, design)
edger_test <- glmQLFTest(edger_fit, coef = 2)

# Get the complete edgeR result table.
edger_results <- topTags(edger_test, n = Inf)$table

# Display the 10 highest-ranked genes.
head(edger_results, 10)

# Save the complete edgeR results.
write.csv(
  edger_results,
  paste0("results/simulation/", dataset_name, "_edgeR.csv"))

# --------------------------------
# limma-voom
# --------------------------------

# Transform the count data with voom.
voom_data <- voom(dge, design, plot = FALSE)

# Fit and moderate the linear model.
limma_fit <- lmFit(voom_data, design)
limma_fit <- eBayes(limma_fit)

# Get the complete limma-voom result table.
limma_results <- topTable(
  limma_fit,
  coef = 2,
  number = Inf,
  sort.by = "P")

# Display the 10 highest-ranked genes.
head(limma_results, 10)

# Save the complete limma-voom results.
write.csv(
  limma_results,
  paste0("results/simulation/", dataset_name, "_limma_voom.csv"))

# --------------------------------
# Calculate performance metrics
# --------------------------------

# Use the genes in the count matrix as the full set of genes.
all_genes <- rownames(counts)

# Get the names of the genes known to be differentially expressed.
true_de_genes <- rownames(counts)[true_de]

# Use adjusted p-value < 0.05 as the significance cutoff for all methods.

# Identify genes called significant by DESeq2.
deseq2_called <- rownames(deseq2_results)[
  !is.na(deseq2_results$padj) & deseq2_results$padj < 0.05]

# Identify genes called significant by edgeR.
edger_called <- rownames(edger_results)[
  !is.na(edger_results$FDR) & edger_results$FDR < 0.05]

# Identify genes called significant by limma-voom.
limma_called <- rownames(limma_results)[
  !is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.05]

# --------------------------------
# DESeq2 performance metrics
# --------------------------------

deseq2_true_positive <- sum(deseq2_called %in% true_de_genes)
deseq2_false_positive <- sum(!(deseq2_called %in% true_de_genes))
deseq2_false_negative <- sum(!(true_de_genes %in% deseq2_called))
deseq2_true_negative <- length(all_genes) -
  deseq2_true_positive - deseq2_false_positive - deseq2_false_negative

deseq2_metrics <- data.frame(
  dataset = dataset_name,
  method = "DESeq2",
  true_positive = deseq2_true_positive,
  false_positive = deseq2_false_positive,
  false_negative = deseq2_false_negative,
  true_negative = deseq2_true_negative,
  precision = deseq2_true_positive /
    (deseq2_true_positive + deseq2_false_positive),
  recall = deseq2_true_positive /
    (deseq2_true_positive + deseq2_false_negative),
  observed_fdr = deseq2_false_positive /
    (deseq2_true_positive + deseq2_false_positive))

# --------------------------------
# edgeR performance metrics
# --------------------------------

edger_true_positive <- sum(edger_called %in% true_de_genes)
edger_false_positive <- sum(!(edger_called %in% true_de_genes))
edger_false_negative <- sum(!(true_de_genes %in% edger_called))
edger_true_negative <- length(all_genes) -
  edger_true_positive - edger_false_positive - edger_false_negative

edger_metrics <- data.frame(
  dataset = dataset_name,
  method = "edgeR",
  true_positive = edger_true_positive,
  false_positive = edger_false_positive,
  false_negative = edger_false_negative,
  true_negative = edger_true_negative,
  precision = edger_true_positive /
    (edger_true_positive + edger_false_positive),
  recall = edger_true_positive /
    (edger_true_positive + edger_false_negative),
  observed_fdr = edger_false_positive /
    (edger_true_positive + edger_false_positive))

# --------------------------------
# limma-voom performance metrics
# --------------------------------

limma_true_positive <- sum(limma_called %in% true_de_genes)
limma_false_positive <- sum(!(limma_called %in% true_de_genes))
limma_false_negative <- sum(!(true_de_genes %in% limma_called))
limma_true_negative <- length(all_genes) -
  limma_true_positive - limma_false_positive - limma_false_negative

limma_metrics <- data.frame(
  dataset = dataset_name,
  method = "limma-voom",
  true_positive = limma_true_positive,
  false_positive = limma_false_positive,
  false_negative = limma_false_negative,
  true_negative = limma_true_negative,
  precision = limma_true_positive /
    (limma_true_positive + limma_false_positive),
  recall = limma_true_positive /
    (limma_true_positive + limma_false_negative),
  observed_fdr = limma_false_positive /
    (limma_true_positive + limma_false_positive))

# Combine the metrics from the three methods into one table.
simulation_metrics <- rbind(
  deseq2_metrics,
  edger_metrics,
  limma_metrics)

# Add this dataset's three method rows to the growing table.
benchmark_metrics <- rbind(benchmark_metrics, simulation_metrics)
}

# Display and save the combined metrics from both datasets.
print(benchmark_metrics)

write.csv(benchmark_metrics,
  "results/simulation/simulation_benchmark_metrics.csv",
  row.names = FALSE)
