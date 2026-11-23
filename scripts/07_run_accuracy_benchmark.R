# Compare DESeq2, edgeR, and limma-voom using simulated data with known results.
# Each dataset has 1,000 genes, including 100 genes known to be differentially
# expressed. The methods' significant-gene calls are compared with that truth.

# Run it from the main project folder so the file paths work.

# Load compcodeR for simulation.
library(compcodeR)

# Load the packages for the three differential-expression methods.
library(DESeq2)
library(edgeR)
library(limma)

# Create folders for the accuracy data and results if needed.
if (!dir.exists("data/accuracy")) {
  dir.create("data/accuracy", recursive = TRUE)}

if (!dir.exists("results/accuracy")) {
  dir.create("results/accuracy", recursive = TRUE)}

if (!dir.exists("results/accuracy/method_results")) {
  dir.create("results/accuracy/method_results", recursive = TRUE)}

# The benchmark tests 14 total samples (7 per group), 50 (25 per group),
# and 100 (50 per group), with 5 replicates at each size.
# 15 datasets and 45 method runs in total.

# Settings for the accuracy benchmark.
benchmark_settings <- data.frame(
  total_samples = c(rep(14, 5), rep(50, 5), rep(100, 5)),
  replicate = rep(1:5, times = 3))

# Start an empty table to collect metrics from each dataset.
benchmark_metrics <- data.frame()

# Repeat the steps below for every row in the settings table.
for (i in 1:nrow(benchmark_settings)) {
  total_samples <- benchmark_settings$total_samples[i]
  samples_per_group <- total_samples / 2
  replicate_id <- benchmark_settings$replicate[i]
  
  # ----------------------------
  # 1. Simulate one dataset
  # ----------------------------
  
  # Include the total sample size and replicate in the dataset name.
  dataset_name <- paste0("accuracy_total", total_samples,"_repl", replicate_id)
  
  # Set a reproducible random seed for this dataset.
  set.seed(samples_per_group * 1000 + replicate_id)
  
  # Generate 1,000 genes, including 100 known DE genes.
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
    output.file = paste0("data/accuracy/", dataset_name, ".rds"))
  
  # Check whether compcodeR considers this dataset valid.
  check_compData(simulation_data)
  
  # Check the number of genes and samples.
  dim(simulation_data@count.matrix)
  
  # Check the number of samples in each group.
  table(simulation_data@sample.annotations$condition)
  
  # Mark which genes are truly differentially expressed.
  true_de <- simulation_data@variable.annotations$differential.expression == 1
  
  # Count the known DE genes; this should be 100.
  sum(true_de)
  
  # Get the count matrix and sample information for the method analyses.
  counts <- simulation_data@count.matrix
  sample_info <- simulation_data@sample.annotations
  
  # -----------------------------------
  # 2. Prepare the sample information
  # -----------------------------------
  
  # Give the two groups readable names.
  sample_info$group <- factor(
    sample_info$condition,
    levels = c(1, 2),
    labels = c("condition_1", "condition_2"))
  
  # Match the sample information to the count-matrix columns.
  rownames(sample_info) <- colnames(counts)
  
  # Describe the comparison between the two groups.
  design <- model.matrix(~ group, data = sample_info)
  
  # ------------------
  # 3. Run DESeq2
  # ------------------
  
  # Create a DESeq2 dataset.
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = ~ group )
  
  # Run the DESeq2 analysis.
  dds <- DESeq(dds)
  
  # Compare condition_2 with condition_1.
  deseq2_results <- results(dds,
                            contrast = c("group", "condition_2", "condition_1"))
  
  # Display the 10 highest-ranked genes by adjusted p-value.
  deseq2_results <- deseq2_results[order(deseq2_results$padj), ]
  head(deseq2_results, 10)
  
  # Save all DESeq2 results for this dataset.
  write.csv(as.data.frame(deseq2_results),
            paste0("results/accuracy/method_results/", dataset_name, "_DESeq2.csv"))
  
  # -------------------
  # 4. Run edgeR
  # -------------------
  
  # Create an edgeR object.
  dge <- DGEList(counts = counts, group = sample_info$group)
  
  # Normalize library sizes and estimate dispersion.
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design)
  
  # Fit the model and test condition_2 against condition_1.
  edger_fit <- glmQLFit(dge, design)
  edger_test <- glmQLFTest(edger_fit, coef = 2)
  
  # Get all edgeR results.
  edger_results <- topTags(edger_test, n = Inf)$table
  
  # Display the 10 highest-ranked genes.
  head(edger_results, 10)
  
  # Save all edgeR results for this dataset.
  write.csv(edger_results,
            paste0("results/accuracy/method_results/", dataset_name, "_edgeR.csv"))
  
  # --------------------
  # 5. Run limma-voom
  # --------------------
  
  # Transform the count data with voom.
  voom_data <- voom(dge, design, plot = FALSE)
  
  # Fit and moderate the linear model.
  limma_fit <- lmFit(voom_data, design)
  limma_fit <- eBayes(limma_fit)
  
  # Get all limma-voom results.
  limma_results <- topTable(
    limma_fit,
    coef = 2,
    number = Inf,
    sort.by = "P")
  
  # Display the 10 highest-ranked genes.
  head(limma_results, 10)
  
  # Save all limma-voom results for this dataset.
  write.csv(limma_results,
            paste0("results/accuracy/method_results/", dataset_name, "_limma_voom.csv"))
  
  # ---------------------------------------------------
  # 6. Compare the method results with the known truth
  # ---------------------------------------------------
  
  # Use all genes in this dataset as the comparison set.
  all_genes <- rownames(counts)
  
  # Get the names of the genes known to be differentially expressed.
  true_de_genes <- rownames(counts)[true_de]
  
  # Use adjusted p-value < 0.05 as the cutoff for all methods.
  # Recall is also called sensitivity, so it is included once in the results.
  
  # Get the significant genes reported by each method.
  deseq2_called <- rownames(deseq2_results)[
    !is.na(deseq2_results$padj) & deseq2_results$padj < 0.05]
  
  edger_called <- rownames(edger_results)[
    !is.na(edger_results$FDR) & edger_results$FDR < 0.05]
  
  limma_called <- rownames(limma_results)[
    !is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.05]
  
  # These counts compare each method's calls with the known DE genes:
  # TP = true DE gene detected; FP = non-DE gene called DE;
  # FN = true DE gene missed; TN = non-DE gene correctly not called DE.
  
  # ---------------------------
  # DESeq2 performance metrics
  # ---------------------------
  
  deseq2_true_positive <- sum(deseq2_called %in% true_de_genes)
  deseq2_false_positive <- sum(!(deseq2_called %in% true_de_genes))
  deseq2_false_negative <- sum(!(true_de_genes %in% deseq2_called))
  deseq2_true_negative <- length(all_genes) -
    deseq2_true_positive - deseq2_false_positive - deseq2_false_negative
  
  deseq2_metrics <- data.frame(
    dataset = dataset_name,
    total_samples = total_samples,
    samples_per_group = samples_per_group,
    replicate = replicate_id,
    method = "DESeq2",
    true_positive = deseq2_true_positive,
    false_positive = deseq2_false_positive,
    false_negative = deseq2_false_negative,
    true_negative = deseq2_true_negative,
    precision = deseq2_true_positive /
      (deseq2_true_positive + deseq2_false_positive),
    recall = deseq2_true_positive /
      (deseq2_true_positive + deseq2_false_negative),
    specificity = deseq2_true_negative /
      (deseq2_true_negative + deseq2_false_positive),
    observed_fdr = deseq2_false_positive /
      (deseq2_true_positive + deseq2_false_positive))
  
  # --------------------------
  # edgeR performance metrics
  # --------------------------
  
  edger_true_positive <- sum(edger_called %in% true_de_genes)
  edger_false_positive <- sum(!(edger_called %in% true_de_genes))
  edger_false_negative <- sum(!(true_de_genes %in% edger_called))
  edger_true_negative <- length(all_genes) -
    edger_true_positive - edger_false_positive - edger_false_negative
  
  edger_metrics <- data.frame(
    dataset = dataset_name,
    total_samples = total_samples,
    samples_per_group = samples_per_group,
    replicate = replicate_id,
    method = "edgeR",
    true_positive = edger_true_positive,
    false_positive = edger_false_positive,
    false_negative = edger_false_negative,
    true_negative = edger_true_negative,
    precision = edger_true_positive /
      (edger_true_positive + edger_false_positive),
    recall = edger_true_positive /
      (edger_true_positive + edger_false_negative),
    specificity = edger_true_negative /
      (edger_true_negative + edger_false_positive),
    observed_fdr = edger_false_positive /
      (edger_true_positive + edger_false_positive))
  
  # -------------------------------
  # limma-voom performance metrics
  # -------------------------------
  
  limma_true_positive <- sum(limma_called %in% true_de_genes)
  limma_false_positive <- sum(!(limma_called %in% true_de_genes))
  limma_false_negative <- sum(!(true_de_genes %in% limma_called))
  limma_true_negative <- length(all_genes) -
    limma_true_positive - limma_false_positive - limma_false_negative
  
  limma_metrics <- data.frame(
    dataset = dataset_name,
    total_samples = total_samples,
    samples_per_group = samples_per_group,
    replicate = replicate_id,
    method = "limma-voom",
    true_positive = limma_true_positive,
    false_positive = limma_false_positive,
    false_negative = limma_false_negative,
    true_negative = limma_true_negative,
    precision = limma_true_positive /
      (limma_true_positive + limma_false_positive),
    recall = limma_true_positive /
      (limma_true_positive + limma_false_negative),
    specificity = limma_true_negative /
      (limma_true_negative + limma_false_positive),
    observed_fdr = limma_false_positive /
      (limma_true_positive + limma_false_positive))
  
  # Combine the three methods' metrics for this dataset.
  simulation_metrics <- rbind(
    deseq2_metrics,
    edger_metrics,
    limma_metrics )
  
  # Add these three rows to the full benchmark table.
  benchmark_metrics <- rbind(benchmark_metrics, simulation_metrics)
}

# -----------------------------------------
# 7. Display and save all benchmark metrics
# -----------------------------------------

# This table has 45 rows: 15 datasets multiplied by 3 methods.
print(benchmark_metrics)

write.csv(
  benchmark_metrics,
  "results/accuracy/accuracy_benchmark_metrics.csv",
  row.names = FALSE)