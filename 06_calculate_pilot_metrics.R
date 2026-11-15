# Calculate performance metrics for the compcodeR pilot data.
# The DE genes labeled by compcodeR are treated as the known truth.

# Read the simulated compcodeR dataset.
pilot_data <- readRDS("data/compcodeR_pilot_n7_repl1.rds")

# Read the complete result table from each method.
deseq2_results <- read.csv("results/pilot/deseq2_pilot_results.csv",
  row.names = 1,
  check.names = FALSE )

edger_results <- read.csv("results/pilot/edger_pilot_results.csv",
  row.names = 1,
  check.names = FALSE)

limma_results <- read.csv("results/pilot/limma_voom_pilot_results.csv",
  row.names = 1,
  check.names = FALSE )

# Identify the genes that compcodeR labeled as truly DE.
# pilot_data is a compData object.
# The @ symbol opens a named slot in that object.
# variable.annotations stores information about each gene.
# The $ symbol selects the differential.expression column.
# A value of 1 means DE, and a value of 0 means not DE.

true_de <- pilot_data@variable.annotations$differential.expression == 1

# Use the count-matrix row names to get the gene names.
true_de_genes <- rownames(pilot_data@count.matrix)[true_de]

# Display the number of known DE genes.
length(true_de_genes)

# Display the first few known DE genes.
head(true_de_genes)

# NA means that an adjusted p-value was not available.
# The is.na() check prevents missing values from being included.

# Identify genes called significant by DESeq2.
deseq2_called <- rownames(deseq2_results)[
  !is.na(deseq2_results$padj) & deseq2_results$padj < 0.05]

# Identify genes called significant by edgeR.
edger_called <- rownames(edger_results)[
  !is.na(edger_results$FDR) & edger_results$FDR < 0.05]

# Identify genes called significant by limma-voom.
limma_called <- rownames(limma_results)[
  !is.na(limma_results$adj.P.Val) & limma_results$adj.P.Val < 0.05]

# Check how many genes each method called significant.
length(deseq2_called)
length(edger_called)
length(limma_called)

# --------------------------------
# DESeq2 performance metrics
# --------------------------------

# List all genes in the simulated dataset.
all_genes <- rownames(pilot_data@count.matrix)

# Compare DESeq2's results with the known compcodeR truth.

# True positives:
# DESeq2 called the gene significant, and compcodeR says it is truly DE.
true_positive <- sum(deseq2_called %in% true_de_genes)

# False positives:
# DESeq2 called the gene significant, but compcodeR says it is not truly DE.
false_positive <- sum(!(deseq2_called %in% true_de_genes))

# False negatives:
# compcodeR says the gene is truly DE, but DESeq2 did not call it significant.
false_negative <- sum(!(true_de_genes %in% deseq2_called))

# True negatives:
# These are genes that are not truly DE and were not called significant.
true_negative <- 
  length(all_genes) - 
  true_positive - 
  false_positive -
  false_negative

# Calculate precision.
precision <- true_positive / (true_positive + false_positive)

# Calculate recall.
recall <- true_positive / (true_positive + false_negative)

# Calculate the observed false discovery rate.
observed_fdr <- false_positive / (true_positive + false_positive)

# Display the DESeq2 counts and metrics.
deseq2_metrics <- data.frame(
  method = "DESeq2",
  true_positive = true_positive,
  false_positive = false_positive,
  false_negative = false_negative,
  true_negative = true_negative,
  precision = precision,
  recall = recall,
  observed_fdr = observed_fdr )

deseq2_metrics

# --------------------------------
# edgeR performance metrics
# --------------------------------

# True positives:
# edgeR called the gene significant, and compcodeR says it is truly DE.
edger_true_positive <- sum(edger_called %in% true_de_genes)

# False positives:
# edgeR called the gene significant, but compcodeR says it is not truly DE.
edger_false_positive <- sum(!(edger_called %in% true_de_genes))

# False negatives:
# compcodeR says the gene is truly DE, but edgeR did not call it significant.
edger_false_negative <- sum(!(true_de_genes %in% edger_called))

# True negatives:
# These genes are not truly DE and were not called significant by edgeR.
edger_true_negative <-
  length(all_genes) -
  edger_true_positive -
  edger_false_positive -
  edger_false_negative

# Calculate precision.
edger_precision <- edger_true_positive /
  (edger_true_positive + edger_false_positive)

# Calculate recall.
edger_recall <- edger_true_positive /
  (edger_true_positive + edger_false_negative)

# Calculate the observed false discovery rate.
edger_observed_fdr <- edger_false_positive /
  (edger_true_positive + edger_false_positive)

# Display the edgeR counts and metrics.
edger_metrics <- data.frame(
  method = "edgeR",
  true_positive = edger_true_positive,
  false_positive = edger_false_positive,
  false_negative = edger_false_negative,
  true_negative = edger_true_negative,
  precision = edger_precision,
  recall = edger_recall,
  observed_fdr = edger_observed_fdr
)

edger_metrics

# --------------------------------
# limma-voom performance metrics
# --------------------------------

# True positives:
# limma-voom called the gene significant, and compcodeR says it is truly DE.
limma_true_positive <- sum(limma_called %in% true_de_genes)

# False positives:
# limma-voom called the gene significant, but compcodeR says it is not truly DE.
limma_false_positive <- sum(!(limma_called %in% true_de_genes))

# False negatives:
# compcodeR says the gene is truly DE, but limma-voom did not call it significant.
limma_false_negative <- sum(!(true_de_genes %in% limma_called))

# True negatives:
# These genes are not truly DE and were not called significant by limma-voom.
limma_true_negative <-
  length(all_genes) -
  limma_true_positive -
  limma_false_positive -
  limma_false_negative

# Calculate precision.
limma_precision <- limma_true_positive /
  (limma_true_positive + limma_false_positive)

# Calculate recall.
limma_recall <- limma_true_positive /
  (limma_true_positive + limma_false_negative)

# Calculate the observed false discovery rate.
limma_observed_fdr <- limma_false_positive /
  (limma_true_positive + limma_false_positive)

# Display the limma-voom counts and metrics.
limma_metrics <- data.frame(
  method = "limma-voom",
  true_positive = limma_true_positive,
  false_positive = limma_false_positive,
  false_negative = limma_false_negative,
  true_negative = limma_true_negative,
  precision = limma_precision,
  recall = limma_recall,
  observed_fdr = limma_observed_fdr
)

limma_metrics

# --------------------------------
# Combine the pilot metrics
# --------------------------------

# Combine the three method-specific metric tables.
pilot_metrics <- rbind(
  deseq2_metrics,
  edger_metrics,
  limma_metrics)

# Display the comparison.
pilot_metrics

# Save the combined metrics table.
write.csv(
  pilot_metrics,
  "results/pilot/pilot_metrics.csv",
  row.names = FALSE)