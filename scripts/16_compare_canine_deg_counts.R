# Compare canine differential-expression results from the three methods.
# Read the full results saved by Script 15 and count genes meeting the paper's cutoffs.

# Read the three results tables.
deseq2_results <- read.csv(
  "results/canine_data/deseq2_canine_results.csv",
  row.names = 1)

edger_results <- read.csv(
  "results/canine_data/edger_canine_results.csv",
  row.names = 1)

limma_results <- read.csv(
  "results/canine_data/limma_voom_canine_results.csv",
  row.names = 1)

# Set the FDR cutoff used for all three methods.
fdr_cutoff <- 0.05

# Set the log2 fold-change cutoff used for this comparison.
# The absolute value includes genes changing in either direction.
log2_fold_change_cutoff <- 1.5

# First, find genes below the FDR cutoff for each method.
deseq2_fdr_genes <- rownames(deseq2_results)[
  !is.na(deseq2_results$padj) &
    deseq2_results$padj < fdr_cutoff]

edger_fdr_genes <- rownames(edger_results)[
  !is.na(edger_results$FDR) &
    edger_results$FDR < fdr_cutoff]

limma_fdr_genes <- rownames(limma_results)[
  !is.na(limma_results$adj.P.Val) &
    limma_results$adj.P.Val < fdr_cutoff]

# Find genes with FDR < 0.05 and an absolute log2 fold change of at least 1.5.
deseq2_log2fc_genes <- rownames(deseq2_results)[
  !is.na(deseq2_results$padj) &
    deseq2_results$padj < fdr_cutoff &
    abs(deseq2_results$log2FoldChange) >= log2_fold_change_cutoff]

edger_log2fc_genes <- rownames(edger_results)[
  !is.na(edger_results$FDR) &
    edger_results$FDR < fdr_cutoff &
    abs(edger_results$logFC) >= log2_fold_change_cutoff]

limma_log2fc_genes <- rownames(limma_results)[
  !is.na(limma_results$adj.P.Val) &
    limma_results$adj.P.Val < fdr_cutoff &
    abs(limma_results$logFC) >= log2_fold_change_cutoff]

# Count genes meeting both cutoffs for each method.
canine_summary <- data.frame(
  method = c("DESeq2", "edgeR", "limma-voom"),
  genes_meeting_paper_cutoffs = c(
    length(deseq2_log2fc_genes),
    length(edger_log2fc_genes),
    length(limma_log2fc_genes)))

# Formatting column header
names(canine_summary) <- c(
  "Method",
  "DEGs (FDR < 0.05 and |log2FC| >= 1.5)")

# Display and save the DEG counts.
print(canine_summary)

write.csv(
  canine_summary,
  "results/canine_data/canine_method_comparison.csv",
  row.names = FALSE)