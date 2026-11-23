# Test the runtime and memory use of DESeq2
# on one saved runtime dataset.
library(DESeq2)

# Read the dataset name passed in by Script 13.
script_arguments <- commandArgs(trailingOnly = TRUE)

# Take the first argument as the dataset name.
dataset_name <- script_arguments[1]

dataset_file <- paste0("data/runtime/", dataset_name, ".rds")

# Read the selected dataset created by Script 09.
simulation_data <- readRDS(dataset_file)

# Get the count matrix and sample information.
counts <- simulation_data@count.matrix
sample_info <- simulation_data@sample.annotations

# Give the two groups readable names.
sample_info$group <- factor(
  sample_info$condition,
  levels = c(1, 2),
  labels = c("condition_1", "condition_2"))

# Match sample information to the count-matrix columns.
rownames(sample_info) <- colnames(counts)

# Measure the time taken by the DESeq2 analysis.
deseq2_time <- system.time({
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = ~ group)
  
  dds <- DESeq(dds)
  
  deseq2_results <- results(
    dds,
    contrast = c("group", "condition_2", "condition_1"))
})

# Show the timing and a few results.
print(dataset_name)
print(deseq2_time)
print(head(deseq2_results))