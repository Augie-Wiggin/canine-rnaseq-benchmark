# Test the runtime and memory use of limma-voom
# on one saved runtime dataset.
library(edgeR)
library(limma)

# Read the dataset name passed in by Script 13.
script_arguments <- commandArgs(trailingOnly = TRUE)
dataset_name <- script_arguments[1]
dataset_file <- paste0("data/runtime/", dataset_name, ".rds")

# Read the selected dataset.
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

# Describe the comparison between the two groups.
design <- model.matrix(~ group, data = sample_info)

# Measure the time taken by the limma-voom analysis.
limma_time <- system.time({
  dge <- DGEList(
    counts = counts,
    group = sample_info$group)
  
  dge <- calcNormFactors(dge)
  
  voom_data <- voom(
    dge,
    design,
    plot = FALSE)
  
  limma_fit <- lmFit(voom_data, design)
  limma_fit <- eBayes(limma_fit)
  
  limma_results <- topTable(
    limma_fit,
    coef = 2,
    number = Inf,
    sort.by = "P")
})

# Show the timing and a few results.
print(dataset_name)
print(limma_time)
print(head(limma_results))