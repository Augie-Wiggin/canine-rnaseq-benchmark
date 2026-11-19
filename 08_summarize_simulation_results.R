# Summarize the results from the 20-dataset accuracy benchmark.
# This script reads the metrics saved by Script 07.
# It does not generate simulations or rerun the methods.

metrics <- read.csv(
  "results/simulation/simulation_benchmark_metrics.csv")

# Check the number of result rows and columns.
dim(metrics)

# Add the sample size for each dataset as a new column.
# Start by labeling every row as 7 samples per group.
metrics$samples_per_group <- 7

# Change rows whose dataset name contains "n15" to 15.
metrics$samples_per_group[grepl("n15", metrics$dataset)] <- 15

# Check how many result rows belong to each sample size.
table(metrics$samples_per_group)

# Calculate the average metrics for each method and sample size.
accuracy_summary <- aggregate(
  cbind(precision, recall, observed_fdr) ~
    samples_per_group + method,
  data = metrics,
  FUN = mean )

# Calculate the standard deviation across replicates for each method 
# and sample size.
accuracy_sd <- aggregate(
  cbind(precision, recall, observed_fdr) ~
    samples_per_group + method,
  data = metrics,
  FUN = sd)

# Match the mean and SD rows by sample size and method.
accuracy_summary_full <- merge(
  accuracy_summary,
  accuracy_sd,
  by = c("samples_per_group", "method"),
  suffixes = c("_mean", "_sd")
)

# Display the combined summary.
accuracy_summary_full

# Save the summary table.
write.csv(
  accuracy_summary_full,
  "results/simulation/simulation_accuracy_summary.csv",
  row.names = FALSE)

# Save the recall plot as an image.
png("figures/accuracy_recall_summary.png",
  width = 1200,
  height = 900
)

# Plot average recall for each method and sample size.
barplot(
  accuracy_summary_full$recall_mean,
  names.arg = paste(
    accuracy_summary_full$method,
    accuracy_summary_full$samples_per_group,
    sep = "\n"),
  ylim = c(0, 1),
  ylab = "Average recall",
  xlab = "Method and samples per group",
  main = "Average recall across simulation datasets" )

# Finish saving the plot to the PNG file.
dev.off()