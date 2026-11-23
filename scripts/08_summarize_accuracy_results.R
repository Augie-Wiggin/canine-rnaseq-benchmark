# Summarize and plot the accuracy results from Script 07.
# This script reads the saved results; it does not rerun the benchmark.

library(ggplot2)

# Read Script 07's results table.
metrics <- read.csv("results/accuracy/accuracy_benchmark_metrics.csv")

# Check the number of rows and columns.
dim(metrics)

# Check how many results there are for each sample size and method.
table(metrics$total_samples, metrics$method)

# Calculate the mean results for each sample size and method.
accuracy_mean <- aggregate(
  cbind(
    true_positive, false_positive, false_negative, true_negative,
    precision, recall, specificity, observed_fdr
  ) ~ total_samples + method, data = metrics, FUN = mean )

# Calculate the standard deviation across the five replicates.
accuracy_sd <- aggregate(
  cbind(
    true_positive, false_positive, false_negative, true_negative,
    precision, recall, specificity, observed_fdr
  ) ~ total_samples + method, data = metrics, FUN = sd )

# Put the means and standard deviations into one table.
accuracy_summary <- merge(
  accuracy_mean,
  accuracy_sd,
  by = c("total_samples", "method"),
  suffixes = c("_mean", "_sd"))

# Display the summary table.
accuracy_summary

# Save the summary table as a CSV file.
write.csv(accuracy_summary, 
          "results/accuracy/accuracy_summary.csv",
          row.names = FALSE)

# ----------------------
# Plot the Mean Precision 
# ----------------------

# Put the sample sizes in order on the x-axis.
precision_plot_data <- accuracy_summary
precision_plot_data$total_samples <- factor(
  precision_plot_data$total_samples,
  levels = c(14, 50, 100))

# Use the same colors as the runtime plots.
method_colors <- c(
  "DESeq2" = "red",
  "edgeR" = "green",
  "limma-voom" = "blue")

# Plot mean precision for each method and sample size.
precision_plot <- ggplot(
  precision_plot_data,
  aes(x = total_samples, y = precision_mean, color = method)
) +
  geom_errorbar(
    aes(
      ymin = precision_mean - precision_sd,
      ymax = precision_mean + precision_sd),
    position = position_dodge(width = 0.5),
    width = 0.15) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3) +
  scale_color_manual(values = method_colors) +
  coord_cartesian(ylim = c(0.75, 1)) +
  labs(
    title = "Mean Precision by Total Sample Size",
    x = "Total Samples",
    y = "Mean Precision",
    color = "Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(precision_plot)

ggsave(
  "figures/accuracy_precision_by_sample_size.png",
  plot = precision_plot,
  width = 8,
  height = 6)

# -------------------
# Plot the Mean Recall
# -------------------

# Put the sample sizes in order on the x-axis.
recall_plot_data <- accuracy_summary
recall_plot_data$total_samples <- factor(
  recall_plot_data$total_samples,
  levels = c(14, 50, 100))

# Plot mean recall for each method and sample size.
recall_plot <- ggplot(
  recall_plot_data,
  aes(x = total_samples, y = recall_mean, color = method)
) +
  geom_errorbar(
    aes(
      ymin = recall_mean - recall_sd,
      ymax = recall_mean + recall_sd),
    position = position_dodge(width = 0.5),
    width = 0.15) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3) +
  scale_color_manual(values = method_colors) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(
    title = "Mean Recall (Sensitivity) by Total Sample Size",
    x = "Total Samples",
    y = "Mean Recall (Sensitivity)",
    color = "Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(recall_plot)

ggsave(
  "figures/accuracy_recall_by_sample_size.png",
  plot = recall_plot,
  width = 8,
  height = 6)
