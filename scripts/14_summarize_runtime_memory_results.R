# Summarize and plot the average runtime and memory results.
# This script reads the CSV made by Script 13; it does not rerun the methods.
# It averages the three repeats at each sample size and saves two plots.

library(ggplot2)

# Read the runtime and memory measurements from Script 13.
runtime_results <- read.csv("results/runtime/runtime_memory_results.csv")

# Check the number of rows and columns in the results table.
dim(runtime_results)

# Add the total sample size based on each dataset's name.
runtime_results$total_samples <- 10

runtime_results$total_samples[
  grepl("samples100_", runtime_results$dataset)] <- 100

runtime_results$total_samples[
  grepl("samples1000_", runtime_results$dataset)] <- 1000

runtime_results$total_samples[
  grepl("samples10000_", runtime_results$dataset)] <- 10000

# Check that all four sample sizes are included.
table(runtime_results$total_samples)

# Calculate the average runtime and memory for each method and sample size.
runtime_summary <- aggregate(
  cbind(run_time_seconds, peak_memory_MB) ~ total_samples + method,
  data = runtime_results,
  FUN = mean)

# Give the summary columns clear names.
names(runtime_summary) <- c(
  "total_samples",
  "method",
  "average_runtime_seconds",
  "average_memory_MB")

# Convert average memory from MB to GB.
runtime_summary$average_memory_GB <-
  runtime_summary$average_memory_MB / 1000

# Display and save the summary table.
runtime_summary

write.csv(
  runtime_summary,
  "results/runtime/runtime_memory_summary.csv",
  row.names = FALSE)

# Use the same color for each method in both plots.
method_colors <- c(
  "DESeq2" = "red",
  "edgeR" = "green",
  "limma-voom" = "blue")

# ------------------------------
# Plot Average Runtime
# ------------------------------

runtime_plot <- ggplot(
  runtime_summary,
  aes(
    x = total_samples,
    y = average_runtime_seconds,
    color = method,
    group = method)
) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = method_colors) +
  scale_x_log10(
    breaks = c(10, 100, 1000, 10000),
    labels = c("10", "100", "1,000", "10,000")) +
  labs(
    title = "Runtime by Total Sample Size",
    x = "Total Sample Size",
    y = "Average Runtime (Seconds)",
    color = "Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(runtime_plot)

ggsave(
  "figures/runtime_by_sample_size.png",
  plot = runtime_plot,
  width = 8,
  height = 6,
  dpi = 300)

# ------------------------------
# Plot Average Memory
# ------------------------------

memory_plot <- ggplot(
  runtime_summary,
  aes(
    x = total_samples,
    y = average_memory_GB,
    color = method,
    group = method)
) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = method_colors) +
  scale_x_log10(
    breaks = c(10, 100, 1000, 10000),
    labels = c("10", "100", "1,000", "10,000")) +
  labs(
    title = "Memory Usage by Total Sample Size",
    x = "Total Sample Size",
    y = "Average Memory (GB)",
    color = "Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(memory_plot)

ggsave(
  "figures/memory_usage_by_sample_size.png",
  plot = memory_plot,
  width = 8,
  height = 6,
  dpi = 300)