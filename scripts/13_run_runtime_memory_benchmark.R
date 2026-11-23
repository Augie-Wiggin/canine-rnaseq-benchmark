# Run the runtime and memory benchmark.
# Scripts 10–12 run the three methods on the datasets made by Script 09.
# This script runs all three methods on all 12 datasets (36 runs total).
# Results are saved in results/runtime/runtime_memory_results.csv.

# Make the names of the 12 datasets created by Script 09.
# The first line repeats each sample size three times.
# The second line adds replicate numbers 1, 2, and 3 for each size.
dataset_names <- paste0(
  "runtime_samples",
  rep(c(10, 100, 1000, 10000), each = 3),
  "_repl",
  rep(1:3, times = 4))

# Make a small table that connects each method to its R script.
method_scripts <- data.frame(
  method = c("DESeq2", "edgeR", "limma-voom"),
  script = c(
    "scripts/10_test_deseq2_runtime.R",
    "scripts/11_test_edger_runtime.R",
    "scripts/12_test_limma_voom_runtime.R"))

# Print the dataset names and method scripts so they can be checked.
print(dataset_names)
print(method_scripts)

# Create the results folder if it is not already there.
if (!dir.exists("results/runtime")) {
  dir.create("results/runtime", recursive = TRUE)}

# Start an empty table. A row will be added after each method finishes.
runtime_memory_results <- data.frame()

# Go through all 12 datasets, one at a time.
for (dataset_name in dataset_names) {
  
  # For the current dataset, run each of the three methods.
  for (i in 1:nrow(method_scripts)) {
    
    # Use the current row of the method table to get its name and script.
    method_name <- method_scripts$method[i]
    method_script <- method_scripts$script[i]
    
    # Print the method and dataset currently being analyzed.
    print(paste("Running", method_name, "on", dataset_name))
    
    # Start the method script in a separate R session.
    # macOS time -l measures the total elapsed time and peak memory for
    # that session. The dataset name is passed to the method script.
    # Messages from R and the macOS time measurements are captured so the
    # time and memory numbers can be saved below.
    method_output <- system2(
      "/usr/bin/time",
      args = c("-l", "Rscript", method_script, dataset_name),
      stdout = TRUE,
      stderr = TRUE)
    
    # Check whether the method script finished successfully.
    # A nonzero status means there was an error, so stop the benchmark.
    run_status <- attr(method_output, "status")
    
    if (!is.null(run_status) && run_status != 0) {
      stop(paste("The", method_name, "run failed for", dataset_name))}
    
    # grep searches the captured output for the lines with elapsed time
    # and peak memory.
    time_line <- grep(" real ", method_output, value = TRUE)
    memory_line <- grep(
      "maximum resident set size",
      method_output, 
      value = TRUE)
    
    # Stop if either measurement is missing from the macOS output.
    if (length(time_line) == 0 || length(memory_line) == 0) {
      stop(paste("Time or memory was missing for",
                 method_name, "on", dataset_name))}
    
    # trimws removes spaces at the start of the line. sub keeps the number
    # before the first space. That number is the elapsed time in seconds.
    process_elapsed_seconds <- as.numeric(sub(" .*", "", trimws(time_line[1])))
    
    # Do the same for the memory line. The number is peak memory in bytes.
    peak_rss_bytes <- as.numeric(sub(" .*", "", trimws(memory_line[1])))
    
    # Convert bytes to megabytes using 1,000,000 bytes per MB.
    peak_memory_MB <- peak_rss_bytes / 1000000
    
    # Put the dataset name, method, time, and memory in one results row.
    method_measurement <- data.frame(
      dataset = dataset_name,
      method = method_name,
      run_time_seconds = process_elapsed_seconds,
      peak_memory_MB = peak_memory_MB )
    
    # Add this row to the table of completed method runs.
    runtime_memory_results <- rbind(runtime_memory_results, method_measurement)
    
    # Save the updated table after every method finishes.
    # If a later run fails, the completed measurements are still in this file.
    write.csv(runtime_memory_results,
              "results/runtime/runtime_memory_results.csv",
              row.names = FALSE)
    
    # Show this run's measurements and the captured macOS/R output.
    print(method_measurement)
    print(method_output)
  }
}

# Print all measurements collected by the full benchmark.
print(runtime_memory_results)
