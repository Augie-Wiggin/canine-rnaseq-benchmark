# Generate the simulated datasets for the runtime and memory benchmark.
# Each dataset has 1,000 genes and 100 known DE genes.
# The total sample sizes are 10, 100, 1,000, and 10,000.
# Three datasets are generated at each sample size.
# The samples are split evenly between two groups.
# This script only generates and saves the datasets.

# Load compcodeR for simulation.
library(compcodeR)

# Create the folder for the simulated datasets if it does not exist.
if (!dir.exists("data/runtime")) {
  dir.create("data/runtime", recursive = TRUE)}

# List the total sample sizes and replicate numbers to generate.
runtime_settings <- data.frame(
  total_samples = rep(c(10, 100, 1000, 10000), each = 3),
  replicate = rep(1:3, times = 4))

# Generate and save each dataset in the settings table.
for (i in 1:nrow(runtime_settings)) {
  
  # Get this dataset's total sample size and replicate number.
  total_samples <- runtime_settings$total_samples[i]
  samples_per_group <- total_samples / 2
  replicate_id <- runtime_settings$replicate[i]
  
  # Name the dataset using its sample size and replicate number.
  dataset_name <- paste0("runtime_samples", total_samples,
                         "_repl", replicate_id)
  
  # Set a different seed for each dataset so it can be recreated.
  set.seed(total_samples * 1000 + replicate_id)
  
  # Simulate 1,000 genes with 100 known DE genes.
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
      "data/runtime/", dataset_name, ".rds") )
  
  # Check the number of genes, samples, and samples in each group.
  print(dataset_name)
  print(dim(simulation_data@count.matrix))
  print(table(simulation_data@sample.annotations$condition))
  
  # Check that compcodeR considers the dataset valid.
  print(check_compData(simulation_data))
}