# This script creates one pilot RNA-seq dataset with compcodeR.
# The simulation has known differentially expressed genes.
# The setup follows the official compcodeR vignette with smaller
# settings for this pilot test.

# Load compcodeR.
library(compcodeR)

# Generate one simulated dataset.
pilot_data <- generateSyntheticData(
  dataset = "pilot_n7",
  n.vars = 1000,
  samples.per.cond = 7,
  n.diffexp = 100,
  repl.id = 1,
  seqdepth = 1e7,
  fraction.upregulated = 0.5,
  between.group.diffdisp = FALSE,
  filter.threshold.total = 1,
  filter.threshold.mediancpm = 0,
  fraction.non.overdispersed = 0,
  output.file = "data/compcodeR_pilot_n7_repl1.rds"
)

# A compData object stores different parts of the simulation
# in separate slots. Show the names of those slots.
slotNames(pilot_data)

# The @ symbol opens a named slot inside the compData object.
# This line checks the number of genes and samples.
dim(pilot_data@count.matrix)

# Check the number of samples in each simulated group.
table(pilot_data@sample.annotations$condition)

# Identify the genes labeled as truly differentially expressed.
# This creates a TRUE/FALSE value for each gene.
true_de <- pilot_data@variable.annotations$differential.expression == 1

# Count the number of TRUE values.
sum(true_de)

# Display the names of the first few known DE genes.
head(rownames(pilot_data@count.matrix)[true_de])

# Display the first few rows of the count matrix.
head(pilot_data@count.matrix)

# Check whether the compData object is valid.
check_compData(pilot_data)

# Create a diagnostic report for the pilot simulation.
summarizeSyntheticDataSet(
  data.set = "data/compcodeR_pilot_n7_repl1.rds",
  output.filename = "compcodeR_pilot_n7_datacheck.html")
