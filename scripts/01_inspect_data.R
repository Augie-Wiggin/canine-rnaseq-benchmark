# Show the folder R is currently working in.
getwd()

# Check that the data file exists. TRUE means R can find it.
file.exists("data/GSE196186_Raw_gene_counts_matrix.tsv")

# Read the tab-separated file and store it as raw_counts.
# Skip the first line because it contains processing information.
# Preserve the original column names.
raw_counts <- read.delim(
  "data/GSE196186_Raw_gene_counts_matrix.tsv",
  skip = 1,
  check.names = FALSE
)

# Show the number of rows and columns. Expected: 30,951 genes and 20 columns.
dim(raw_counts)

# Show all column names.
names(raw_counts)

# Display the first six genes and counts from four samples.
head(raw_counts[, c("Geneid", "L1", "L2", "N1", "N2")])

# List the samples in the order we want:
# seven lymphoma samples, followed by seven non-lymphoma samples.
# c() combines these names into one list of values.
sample_names <- c(
  "L1", "L2", "L3", "L4", "L5", "L6", "L7",
  "N1", "N2", "N3", "N4", "N5", "N6", "N7"
)

# Select only the sample-count columns in the order listed above.
# Convert them into a matrix: a rectangular table of numbers.
counts <- as.matrix(raw_counts[, sample_names])

# Label each row with its gene ID. $ selects a named column from a table.
rownames(counts) <- raw_counts$Geneid

# Check the matrix size. Expected: 30,951 genes and 14 samples.
dim(counts)

# Check whether any counts are missing. Expected: FALSE.
anyNA(counts)

# Check that every count is zero or greater. Expected: TRUE.
all(counts >= 0)

# Check that every count is a whole number. Expected: TRUE.
# floor() rounds down; whole numbers remain unchanged.
all(counts == floor(counts))

# Create a table describing which group each sample belongs to.
# colnames() retrieves the sample names from the count matrix.
# rep() repeats each group label seven times.
# "Non-lymphoma" does not necessarily mean healthy.
sample_info <- data.frame(
  sample = colnames(counts),
  group = c(
    rep("lymphoma", 7),
    rep("non_lymphoma", 7)
  )
)

# Label the rows of the sample-information table with sample names.
rownames(sample_info) <- sample_info$sample

# Print the sample-information table.
print(sample_info)

# Count the samples in each group. Expected: seven per group.
table(sample_info$group)

# Confirm that sample names and their order match between the tables.
# Expected: TRUE. This prevents assigning counts to the wrong group.
identical(colnames(counts), rownames(sample_info))