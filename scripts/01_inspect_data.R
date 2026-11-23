# Check the project folder R is using.
getwd()

# Check that the count file is in the data folder.
count_file <- "data/GSE196186_Raw_gene_counts_matrix.tsv"
file.exists(count_file)

# Read the tab-separated file.
# Skip the first line because it contains processing information.
raw_counts <- read.delim(
  count_file,
  skip = 1,
  check.names = FALSE )

# Check the number of rows and columns.
dim(raw_counts)

# Check the column names.
names(raw_counts)

# Display a few gene counts from four samples.
head(raw_counts[, c("Geneid", "L1", "L2", "N1", "N2")])

# List the 14 study samples in the order used for analysis.
sample_names <- c(
  "L1", "L2", "L3", "L4", "L5", "L6", "L7",
  "N1", "N2", "N3", "N4", "N5", "N6", "N7"
)

# Select those sample columns and convert them to a numeric matrix.
counts <- as.matrix(raw_counts[, sample_names])

# Use gene IDs as the row names of the count matrix.
rownames(counts) <- raw_counts$Geneid

# Check the count matrix dimensions.
dim(counts)

# Check that the matrix has no missing or negative values
# and that the counts are whole numbers.
anyNA(counts)
all(counts >= 0)
all(counts == floor(counts))

# Make a table connecting each sample to its group.
# Non-lymphoma means the comparison group; it does not necessarily mean healthy.
sample_info <- data.frame(
  sample = colnames(counts),
  group = factor(
    c(rep("lymphoma", 7), rep("non_lymphoma", 7)),
    levels = c("non_lymphoma", "lymphoma"))
)

# Use sample names to label the rows of the sample-information table.
rownames(sample_info) <- sample_info$sample

# Check the sample information and group sizes.
print(sample_info)
table(sample_info$group)

# Confirm that the sample names match and are in the same order.
identical(colnames(counts), rownames(sample_info))