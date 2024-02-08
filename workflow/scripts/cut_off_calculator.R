# Load required libraries if necessary
library("data.table")

# Get the trailing (non-option) command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if at least one argument is provided
if (length(args) < 1) {
  cat("Usage: Rscript script.R sample.coverage.per-base.bed.gz")
  quit(status = 1)
}	

# Read the coverage data from a file
coverage_data <- fread(args[1])

# Calculate mean coverage
mean_coverage <- mean(coverage_data$V4)

# Calculate min cov and min supp reads
cov_cutoff <- ceiling(mean_coverage/3)
supp_reads_cutoff <- ceiling(cov_cutoff*0.03)

# Print the results
cat("cov_cutoff",cov_cutoff,"\n", sep="\t")
cat("supp_reads_cutoff",supp_reads_cutoff,"\n", sep="\t")

