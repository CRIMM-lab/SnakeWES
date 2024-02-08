args <- (commandArgs(trailingOnly = TRUE))

library(PureCN)
library(R.utils)

interval <- args[1]
bam <- args[2]
coverage <- args[3]
coverage_loess <- args[4]
qc <- args[5]

calculateBamCoverageByInterval(bam.file = bam, interval.file = interval, output.file = coverage)
correctCoverageBias(coverage, interval, output.file = coverage_loess, plot.bias = FALSE, output.qc.file= qc)