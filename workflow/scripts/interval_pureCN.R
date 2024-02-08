args <- (commandArgs(trailingOnly = TRUE))

library(PureCN)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

intervals <- import(args[1])
reference <- args[2]
mappability <- import(args[3])
res <- args[4]

outTarget <- preprocessIntervals(intervals, reference, mappability = mappability)

outTarget <- annotateTargets(outTarget, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db)
PureCN:::.writeIntervals(outTarget, res)
