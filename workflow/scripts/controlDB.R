args <- (commandArgs(trailingOnly = TRUE))

library(PureCN)
library(R.utils)

customPoN <- "results/customPoN.PureCN.vcf.gz"
ctr_list <- "alignments/control_coverages.list"
controlDB <- "data/controlDB.rds"
mapping_bias <- "data/mapping_bias.rds"

cov_loess_list <- readLines(ctr_list)

normalDB <- createNormalDatabase(cov_loess_list)
saveRDS(normalDB, file = controlDB)

# Require a PoN with a minimum of 5 controls
bias <- calculateMappingBiasVcf("customPoN.vcf.gz", genome = "hg38", min.normals.betafit = 5, min.normals.position.specific.fit= 5)	
saveRDS(bias, file = mapping_bias)
