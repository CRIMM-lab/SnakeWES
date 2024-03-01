# FORMATTING HEADER VEP   
args <- (commandArgs(trailingOnly = TRUE))

library(data.table)
library(readr)
library(dplyr)


vcf <- args[1]
fi_tsv <- args[2]
fo_tsv <- args[3]

cmd <- paste("bcftools query -l", vcf, sep= " ")
samples <- system(cmd, inter = TRUE)
samples <- gsub("_tumor", "", samples)
samples_depth <- paste0(rep(samples, each=3), c("_RD", "_AD", "_AF"))
samples_gt <- paste0(samples, "_GT")
samples_header <- c(samples_depth, samples_gt)

vepCSQ_header <- readLines("resources/header.vepCSQ.txt")
vepCSQ_header <- unlist(strsplit(vepCSQ_header, "\t"))

if (length(samples) > 1) {

  header <- c(vepCSQ_header[1:4], "samples_ID", "No_Samples", vepCSQ_header[5:length(vepCSQ_header)], samples_header)  

} else {

  header <- c(vepCSQ_header, samples_header)  
}

df <- fread(fi_tsv, header=FALSE) %>% 
  setNames(header) %>%
  select(-matches("Allele"))

write_tsv(df, fo_tsv)