#Rscript

library(ggplot2)
rm=(list=ls())

args <- commandArgs(trailingOnly = TRUE)

allVAF <- read.table(file.path(args[1]), header = F, col.names = "VAF")

p <- ggplot(allVAF, aes(x = VAF)) + 
  geom_density() + geom_vline(aes(xintercept = median(VAF)), linetype = "dashed")
ggsave(plot = p, filename = file.path(args[2]))
