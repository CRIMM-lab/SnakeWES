# TEST filtering algorith

## load sample bed 
library(data.table)
library(tidyr)
library(stringr)
library(dplyr)

# Load and parsing variants bed file
#system("awk 'BEGIN{FS=OFS=\"\t\"} {for(i=1;i<=NF;i++) if($i==\".\") $i=\"NA\"}1' multisample.varscan.vep.tsv > temp.tsv && mv temp.tsv multisample.varscan.vep.tsv")
df <- fread("multisample.varscan.vep.tsv", header = FALSE) %>% 
  select(-V45:-V54)

# Define a function to replace missing point values with "NA"
#replace_period_with_NA <- function(column) {
#  ifelse(column == ".", "NA", column)
#}
#df <- as.data.frame(apply(df, 2, replace_period_with_NA))
cmd <- "bcftools query -l multisample.varscan.vep.vcf.gz"
samples <- system(cmd, inter = TRUE)
sampleAD <- paste(samples, "_AD") 
sampleAF <-  paste(samples, "_AF") 
header <- c("CHROM","POS","REF","ALT","sampleID","No.samples","Consequence","IMPACT","SYMBOL","Gene","EXON","HGVSc","Amino_acids","Codons","Existing_variation","SIFT","PolyPhen","DOMAINS","EUR_AF","gnomADe_AF","gnomADe_NFE_AF","gnomADg_AF","gnomADg_NFE_AF","CLIN_SIG",sampleAD,sampleAF)

# Define population af columns
popAF_cols <- c("EUR_AF","gnomADe_AF","gnomADe_NFE_AF","gnomADg_AF","gnomADg_NFE_AF")

colnames(df) <- header
# Split vep Existing_variant column 
df <- separate(df, Existing_variation, into = c("dbSNP", "COSMIC"), sep = "&", extra = "merge") %>%
  replace_na(replace = list(dbSNP = "NA", COSMIC = "NA"))

# Step 1a - label ad germline if present in at least one pop db (gnomAD. 1kg phase 3) with a MAF > 0.01
is_germline <- apply(df[,popAF_cols], 1, function(row) any(!is.na(row) & row >= 0.01))
df$label <- ifelse(is_germline, "germline", "")

df <- filter(df, label != "germline")

# Step 1b - label ad germline if present in two or more pop db (gnomAD. 1kg phase 3) with a MAF > 0.0001
# Check if the is_germline_lowfreq condition is met
is_germlineLowfreq <- rowSums(df[, popAF_cols] >= 0.001, na.rm = TRUE) >= 2

# Create a function to generate the variant flag based on COSV and clinvar IDs
get_germline_flag <- function(row) {
  
  if (is_germlineLowfreq[row]) {

    flag <- "germlineLowFreq"
    
    if (df$COSMIC[row] != "NA" && df$CLIN_SIG[row] != "") {

      flag <- paste(flag, df$COSMIC[row], df$CLIN_SIG[row], sep = "_")
    
    } else if (df$COSMIC[row] != "NA") {
    
      flag <- paste(flag, df$COSMIC[row], sep = "_")
    
    } else if (df$CLIN_SIG[row] != "") {
    
      flag <- paste(flag, df$CLIN_SIG[row], sep = "_")
    
    }
  
  } else {
  
    flag <- df$label[row]  
  }
  
  return(flag)
}

# Apply the function to create the variant flags
df$label <- sapply(1:nrow(df), get_germline_flag)



# Step 2a - 
#is_somatic <- rowSums(is.na(df[, popAF_cols]), na.rm = TRUE) == 5 & is_germline_lowfreq

is_somatic <- rowSums(df[, popAF_cols] >= 0.001, na.rm = TRUE) < 2

#Create a function to generate the variant flag based on COSV and clinvar IDs
get_somatic_flag <- function(row) {
  
  if (is_somatic[row]) {

    flag <- "somatic"
    
    if (df$COSMIC[row] != "NA" && df$CLIN_SIG[row] != "") {

      flag <- paste(flag, df$COSMIC[row], df$CLIN_SIG[row], sep = "_")
    
    } else if (df$COSMIC[row] != "NA") {
    
      flag <- paste(flag, df$COSMIC[row], sep = "_")
    
    } else if (df$CLIN_SIG[row] != "") {
    
      flag <- paste(flag, df$CLIN_SIG[row], sep = "_")
    
    }
  
  } else {
  
    flag <- df$label[row]
  
  }
  
  return(flag)
}

df$label <- sapply(1:nrow(df), get_somatic_flag)


fwrite(df, "multisample.varscan.vep.lab.tsv", sep = "\t", quote = FALSE)
