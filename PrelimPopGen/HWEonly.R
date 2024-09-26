##################################################################################
# REMOVE SNPs DEVIATING FROM HWE (PP05, PP08, PP19)
##################################################################################
rm(list=ls())
# library(adegenet)
# library(hierfstat)
# library(pegas)
# library(ggplot2)
library(dplyr)
genotype_data <- read.table("SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final.tab", header = TRUE, sep = "\t")

# Remove SNPs deviating from HWE
genotype_data <- genotype_data %>%
  select(-PP05, -PP05.1, -PP08, -PP08.1, -PP19, -PP19.1)

# Remove the ".1" from column names
colnames(genotype_data) <- gsub("\\.1$", "", colnames(genotype_data))

# Check data
head(genotype_data)

# write file
write.table(genotype_data, file = "SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final_HWEonly.tab", sep = "\t", row.names = FALSE, quote = FALSE)

# Note: need to remove any ".1" in column names and remove quotes around zeros to be consistent with original dataset. Can do this with the following code: 
# `awk 'BEGIN {FS=OFS="\t"} NR==1 {gsub(/.1/, "", $0)} NR>1 {gsub(/"0"/, 0, $0); gsub(/".1"/, "", $0)} 1' "Puma_SNP_genotypes_cut_AZ_final_HWEonly.tab" > "Puma_SNP_genotypes_cut_AZ_final_HWEonly.tab"`

