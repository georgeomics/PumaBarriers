#################################################
# GRAB DATA AND GENERATE GENIND OBJECT
#################################################
rm(list=ls())
library(adegenet)
getwd()

genotype_matrix <- read.table("SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final.tab", header = TRUE, sep = "\t") # nolint
genotype_data <- genotype_matrix[, 5:ncol(genotype_matrix)]

# IMPORTANT!: Need to combine collumns so you have one locus in each two column and two alleles # nolint
combined_df <- data.frame(matrix(nrow = nrow(genotype_data), ncol = 0))
for (i in seq(1, ncol(genotype_data), by = 2)) {
  col_name <- gsub("\\.", "_", names(genotype_data)[i])
  combined_df[col_name] <- paste(genotype_data[[i]], genotype_data[[i + 1]], sep = "")
}
# Print the combined data frame
print(combined_df)

# Assuming 'genotype_data' is your data frame
combined_df[] <- lapply(combined_df, function(x) ifelse(x == "00", NA, x))

# Generate genind object
## IMPORTANT: `ncode=1` since each column has two alleles, but only one is the actual reference
datagen <- df2genind(combined_df, ploidy = 2, ncode = 1, ind.names = NULL, loc.names=NULL, check.ploidy=TRUE) 


