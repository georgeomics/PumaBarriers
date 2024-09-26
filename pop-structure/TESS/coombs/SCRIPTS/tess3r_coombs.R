### Load Packages
rm(list=ls())
library(devtools)
library(pkgbuild)
library(RcppEigen)
library(tess3r)
library(adegenet)
library(data.table)

#########################
# Grab SNP data and generate a genind object (also need this block as we reference `genotype_matrix` later)
#########################
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
pumadata.gtm <- datagen@tab

#########################
# TESS3
#########################
# Grab coordinate data
coo <- genotype_matrix[,c("Latitude","Longitude")] # grabbing just lat/long
## Note: check if X/Y AKA lat long are correctly assigned
# convert to double (IMPORTANT!)
coo <- cbind(coo$Latitude, coo$Longitude) 


# running t3ssr from K=1 to 10 (may want to add 'rep' option (to rep each K 10 times like with structure, and max.iteration higher than 200 (default) to match structure))
maxk=10
tess3.obj <- tess3(X=pumadata.gtm,
                   coord=as.matrix(coo),
                   K=1:maxk, # number of ancestral pops we want to test
                   rep=10,
                   max.iteration=500000, 
                   method="projected.ls", # projection method
                   ploidy=2, # ploidy level
                   openMP.core.num=16, # how many cores for analysis
                   keep="all") # which result to keep for each value of k. If "best", only the result with the lowest rmse score will be kept for each value of K. If "all", all results will be kept and returned for each value of K.

# Save tess3.obj
save(tess3.obj, file = "pop-structure/TESS/coombs/SCRIPTS/Routput/tess3r.RData", version = 2)


