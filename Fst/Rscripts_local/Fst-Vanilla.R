##################################################################################
# GET DATA, GENERATE MERGED GENOTYPE DATA AND COMBINE PATCH+GENOTYPE DATA 
# Used for all patch types
##################################################################################
rm(list=ls())
library(adegenet)
library(hierfstat)
library(pegas)
library(ggplot2)
library(dplyr)
patch_data <- read.csv("GIS/DATA/Created/Puma_patches.csv", header=T)
genotype_data <- read.table("SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final_HWEonly.tab", header = TRUE, sep = "\t")

# Order by puma_ID for consistency across datasets
# ------------------------------------------------------------------------------------------------------------
patch_data <- patch_data %>%
  arrange(Puma_ID)
genotype_data <- genotype_data %>%
  arrange(Puma_ID)
# ------------------------------------------------------------------------------------------------------------

# IMPORTANT!: Need to combine collumns so you have one locus in each two column and two alleles
genotype_data_NN <- genotype_data[, c("Puma_ID", "Latitude", "Longitude", "Datset")]
for (i in seq(5, ncol(genotype_data), by = 2)) { # Start from the 5th column to skip Puma_ID, Latitude, Longitude, and Datset
  col_name <- gsub("\\.", "_", names(genotype_data)[i]) # create  ew column name based on allele names
  genotype_data_NN[col_name] <- paste(genotype_data[[i]], genotype_data[[i + 1]], sep = "")
}

# Combine condensed/merged genotype data and patch data
data <- cbind(genotype_data_NN, patch_data[, c("rivers", "roads", "US_L3CODE", "US_L3NAME")])

# The following blocks are used to generate observed Fst data for ecoregions, roads, rivers, etc. as needed

##################################################################################
# LEVEL III ECOREGIONS
##################################################################################
# PREP DATA FOR PAIRWISE COMPARISONS
# ------------------------------------------------------------------------------------------------------------
# SNIPPET FOR ECOREGIONS (ADDS COLUMN WITH ECOREGION REPRESENTED BY NUMBER)
# Update the US_L3NAME for Puma_IDs UA00046759 and UA00013730
data <- data %>%
  mutate(US_L3NAME = ifelse(Puma_ID == "UA00046759", "Madrean Archipelago", US_L3NAME))
data <- data %>%
  mutate(US_L3NAME = ifelse(Puma_ID == "UA00013730", "Arizona/New Mexico Plateau", US_L3NAME))
# Sort and get unique values from 'US_L3NAME'
US_L3NAME_num <- sort(unique(data$US_L3NAME))
# Convert 'US_L3NAME' to a factor and then to numeric values
data$US_L3NAME_num = as.numeric(factor(data$US_L3NAME, levels = US_L3NAME_num))
# ------------------------------------------------------------------------------------------------------------

# SELECT WHAT PATCH TYPE YOUR USING TO GENERATE PERMUTATIONS
patch <- sort(unique(data$US_L3NAME_num)) # change to what you're setting as boundary

# Repeats each number by a vector (max(regs):1)+1, which repeats each number in reverse order. Unlist() converts the matrix into a vector.
a <- rep(1:max(patch), times = (max(patch):1))
a <- unlist(a) 

# use sapply() to apply function function(x) x:max(regs) to each element (x) of regs. Use unlist to turn into vector
b <- sapply(patch, function(x) x:max(patch))
b <- unlist(b)

# pair a and b
pair <- data.frame(a,b)
pair <- subset(pair, pair[,1] != pair[,2]) # don't want to compare same region to itself (e.g. 1v1, 2v2, etc.) if not sub-setting like before
row.names(pair) <- NULL
rownames(pair) <- 1:nrow(pair)

# make df for output data
outputdf <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(outputdf) <- c("Fst", "Patch")

# GENERATE VANILLA FST
# Note: It looks similar to permutation code minus a few steps (e.g. not iterating 1000+ times)
#results <- NULL # create empty df to store results if needed
results <- data.frame(Fst=numeric(), pair=character())

for (i in 1:nrow(pair)) {
  d1 <- subset(data, US_L3NAME_num %in% pair[i,])
  d2 <- d1[,!names(data) %in% c("Puma_ID", "Latitude", "Longitude","Datset","rivers","roads","US_L3CODE","US_L3NAME")] # removes unnecessary columns
  d2 <- d2[, c(ncol(d2), 1:(ncol(d2)-1))] # moving column from end of df to beginning
  dg1 <- df2genind(d2[,-1], ploidy=2, ncode=1, pop=d2$US_L3NAME_num) # IMPORTANT: Need to make sure first column has population info and that the 'df2genind' object ignores the first column (i.e., d2[,-1]) and instead specifies population separately (i.e., pop=d2$eco3)
  matFst <- genet.dist(dg1, method = "WC84") 
  output <- data.frame(Fst=matFst[1], pair=pair[i,])
  results <- rbind(results, output) # store the result in the df
}

results$paircombo <- paste(results$pair.a, results$pair.b, sep="") # create new column that shows which two regions were compaired (useful for plotting)

# CSV OF FST DATA
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-US_L3NAME_num.csv")

##################################################################################
# RIVERS
##################################################################################
# PREP DATA FOR PAIRWISE COMPARISONS
# ------------------------------------------------------------------------------------------------------------
# manually assign pumas outisde of state boundary to river region 2
data <- data %>%
  mutate(rivers = ifelse(Puma_ID == "UA00031261", 2, rivers))
data <- data %>%
  mutate(rivers = ifelse(Puma_ID == "UA00047644", 2, rivers))
# ------------------------------------------------------------------------------------------------------------

# SELECT WHAT PATCH TYPE YOUR USING TO GENERATE PERMUTATIONS
patch <- sort(unique(data$rivers)) # change to what you're setting as boundary

# Repeats each number by a vector (max(regs):1)+1, which repeats each number in reverse order. Unlist() converts the matrix into a vector.
a <- rep(1:max(patch), times = (max(patch):1))
a <- unlist(a) 

# use sapply() to apply function function(x) x:max(regs) to each element (x) of regs. Use unlist to turn into vector
b <- sapply(patch, function(x) x:max(patch))
b <- unlist(b)

# pair a and b
pair <- data.frame(a,b)
pair <- subset(pair, pair[,1] != pair[,2]) # don't want to compare same region to itself (e.g. 1v1, 2v2, etc.) if not sub-setting like before
row.names(pair) <- NULL
rownames(pair) <- 1:nrow(pair)

# make df for output data
outputdf <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(outputdf) <- c("Fst", "Patch")

# GENERATE VANILLA FST
# Note: It looks similar to permutation code minus a few steps (e.g. not iterating 1000+ times)
#results <- NULL # create empty df to store results if needed
results <- data.frame(Fst=numeric(), pair=character())

for (i in 1:nrow(pair)) {
  d1 <- subset(data, rivers %in% pair[i,])
  d2 <- d1[,!names(data) %in% c("Puma_ID", "Latitude", "Longitude","Datset","roads","US_L3CODE","US_L3NAME", "US_L3NAME_num")] # removes unnecessary columns
  d2 <- d2[, c(ncol(d2), 1:(ncol(d2)-1))] # moving column from end of df to beginning
  dg1 <- df2genind(d2[,-1], ploidy=2, ncode=1, pop=d2$rivers) # IMPORTANT: Need to make sure first column has population info and that the 'df2genind' object ignores the first column (i.e., d2[,-1]) and instead specifies population separately (i.e., pop=d2$eco3)
  matFst <- genet.dist(dg1, method = "WC84") 
  output <- data.frame(Fst=matFst[1], pair=pair[i,])
  results <- rbind(results, output) # store the result in the df
}

results$paircombo <- paste(results$pair.a, results$pair.b, sep="") # create new column that shows which two regions were compaired (useful for plotting)

# CSV OF FST DATA
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-rivers.csv")

##################################################################################
# ROADS
##################################################################################
# PREP DATA FOR PAIRWISE COMPARISONS
# ------------------------------------------------------------------------------------------------------------
# manually assign pumas outisde of state boundary to road region 1
data <- data %>%
  mutate(roads = ifelse(Puma_ID == "UA00031261", 1, roads))
data <- data %>%
  mutate(roads = ifelse(Puma_ID == "UA00047644", 1, roads))
# ------------------------------------------------------------------------------------------------------------

# SELECT WHAT PATCH TYPE YOUR USING TO GENERATE PERMUTATIONS
patch <- sort(unique(data$roads)) # change to what you're setting as boundary

# Repeats each number by a vector (max(regs):1)+1, which repeats each number in reverse order. Unlist() converts the matrix into a vector.
a <- rep(1:max(patch), times = (max(patch):1))
a <- unlist(a) 

# use sapply() to apply function function(x) x:max(regs) to each element (x) of regs. Use unlist to turn into vector
b <- sapply(patch, function(x) x:max(patch))
b <- unlist(b)

# pair a and b
pair <- data.frame(a,b)
pair <- subset(pair, pair[,1] != pair[,2]) # don't want to compare same region to itself (e.g. 1v1, 2v2, etc.) if not sub-setting like before
row.names(pair) <- NULL
rownames(pair) <- 1:nrow(pair)

# make df for output data
outputdf <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(outputdf) <- c("Fst", "Patch")

# GENERATE VANILLA FST
# Note: It looks similar to permutation code minus a few steps (e.g. not iterating 1000+ times)
#results <- NULL # create empty df to store results if needed
results <- data.frame(Fst=numeric(), pair=character())

for (i in 1:nrow(pair)) {
  d1 <- subset(data, roads %in% pair[i,])
  d2 <- d1[,!names(data) %in% c("Puma_ID", "Latitude", "Longitude","Datset","rivers","US_L3CODE","US_L3NAME", "US_L3NAME_num")] # removes unnecessary columns
  d2 <- d2[, c(ncol(d2), 1:(ncol(d2)-1))] # moving column from end of df to beginning
  dg1 <- df2genind(d2[,-1], ploidy=2, ncode=1, pop=d2$roads) # IMPORTANT: Need to make sure first column has population info and that the 'df2genind' object ignores the first column (i.e., d2[,-1]) and instead specifies population separately (i.e., pop=d2$eco3)
  matFst <- genet.dist(dg1, method = "WC84") 
  output <- data.frame(Fst=matFst[1], pair=pair[i,])
  results <- rbind(results, output) # store the result in the df
}

results$paircombo <- paste(results$pair.a, results$pair.b, sep="") # create new column that shows which two regions were compaired (useful for plotting)

# CSV OF FST DATA
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-roads.csv")

