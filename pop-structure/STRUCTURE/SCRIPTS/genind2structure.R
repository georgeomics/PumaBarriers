########## GENIND TO STRUCTURE FUNCTION: CREATE FUNCTION ##########
# From: https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R 
# Function to export to STRUCTURE format from genind object.
# genind objects are created in the R package adegenet.  The function below is an R function.
# Lindsay V. Clark, 26 July 2015
# obj: genind object
# file: file name to write
# pops: whether to include population info in the file
# Function is flexible with regards to ploidy, although genotypes are
# considered to be unambiguous.
# Missing data must be recorded as NA in obj@tab.

genind2structure <- function(obj, file="", pops=FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- adegenet::nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- adegenet::locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                           dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}



########## GENIND TO STRUCTURE FUNCTION: GENERATE GENIND OBJECT ##########
library(adegenet)
getwd() 

genotype_matrix <- read.table("SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final_HWEonly.tab", header = TRUE, sep = "\t") # nolint
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

# NOTE: THE FOLLOWING ADDED TO THIS  `genind2structure.R` SCRIPT SPECIFICALLY
# Get puma IDs/names
#names <- genotype_matrix$Puma_ID
rownames(combined_df) <- genotype_matrix$Puma_ID 

# Generate genind object
## IMPORTANT: `ncode=1` since each column has two alleles, but only one is the actual reference
datagen <- df2genind(combined_df, ploidy = 2, ncode = 1, ind.names = NULL, loc.names=NULL, check.ploidy=TRUE) 



########## GENIND TO STRUCTURE FUNCTION: RUN FUNCTION ##########
# example use: 
# data(nancycats)
# genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)
genind2structure(datagen, file="pop-structure/STRUCTURE/SCRIPTS/Routput/Puma_SNP_genotypes_cut_AZ_final_HWEonly.str") # Note no `pops=TRUE`
# IMPORTANT: 
# Afterwards, make sure your `.str` file's header only has allele names. If you have a column name/header for individuals (e.g "ind"), populations (e.g. "pop"), etc. You may have to manually delete this after running the `genindstructure.R ` script.
