rm(list=ls())
library(adegenet)

#################################################
# GRAB DATA AND GENERATE GENIND OBJECT
#################################################
# genotype_matrix <- read.table("SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final.tab", header = TRUE, sep = "\t") # nolint
# genotype_data <- genotype_matrix[, 5:ncol(genotype_matrix)]

# # IMPORTANT!: Need to combine collumns so you have one locus in each two column and two alleles # nolint
# combined_df <- data.frame(matrix(nrow = nrow(genotype_data), ncol = 0))
# for (i in seq(1, ncol(genotype_data), by = 2)) {
#   col_name <- gsub("\\.", "_", names(genotype_data)[i])
#   combined_df[col_name] <- paste(genotype_data[[i]], genotype_data[[i + 1]], sep = "")
# }
# # Print the combined data frame
# print(combined_df)

# # Assuming 'genotype_data' is your data frame
# combined_df[] <- lapply(combined_df, function(x) ifelse(x == "00", NA, x))

# # Generate genind object
# ## IMPORTANT: `ncode=1` since each column has two alleles, but only one is the actual reference
# datagen <- df2genind(combined_df, ploidy = 2, ncode = 1, ind.names = NULL, loc.names=NULL, check.ploidy=TRUE) 


#################################################
# Mantel Test
#################################################
# From: https://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf
# The mantel test is useful because it doesn't just focus on two or a few principle components (like Jombart et al. 2008's sPCA adjacent global and local tests), which tests spatial structures in the whole data by assessing correlation betwen genetic distances and geographic distances.
library(adegenet)

# Get matrix of scaled allele frequenceis. 'scaleGen' returns a matrix of scaled allele frequencies with genotypes (genind) or populations in (genpop) in rows and alleles in columns. Input is a genind object
datagen.X <- scaleGen(datagen) 

# Calculate distance matrix of coordinates using geosphere package
library(geosphere)
coo <- genotype_matrix [,c(1,2,3)]
coo_dist <- distm(coo[,c("Longitude", "Latitude")]) # Note that order is this line is Longitude, then Latitude (unlike input dataset)
rownames(coo_dist) <- coo$Puma_ID
colnames(coo_dist) <- coo$Puma_ID

# Permform mantel test
mtest <- mantel.randtest(dist(datagen.X), dist(coo_dist)) # pairwise Euclidean distances are computed using 'dist'
pdf("EDA/SCRIPTS/Rscripts/Routput/Mantel-plot.pdf")
plot(mtest) # figure shows a histogram of permuted test statistics and indicates the observed statistics by a black dot and segment. Basically, in line/do is quite a bit higher, shows that the observed statistics is larger than most simulated values. This means you can reject the null hypothesis of absence of spatial structure, which would be encouraging for looking for spatial structure.
dev.off()


#################################################
# SPATIAL PRINCIPAL COMPONENT ANALYSIS (sPCA)
#################################################
library(adegenet)
library(pegas)
# may also need package "spdep" installed

########## SPATIAL PRINCIPAL COMPONENT ANALYSIS (sPCA) ##########
genotype_matrix <- read.table("SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final.tab", header = TRUE, sep = "\t")
xy <- genotype_matrix[,c(2,3)] # make sure you actually have this from the Generate-Genind.R script

# Need spdep package
#install.packages("spdep")
library(spdep)

mySpca <- spca(datagen, xy=xy, type=6, k=10) # 'type' refers to connection network. 'type=6' is K nearest neighbors with k=10 neighbors
# Graphical display of spca results
# Note that up to 3 axes can be chosen

# "Positive eigenvalues (on the left) correspond to global structures, while negative eigenvalues (on the right) indicate local patterns"
pdf("EDA/SCRIPTS/Rscripts/Routput/sPCA-eigenvalues.pdf")
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,25)))
dev.off()

# From Adegenet sPCA tutorial: "A structure with a low spatial autocorrelation can barely be interpreted as a spatial pattern... A structure with a low variance would likely not reflect any genetic structure"
pdf("EDA/SCRIPTS/Rscripts/Routput/sPCA-screeplot.pdf")
screeplot(mySpca)
dev.off()

# 6 plots that summarize a lot (see page 18 of adegenet spca tutorial for descriptions):
pdf("EDA/SCRIPTS/Rscripts/Routput/sPCA-6plot-summary.pdf")
plot(mySpca)
dev.off()
# "Here, [plots 2, 3, and 4] show that genotypes are split into two genetic clusters, one in the west[ish] (or left) and one in the east[ish] (right)"
# Sparse panthers/data in the SW may be contributing to NE skew


########## GLOBAL AND LOCAL TESTS ##########
# Note: global patterns are more indicative of neutral processes, whereas local patterns may be more indicative of selective (or other non-neutral) processes.
  # "The global and local tests proposed in [1] can be used to reinforce the decision of interpreting or not interpreting global and local structures"

# Global test (testing for effects of neutral processes)
myGtest <- global.rtest(mySpca$tab,mySpca$lw,nperm=9999) # note: not sure where the "obj$tab" came from in the tutorial. Takes a couple of mins to run
myGtest
pdf("EDA/SCRIPTS/Rscripts/Routput/sPCA-Gtest.pdf")
plot(myGtest)
dev.off()

# Local test (testing for effects of non-neutral processes)
myLtest <- local.rtest(mySpca$tab,mySpca$lw,nperm=9999)
myLtest
pdf("EDA/SCRIPTS/Rscripts/Routput/sPCA-Ltest.pdf")
plot(myLtest)
dev.off()
















##########################
# SANDBOX
##########################
# Basic PCA test
# # Not really a lot of useful info


# # Convert genind to genlight
# library(dartR)
# datagl <- gi2gl(datagen)
# 
# pca1 <- glPca(datagl, nf=25) # 'nf' is number of axes to maintain
# 
# # Make PCA plot
# plot(pca1$scores[,1], pca1$scores[,2], 
#      cex=2, pch=20, col=datagl$pop, 
#      xlab="Principal Component 1", 
#      ylab="Principal Component 2", 
#      main="PCA on SNP data")
# legend("topleft", 
#        legend=unique(datagl$pop), 
#        pch=20, 
#        col=c("black", "red"))
# 




