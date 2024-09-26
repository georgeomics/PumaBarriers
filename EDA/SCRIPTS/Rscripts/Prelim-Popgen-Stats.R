##########################################
# Observed and Expected Heterozygosity
##########################################
# From: https://popgen.nescent.org/StartSNP.html

# REQUIRED: Generated genind object (e.g. "datagen")

# Generate summary stats including He and Ho
div <- summary(datagen)
div

# Plots
names(div)
pdf("EDA/SCRIPTS/Rscripts/Routput/Ho-per-locus.pdf")
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
dev.off()

pdf("EDA/SCRIPTS/Rscripts/Routput/Ho-vs-He.pdf")
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")
dev.off()

# Test if there is a difference in the variance of He and Ho
#bartlett.test(list(div$Hexp, div$Hobs)) # may not need this

##########################################
# Testing Hardy-Weinberg Equilibrium
##########################################
library(pegas)
library(adegenet)

# REQUIRED: Generated genind object (e.g. "datagen")
pumadata.gtm <- datagen@tab # Actual genotype matrix from genind object

# HWE
hwe_result_exact <- hw.test(datagen, B=1000) # Exact test using MCMC
hwe_result_Chi <- hw.test(datagen, B=0) # Classic Chi-square
# Bonferroni correction
p_adjusted <- p.adjust(hwe_result_exact[,4], method = "bonferroni") # can also use the following methods: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"