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
pdf("PrelimPopgen/Routput/Ho-per-locus.pdf")
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
dev.off()

pdf("PrelimPopgen/Routput/Ho-vs-He.pdf")
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
     # Fitak et al. 2016: "Only three loci (PP05, PP08, and PP19) showed significant deviations from Hardy–Weinberg equilibrium after Bonferroni correction."
     # This is similar to what I observe for p_adjusted

##########################################
# Effective Population Size?
##########################################
# BiocManager::install("SNPRelate")
# install.packages("dartR")
## https://search.r-project.org/CRAN/refmans/dartR/html/gl.LDNe.html
library(dartR)

# Covert genind to genlight
datagl <- gi2gl(datagen, parallel = FALSE, verbose = NULL)


# # Options
# gl.LDNe(
#   x, # name of the genlight object containing SNPdata
#   outfile = "genepopLD.txt", # File name of the output file with all results from Neestimator 2 [default 'genepopLD.txt'].
#   outpath = tempdir(),
#   neest.path = getwd(), # Path to the folder of the NE2-1 (Waples and friends NeEstimator software)
#   critical = 0,
#   singleton.rm = TRUE,
#   mating = "random",
#   plot.out = TRUE,
#   plot_theme = theme_dartR(),
#   plot_colors_pop = discrete_palette,
#   save2tmp = FALSE,
#   verbose = NULL
# )

# # Example
# nes <- gl.LDNe(pops, outfile="popsLD.txt", outpath=tempdir(), 
# neest.path = "./path_to Ne-21", critical=c(0,0.05), 
# singleton.rm=TRUE, mating='random')
# nes

# Generate Ne
outpath <- getwd() # weird but only works for current working directory for this function
nes <- gl.LDNe(datagl, outfile="NeEstLD_output.txt", outpath=outpath, # not sure what the issue with the outfile is ("genepopLD.txt"), but this should be provided by the NeEstimator software
neest.path = "programs/NeEstimator/", critical=c(0,0.05), 
singleton.rm=TRUE, mating='random')
