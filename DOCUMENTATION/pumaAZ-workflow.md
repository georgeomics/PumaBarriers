# Arizona Puma Project Pipeline
Manuscript *"Evaluating the impact of landscape features on genetic structure in pumas using a hypothesis-based Bayesian approach"* available at: https://doi.org/10.1016/j.biocon.2026.111877

# Overview
Below is a brief outline regarding the steps completed for this project:

* Data Download
  * Processing/Clean-up of SNP data
  * Adding Patch Info in GIS
* Preliminary PopGen Stats
  * Creating a Genind Object
  * Observed and Expected Heterozygosity
  * Testing Hardy-Weinberg Equilibrium
* Exploratory Data Analysis
  * Mantel Test
  * sPCA
  * Global and Local Tests of Landscape Patterns of Genetic Variation
* Estimating Population Structure
  * STRUCTURE
* Fixation Index (Fst) Analysis
  * Fst Permutations
  * Traditional Fst
  * P-Value and Z-score Estimation
* Modeling
  * Retrieving Ancestry Coefficients (Q) from STRUCTURE Output
  * R Script for Modeling
  * Stan Scripts Specifying Model and Parameters
  * Model Summaries
  * STRUCTURE-like Barplots

# Data Download
SNP data was downloaded from Fitak et al. 2016's dataset available at https://doi.pangaea.de/10.1594/PANGAEA.835154.

Landscape feature data included the following

* Data for major roads downloaded from ESRI's USA Freeway System layer package: https://hub.arcgis.com/maps/esri::usa-freeway-system.
* River data downloaded from the US Department of Agriculture's US Major Rivers layer package: https://koordinates.com/layer/12243-us-major-rivers-national/.
* A shapefile delineating level III ecoregions from the US EPA ecosystems research website: https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states.

*Note: For analyses, all geospatial data was projected to the UTM Zone 12N coordinate system (Arizona)*

* The Little Colorado River was added to the major rivers dataset as it is a major tributary of the Colorado river and forms the Little Colorado River gorge in its lower 50 miles as it approaches the Grand Canyon. A new vector line was constructed to link the sources of the Little Colorado, Salt, and San Francisco (tributary of the Gila River) rivers to delineate distinct geographic regions based on river courses in Arizona. These river sources were located within 30 km of each other near the Mount Baldy Wilderness in the Apache-Sitgreaves National Forest.
* The EPA level III ecoregions within Arizona were classified into three groups: group 1 included the Colorado Plateau, Arizona/New Mexico Plateau, and Arizona/New Mexico Mountains ecoregions; group 2 included the Sonoran Basin and Range and Mojave Basin and Range ecoregions; and group 3 consisted of the Madrean Archipelago ecoregion. These ecoregions generally define the geological structure in Arizona, including the Mogollon Rim escarpment, and aid in accounting for reduced gene flow through topographically complex areas.

## Processing/Clean-up of SNP data

* Raw downloaded `Puma_SNP_genotypes.tab` data contains many extra empty columns. Can clean up using: `cut -f1-54 Puma_SNP_genotypes.tab > Puma_SNP_genotypes_cut.tab`
* A 1 degree (~100km) buffer was created in QGIS to include points outside of the state boundaries of Arizona. Data points were extracted if the fell within this region. The resulting dataset is called `Puma_SNP_genotypes_cut_AZ.tab`
* Due to QGIS export formatting, also need to remove any "_1" in column names and remove quotes around zeros to be consistent with original dataset. Can do this with the following code: `awk 'BEGIN {FS=OFS="\t"} NR==1 {gsub(/_1/, "", $0)} NR>1 {gsub(/"0"/, 0, $0); gsub(/"_1"/, "", $0)} 1' "Puma_SNP_genotypes_cut_AZ.tab" > "Puma_SNP_genotypes_cut_AZ_final.tab"`

## Adding Patch Info in GIS

For this project, habitat patches were defined using the following geographic datasets:

* USA Major Roads
* USA Major Rivers
* EPA Level III Ecoregions

Both the puma location data as well as the shapefiles for defining patches were loaded into QGIS. Additionally, the following steps were performed to convert linear shapefiles to polygon shapefiles (using USA_Major_Rivers line vector dataset as example):

1. Add downloaded dataset (e.g. USA_Major_Rivers)
2. "Extract/clip by extent" to Arizona (Requires extracting Arizona polygon from US state boundary data)
3. Edit rivers/vertices and add additional rivers as needed
4. Little Colorado River was manually added to major rivers shapefile
   * Dissolve rivers into one line vector
5. "Merge vector layers" of AZ and rivers
6. Use "Polygonize" to create polygons based on merged line vector
7. Export to save as permanent ESRI shapefile
8. Re-add as layer and use "Extract by expression" to keep only main polygons (40,41,1,42,43)
   * Many other tiny random polygons likely came from western state boundary
9. Use attribute table editing to rename polygons according to desired number scheme
   * Use "Flash feature" right click option to see which polygon the FID corresponds to
10. Export to save as permanent ESRI shapefile

Finally, the following was performed to join puma location data with the patch (i.e. polygon) where they were located

11. Use "Join attributes by location" with point data as base layer and polygons as join layer to add polygon info to puma data
    Export file with polygon column added as CSV


# Preliminary PopGen Stats

## Creating a Genind Object

Genind objects contain individual genotypes. They also enable easy conversion between genind and other formats.

* Creating a genind object is done many times throughout these analysis. Below is the script for generating such an object from the `` dataset.

### `Generate-Genind.R`
```R
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
```

## Observed and Expected Heterozygosity

From: https://popgen.nescent.org/StartSNP.html

### `Prelim-Popgen-Stats.R` (part 1)
```R
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
```

## Testing Hardy-Weinberg Equilibrium

### `Prelim-Popgen-Stats.R` (part 2)
```R
library(pegas)
library(adegenet)

# REQUIRED: Generated genind object (e.g. "datagen")
pumadata.gtm <- datagen@tab # Actual genotype matrix from genind object

# HWE
hwe_result_exact <- hw.test(datagen, B=1000) # Exact test using MCMC
hwe_result_Chi <- hw.test(datagen, B=0) # Classic Chi-square
# Bonferroni correction
p_adjusted <- p.adjust(hwe_result_exact[,4], method = "bonferroni")
```


# Exploratory Data Analysis

## Mantel Test
The Mantel test assesses spatial structure across the entire dataset by assessing correlation between two matrices (e.g. genetic distances vs geographic distances). Calculates a correlation coefficient (often Pearson's correlation), helping determine if there is a significant relationship between the two sets of data.

* In `R`, the X-axis contains simulated r value by permutation and actual estimated r values represented by black vertical line.
* The mantel test is useful because it doesn't just focus on two or a few principle components (like Jombart et al. 2008's sPCA adjacent global and local tests), but rather tests spatial structures in the whole data by assessing correlation betwen genetic distances and geographic distances
* The data/figure shows a histogram of permuted test statistics and indicates the observed statistics by a black dot and segment. Basically, in line/do is quite a bit higher, shows that the observed statistics is larger than most simulated values. This means you can reject the null hypothesis of absence of spatial structure, which would be encouraging if you're looking for spatial structure.

### `Mantel-sPCA.R` (part 1)
```R
rm(list=ls())

# REQUIRED: Generated genind object (e.g. "datagen")

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
pdf("/home/ge199066/projects/puma-AZ/EDA/SCRIPTS/Rscripts/Routput/Mantel-plot.pdf")
plot(mtest) # figure shows a histogram of permuted test statistics and indicates the observed statistics by a black dot and segment. Basically, in line/do is quite a bit higher, shows that the observed statistics is larger than most simulated values. This means you can reject the null hypothesis of absence of spatial structure, which would be encouraging for looking for spatial structure.
dev.off()
```

## sPCA
* For the sPCA eigenvalues plot: "positive eigenvalues (on the left) correspond to global structures, while negative eigenvalues (on the right) indicate local patterns"
* For the sPCA screeplot: From Adegenet sPCA tutorial: "A structure with a low spatial autocorrelation can barely be interpreted as a spatial pattern... A structure with a low variance would likely not reflect any genetic structure"
* For the 6 plot summary: "[plots 2, 3, and 4] show that genotypes are split into two genetic clusters, one in the west[ish] (or left) and one in the east[ish] (right)". Note that sparse panthers/data in the SW may be contributing to NE skew.

### `Mantel-sPCA.R` (part 2)
```R
library(adegenet)
library(pegas)
#may also need package "spdep" installed (i.e. install.packages("spdep"))
#library(spdep)

mySpca <- spca(datagen, xy=xy, type=6, k=10) # 'type' refers to connection network. 'type=6' is K nearest neighbors with k=10 neighbors
# Graphical display of spca results
# Note that up to 3 axes can be chosen

pdf("EDA/SCRIPTS/Rscripts/Routput/sPCA-eigenvalues.pdf")
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,25)))
dev.off()

pdf("EDA/SCRIPTS/Rscripts/Routput/sPCA-screeplot.pdf")
screeplot(mySpca)
dev.off()

# 6 plots that summarize a lot (see page 18 of adegenet spca tutorial for descriptions):
pdf("EDA/SCRIPTS/Rscripts/Routput/sPCA-6plot-summary.pdf")
plot(mySpca)
dev.off()
# "Here, [plots 2, 3, and 4] show that genotypes are split into two genetic clusters, one in the west[ish] (or left) and one in the east[ish] (right)"
# Sparse panthers/data in the SW may be contributing to NE skew
```

## Global and Local Tests of Landscape Patterns of Genetic Variation
Global patterns (in chromosomal genetic variation) tend to be more indicative of neutral processes, whereas local patterns (of smaller regions genetic variation) may be more indicative of selective (or other non-neutral) processes.
* The global and local tests proposed can be used to reinforce the decision of whether or not to further interpret (or not interpret) global and local structures of genetic variation
* If global tests are significant, suggests that some neutral process is likely affecting landscape patterns of genetic variation
* If local tests are significant, suggests that adaptation is likely affecting landscape patterns of genetic variation (I don't see this with the pumas however, which is reasonable as theirs only 25 SNPs)

### `Mantel-sPCA.R` (part 3)
```R
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
```


# Estimating Population Structure

## STRUCTURE
First, need to convert genotype data to genind object, and then to structure format (can do this locally).
* Note that "coombs" refers to the University of Central Florida's Genomics and Bioinformatics Cluster high performance computing cluster (HPC)

### `genind2structure.R`
```R
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
genind2structure(datagen, file="pop-structure/STRUCTURE/SCRIPTS/Routput/Puma_SNP_genotypes_cut_AZ_final.str") # Note no `pops=TRUE`
# IMPORTANT: 
# Afterwards, make sure your `.str` file's header only has allele names. If you have a column name/header for individuals (e.g "ind"), populations (e.g. "pop"), etc. You may have to manually delete this after running the `genindstructure.R ` script.
```

### Important note on STRUCTURE input file (`.str`) headers
IMPORTANT: Make sure your `.str` file's header only has allele names. If you have a column name/header for individuals (e.g "ind"), populations (e.g. "pop"), etc. You may have to manually delete this after running the `genindstructure.R ` script. Example:

* Correct:
  ```bash
  PP01  PP02  PP03  PP04
  ind1  1 1 1 -9
  ind2  1 1 1 -9
  ind3  1 1 1 1
  ```
* Incorrect:
  ```bash
  ind PP01  PP02  PP03  PP04
  ind1  1 1 1 -9
  ind2  1 1 1 -9
  ind3  1 1 1 1
  ```

### STRUCTURE Param Files

Param files are required to specify parameters used to run STRUCTURE
* `mainparams` file with parameters used for STRUCTURE run
* `extraparams` file with optional parameters used for STRUCTURE run

### `mainparams`
See STRUCTURE manual here for info on each parameter (Section 7.2): https://www.ccg.unam.mx/~vinuesa/tlem09/docs/structure_doc.pdf 
* Additional 'mainparams'file info: http://wiki.univ-reunion.fr/ccur/index.php/Structure
* IMPORTANT!: Be careful of confusing location (e.g. LOCDATA) with loci (e.g. NUMLOCI)

```bash
#define BURNIN 100000
#define NUMREPS 500000

#define INFILE Puma_SNP_genotypes_cut_AZ_final.str

#define NUMINDS 536
#define NUMLOCI 25
#define PLOIDY 2
#define MISSING -9
#define ONEROWPERIND 0
#define LABEL 1
#define POPDATA 0
#define POPFLAG 0
#define LOCDATA 0
#define PHENOTYPE 0
#define EXTRACOLS 0
#define MARKERNAMES 1
#define RECESSIVEALLELES 0
#define MAPDISTANCES 0

#define PHASED 0
#define PHASEINFO 0
#define MARKOVPHASE 0
```

### `extraparams`
```bash
#define NOADMIX 0
#define LINKAGE 0
#define USEPOPINFO 0
#define LOCPRIOR 0
#define FREQSCORR 1
#define ONEFST 0
#define INFERALPHA 1
#define POPALPHAS 0
#define ALPHA 1.0
#define INFERLAMBDA 0
#define POPSPECIFICLAMBDA 0
#define LAMBDA 1

#define UNIFPRIORALPHA 1
#define ALPHAMAX 10.0
#define ALPHAPROPSD 0.025

#define FPRIORMEAN 0.01
#define FPRIORSD 0.05
```

### Running STRUCTURE (on HPC)
Finally, to run structure, you need three files:
* `.str`file with STRUCTURE formatted genotype matrix
* `mainparams`
* `extraparams`

Helpful Tips:
* Suggestions on quality control: https://www.youtube.com/watch?v=ORoSZrqiJ1E
* Guide to running structure in HPC: https://cahuparo.com/2019/10/12/how-to-run-structure-in-henry2/

NOTE: When running structure, need to state the number of pops to run using the `#define MAXPOPS N` line in the `mainparams` file. You do not need this if using the `-K` option in a loop (e.g., `structure -K "$i" -0 outputfile"$j"`)

```bash
#!/bin/bash
#author: ge199066
#SBATCH -J structure
#SBATCH -o structure.out
#SBATCH -e structure.err
#SBATCH -p normal
#SBATCH --cpus-per-task=16
#SBATCH -t 7-00:00:00
#SBATCH --mem-per-cpu=16000

cd /home/ge199066/projects/puma/sumstats/STRUCTURE

module load structure

for i in {1..10};
do
  for j in {1..10};
  do
    structure -K "$i" -0 outputfile"$j"
  done
done
```

### Interpreting STRUCTURE Output With `pophelper` R Package
Pophelper is useful for processing and visualizing output of STRUCTURE, as well as running evanno tests

### `pophelper.R`
```R
##############################################################################################################
# This script is for processing and visualizing output of STRUCTURE, as well as running evanno tests
##############################################################################################################
rm(list=ls())
library(pophelper)
library(ggplot2)
library(gridExtra)

##################################################
# list STRUCTURE output files in character vector
##################################################
sfiles <- list.files("pop-structure/STRUCTURE/SCRIPTS/Routput/coombs/output", full.names=TRUE) # need 'full.names=TRUE' to get correct path to folder with files
sfiles <- sfiles[sfiles != "pop-structure/STRUCTURE/SCRIPTS/Routput/coombs/output/seed.txt"] # remove seed file (otherwise can't convert to qlist)

# convert run files (q-matrices) to qlist
slist <- readQ(sfiles, filetype="structure", indlabfromfile = TRUE)
  
# some basic summary stats
tabulateQ(qlist=readQ(sfiles))
head(summariseQ(tabulateQ(slist)))

##################################################
# EVANNO METHOD
##################################################
# Used to estimate the optimal number of K. The summarised runs table output fromsummariseQ() function can be input to evannoMethodStructure(). 
sr1 <- summariseQ(tabulateQ(slist))
evannoMethodStructure(data=sr1)

# Evanno graphs
# See 'https://www.rdocumentation.org/packages/pophelper/versions/2.3.1/topics/evannoMethodStructure' for options
p <- evannoMethodStructure(data=sr1,
                           exportplot=F,
                           returnplot=T,
                           returndata=F,
                           basesize=12,
                           linesize=0.7,
                           xaxisbreaks=c(1,2,3,4,5,6,7,8,9,10),
                           xaxislabels=c(1,2,3,4,5,6,7,8,9,10))

pdf("pop-structure/STRUCTURE/SCRIPTS/Routput/Rplot_STRUCTURE-evanno-graphs.pdf")
# Actual plots (i.e. arrange the grid of plots)
grid.arrange(p)
# Close the PDF device
dev.off()

# Barplots
# 'alignK' aligns/orders 'slist' names for easy grabbing of files for plotting
slist1 <- alignK(slist[c(31:33, 51:53)]) # from list, choose runs/reps of same K

# see 'https://www.rdocumentation.org/packages/pophelper/versions/2.3.1/topics/plotQ' for options
# use `splab` to label strip panels, e.g. `splab=c("test1","test2","test3","test4","test5","test6")`
p1 <- plotQ(slist1, imgoutput="join", returnplot=T, exportplot=F, basesize=11)
pdf("pop-structure/STRUCTURE/SCRIPTS/Routput/Rplot_STRUCTURE-barplots.pdf")
grid.arrange(p1$plot[[1]])
dev.off()
```

# Fixation Index (Fst) Analysis
Pumas were assigned to subpopulations based on patch type. Then, Fst was calculate by permutating individuals across patches and calculating pair-wise Fst between patches. This was repeated over several iterations to generate a null distribution of Fst values. Afterwards, the Fst calculated using the observed data ("vanilla Fst", using WC84 method) was compared to the null distribution, then a p-value and Z-score were generated.

## Fst Permutations
The following script was used to perform Fst permutations.

### `Fst-Permutation.R`
```R
##################################################################################
# GET DATA, GENERATE MERGED GENOTYPE DATA AND COMBINE PATCH+GENOTYPE DATA 
# Used for all patch types
##################################################################################
rm(list=ls())
library(adegenet)
library(hierfstat)
library(pegas)
library(ggplot2)
patch_data <- read.csv("GIS/DATA/Created/Puma_patches.csv", header=T)
genotype_data <- read.table("SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final.tab", header = TRUE, sep = "\t")
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

# The following blocks are used to generate permutation data for ecoregions, roads, rivers, etc. as needed

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

# BEGIN PERMUTATION LOOP
#results <- NULL # create empty df to store results
results <- data.frame(Fst=numeric(), pair=character())

for (j in 1:nrow(pair)) {
  for (i in 1:1000) { # set number of reps (10, 100, 1000, etc.)
    d1 <- subset(data, US_L3NAME_num %in% pair[j,])
    d1$US_L3NAME_num <- sample(d1$US_L3NAME_num) # shuffling population assigned to each individual
    d2 <- d1[,!names(data) %in% c("Puma_ID", "Latitude", "Longitude","Datset","rivers","roads","US_L3CODE","US_L3NAME")] # removes unnecessary columns
    d2 <- d2[, c(ncol(d2), 1:(ncol(d2)-1))] # moving column from end of df to beginning
    dg1 <- df2genind(d2[,-1], ploidy=2, ncode=1, pop=d2$US_L3NAME_num) # IMPORTANT: Need to make sure first column has population info and that the 'df2genind' object ignores the first column (i.e., d2[,-1]) and instead specifies population separately (i.e., pop=d2$eco3)
    matFst <- genet.dist(dg1, method = "WC84") 
    output <- data.frame(Fst=matFst[1], pair=pair[j,])
    results <- rbind(results, output) # store the result in the df
  }
}

results$paircombo <- paste(results$pair.a, results$pair.b, sep="") # create new column that shows which two regions were compaired (useful for plotting)

# CSV OF PERMUTATION DATA
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Permutation-US_L3NAME_num_04092024.csv")

# PLOTS
# quartz()
pdf("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-US_L3NAME_num_04092024-violin-plot.pdf")
ggplot(data=results, aes(paircombo, Fst, group=paircombo)) + 
  geom_violin()
dev.off()

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

# BEGIN PERMUTATION LOOP
#results <- NULL # create empty df to store results
results <- data.frame(Fst=numeric(), pair=character())

for (j in 1:nrow(pair)) {
  for (i in 1:1000) { # set number of reps (10, 100, 1000, etc.)
    d1 <- subset(data, rivers %in% pair[j,])
    d1$rivers <- sample(d1$rivers) # shuffling population assigned to each individual
    d2 <- d1[,!names(data) %in% c("Puma_ID", "Latitude", "Longitude","Datset","roads","US_L3CODE","US_L3NAME", "US_L3NAME_num")] # removes unnecessary columns
    d2 <- d2[, c(ncol(d2), 1:(ncol(d2)-1))] # moving column from end of df to beginning
    dg1 <- df2genind(d2[,-1], ploidy=2, ncode=1, pop=d2$rivers) # IMPORTANT: Need to make sure first column has population info and that the 'df2genind' object ignores the first column (i.e., d2[,-1]) and instead specifies population separately (i.e., pop=d2$eco3)
    matFst <- genet.dist(dg1, method = "WC84") 
    output <- data.frame(Fst=matFst[1], pair=pair[j,])
    results <- rbind(results, output) # store the result in the df
  }
}

results$paircombo <- paste(results$pair.a, results$pair.b, sep="") # create new column that shows which two regions were compaired (useful for plotting)

# CSV OF PERMUTATION DATA
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Permutation-rivers_04102024.csv")

# PLOTS
# quartz()
pdf("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-rivers_04102024-violin-plot.pdf")
ggplot(data=results, aes(paircombo, Fst, group=paircombo)) + 
  geom_violin()
dev.off()

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

# BEGIN PERMUTATION LOOP
#results <- NULL # create empty df to store results
results <- data.frame(Fst=numeric(), pair=character())

for (j in 1:nrow(pair)) {
  for (i in 1:1000) { # set number of reps (10, 100, 1000, etc.)
    d1 <- subset(data, roads %in% pair[j,])
    d1$roads <- sample(d1$roads) # shuffling population assigned to each individual
    d2 <- d1[,!names(data) %in% c("Puma_ID", "Latitude", "Longitude","Datset","rivers","US_L3CODE","US_L3NAME", "US_L3NAME_num")] # removes unnecessary columns
    d2 <- d2[, c(ncol(d2), 1:(ncol(d2)-1))] # moving column from end of df to beginning
    dg1 <- df2genind(d2[,-1], ploidy=2, ncode=1, pop=d2$roads) # IMPORTANT: Need to make sure first column has population info and that the 'df2genind' object ignores the first column (i.e., d2[,-1]) and instead specifies population separately (i.e., pop=d2$eco3)
    matFst <- genet.dist(dg1, method = "WC84") 
    output <- data.frame(Fst=matFst[1], pair=pair[j,])
    results <- rbind(results, output) # store the result in the df
  }
}

results$paircombo <- paste(results$pair.a, results$pair.b, sep="") # create new column that shows which two regions were compaired (useful for plotting)

# CSV OF PERMUTATION DATA
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Permutation-roads_04102024.csv")

# PLOTS
# quartz()
pdf("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-roads_04102024-violin-plot.pdf")
ggplot(data=results, aes(paircombo, Fst, group=paircombo)) + 
  geom_violin()
dev.off()
```

## Traditional Fst
The following script was used to calculate Fst using traditional methods.

### `Fst-Traditional.R`
```R
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
genotype_data <- read.table("SNPs/DATA/Puma_SNP_genotypes_cut_AZ_final.tab", header = TRUE, sep = "\t")

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
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-US_L3NAME_num_04102024.csv")

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
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-rivers_04102024.csv")

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
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-roads_04102024.csv")
```

## P-Value and Z-score Estimation
Because a single value (Fst vanilla) is being compared to a null distribution generated through permutation, a Z-score can be used as a measure of effect size (as opposed to Cohen's D, which requires comparison of two distributions). P-values were estimated by dividing the proportion of null distribution values greater than or equal to alpha (0.95) by the total number of values in the null distribution.

### `Pvalue-Zscore.R`
```R
rm(list=ls())
library(tidyr)
library(dplyr)

##########################################
# LEVEL III ECOREGIONS
##########################################
# RETRIEVE DATA
ecoreg_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-US_L3NAME_num_04092024.csv")
# rivers_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-rivers_04102024.csv")
# roads_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-roads_04102024.csv")

ecoreg_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-US_L3NAME_num_04102024.csv")
# rivers_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-rivers_04102024.csv")
# roads_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-roads_04102024.csv")


# GENERATE P-VALUES AND Z-SCORES
# FOR ECOREGIONS:
# 1: "Arizona/New Mexico Mountains" 
# 2: "Arizona/New Mexico Plateau"  
# 3: "Chihuahuan Deserts"           
# 4: "Colorado Plateaus"           
# 5: "Madrean Archipelago"         
# 6: "Mojave Basin and Range"      
# 7: "Sonoran Basin and Range"     
# Generate results for p_values, z_scores
# Assuming ecoreg_1Fst is already distinct by paircombo or each paircombo is unique
# Otherwise use the following line: ecoreg_1Fst <- ecoreg_1Fst %>% distinct(paircombo, .keep_all = TRUE)

# Initialize an empty data frame for results
results <- tibble(paircombo = integer(),
                  obs_stat = double(),
                  null_mean = double(),
                  null_sd = double(),
                  p_value = double(),
                  lower_conf = double(),
                  upper_conf = double(),
                  z_score = double())

# Iterate over each paircombo
for(pc in unique(ecoreg_1Fst$paircombo)) {
  
  # Extract observed Fst
  obs_stat <- filter(ecoreg_1Fst, paircombo == pc) %>% pull(Fst) %>% .[1]
  
  # Calculate null distribution statistics
  null_distribution <- filter(ecoreg_null, paircombo == pc) %>% pull(Fst)
  null_mean <- mean(null_distribution)
  null_sd <- sd(null_distribution)
  p_value <- sum(null_distribution >= obs_stat) / length(null_distribution)
  conf_interval <- quantile(null_distribution, probs = c(0.025, 0.975))
  z_score <- (obs_stat - null_mean) / null_sd
  
  # Append results
  results <- bind_rows(results, tibble(paircombo = pc,
                                       obs_stat = obs_stat,
                                       null_mean = null_mean,
                                       null_sd = null_sd,
                                       p_value = p_value,
                                       lower_conf = conf_interval[1],
                                       upper_conf = conf_interval[2],
                                       z_score = z_score))
}

# Z = (x-u)/s
# X is observed test statistic (e.g. vanilla Fst)
# u is mean of null (permutated) distribution
# s is standard deviation of the null distribution
# assumption about population (null) normality is OK because of central limit theorem (more subsampling will evntually result in a normal distribution)

# View results
print(results, n = nrow(ecoreg_1Fst), width = Inf)



##########################################
# RIVERS
##########################################
# RETRIEVE DATA
rm(list=ls())
rivers_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-rivers_04102024.csv")
rivers_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-rivers_04102024.csv")

# Initialize an empty data frame for results
results <- tibble(paircombo = integer(),
                  obs_stat = double(),
                  null_mean = double(),
                  null_sd = double(),
                  p_value = double(),
                  lower_conf = double(),
                  upper_conf = double(),
                  z_score = double())

# Iterate over each paircombo
for(pc in unique(rivers_1Fst$paircombo)) {
  
  # Extract observed Fst
  obs_stat <- filter(rivers_1Fst, paircombo == pc) %>% pull(Fst) %>% .[1]
  
  # Calculate null distribution statistics
  null_distribution <- filter(rivers_null, paircombo == pc) %>% pull(Fst)
  null_mean <- mean(null_distribution)
  null_sd <- sd(null_distribution)
  p_value <- sum(null_distribution >= obs_stat) / length(null_distribution)
  conf_interval <- quantile(null_distribution, probs = c(0.025, 0.975))
  z_score <- (obs_stat - null_mean) / null_sd
  
  # Append results
  results <- bind_rows(results, tibble(paircombo = pc,
                                       obs_stat = obs_stat,
                                       null_mean = null_mean,
                                       null_sd = null_sd,
                                       p_value = p_value,
                                       lower_conf = conf_interval[1],
                                       upper_conf = conf_interval[2],
                                       z_score = z_score))
}

# Z = (x-u)/s
# X is observed test statistic (e.g. vanilla Fst)
# u is mean of null (permutated) distribution
# s is standard deviation of the null distribution
# assumption about population (null) normality is OK because of central limit theorem (more subsampling will evntually result in a normal distribution)

# View results
print(results, n = nrow(rivers_1Fst), width = Inf)


##########################################
# ROADS
##########################################
# RETRIEVE DATA
rm(list=ls())
roads_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-roads_04102024.csv")
roads_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-roads_04102024.csv")

# Initialize an empty data frame for results
results <- tibble(paircombo = integer(),
                  obs_stat = double(),
                  null_mean = double(),
                  null_sd = double(),
                  p_value = double(),
                  lower_conf = double(),
                  upper_conf = double(),
                  z_score = double())

# Iterate over each paircombo
for(pc in unique(roads_1Fst$paircombo)) {
  
  # Extract observed Fst
  obs_stat <- filter(roads_1Fst, paircombo == pc) %>% pull(Fst) %>% .[1]
  
  # Calculate null distribution statistics
  null_distribution <- filter(roads_null, paircombo == pc) %>% pull(Fst)
  null_mean <- mean(null_distribution)
  null_sd <- sd(null_distribution)
  p_value <- sum(null_distribution >= obs_stat) / length(null_distribution)
  conf_interval <- quantile(null_distribution, probs = c(0.025, 0.975))
  z_score <- (obs_stat - null_mean) / null_sd
  
  # Append results
  results <- bind_rows(results, tibble(paircombo = pc,
                                       obs_stat = obs_stat,
                                       null_mean = null_mean,
                                       null_sd = null_sd,
                                       p_value = p_value,
                                       lower_conf = conf_interval[1],
                                       upper_conf = conf_interval[2],
                                       z_score = z_score))
}

# Z = (x-u)/s
# X is observed test statistic (e.g. vanilla Fst)
# u is mean of null (permutated) distribution
# s is standard deviation of the null distribution
# assumption about population (null) normality is OK because of central limit theorem (more subsampling will evntually result in a normal distribution)

# View results
print(results, n = nrow(roads_1Fst), width = Inf)
```


# Modeling
To develop an approach for testing hypothesis regarding the habitat patches (defined by surrounding rivers, roads, and ecoregions) that best explained observed ancestry coefficients (Q), the following types of models were generated using stan in R:

* roads
* rivers
* ecoreg
* roads + rivers
* roads + ecoreg
* rivers + ecoreg
* roads + interstates + ecoreg
* intercept
* spatial (NNGP)

## Retrieving Ancestry Coefficients (Q) from STRUCTURE Output
Convert STRUCTURE output files to qlist to use in models downstream.

### `pophelper_qlist.R`
```R
##############################################################################################################
# STRUCTURE
##############################################################################################################
# This script is for processing and visualizing output of STRUCTURE, as well as running evanno tests
rm(list=ls())
library(pophelper)
library(ggplot2)
library(gridExtra)

##################################################
# list STRUCTURE output files in character vector (from coombs HPC)
##################################################
sfiles <- list.files("pop-structure/STRUCTURE/coombs/DATA/output", full.names=TRUE) # need 'full.names=TRUE' to get correct path to folder with files
sfiles <- sfiles[sfiles != "pop-structure/STRUCTURE/coombs/DATA/output/seed.txt"] # remove seed file (otherwise can't convert to qlist)

# convert run files (q-matrices) to qlist
slist <- readQ(sfiles, filetype="structure", indlabfromfile = TRUE)
slistK <- alignK(slist)
slist1 <- slist[c(31:33, 51:53)] # test line
slist2 <- alignK(slist[c(31:33, 51:53)]) # test line

# Save slistK for downstream analyses
saveRDS(slistK, file = "pop-structure/STRUCTURE/SCRIPTS/Routput/slistK.rds")
# slistK <- readRDS("slistK.rds") # read file later
  
# some basic summary stats
tabulateQ(qlist=readQ(sfiles))
head(summariseQ(tabulateQ(slistK)))
```

This qlist can then be convered to a dataframe to use in R
### `qlist2df.R`
```R
rm(list=ls())
library(pophelper)
library(ggplot2)
library(gridExtra)

###################################################################################################
# STRUCTURE
###################################################################################################
# Load qlist
slistK <- readRDS("pop-structure/STRUCTURE/SCRIPTS/Routput/slistK.rds") 

####################
# GETTING K=4 DATA
# Extract data and calculate means
cluster1 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K4_%d_f", i)]][["Cluster1"]]
})
cluster1_means <- rowMeans(cluster1, na.rm = TRUE)

cluster2 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K4_%d_f", i)]][["Cluster2"]]
})
cluster2_means <- rowMeans(cluster2, na.rm = TRUE)

cluster3 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K4_%d_f", i)]][["Cluster3"]]
})
cluster3_means <- rowMeans(cluster3, na.rm = TRUE)

cluster4 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K4_%d_f", i)]][["Cluster4"]]
})
cluster4_means <- rowMeans(cluster4, na.rm = TRUE)

K4_data <- cbind(cluster1_means, cluster2_means, cluster3_means, cluster4_means)
STRUCTURE_data_K4 <- as.data.frame(K4_data)
colnames(STRUCTURE_data_K4) <- c("Pop1", "Pop2", "Pop3", "Pop4")

# Save to file for downstream analyses
saveRDS(STRUCTURE_data_K4, file = "modeling/input_data/STRUCTURE_data_K4.rds")
```


## R Script for Modeling
The following is the R script use to run the most supported model (roads and rivers)
### `K4_roads_riv.R`
```R
########################################################
# RUNNING THROUGH `rstan`
########################################################
rm(list=ls())
library(rstan)
library(dplyr)
library(parallel)
# library(FNN)
library(geosphere)
library(sp)

########## PATH INFO ##########
# Base directory
base_dir <- "" # normally just repo
# Patch data directory
patch_dir <- paste0(base_dir, "GIS/DATA/Created/")
# Ancestry data directory
Q_dir <- paste0(base_dir, "modeling/input_data/")
# Stan directory
stan_dir <- paste0(base_dir, "modeling/bayesian/local/stan/")
# Output directory
output_dir <- paste0(base_dir, "modeling/bayesian/local/Routput/")
# R directory
R_dir <- paste0(base_dir, "modeling/bayesian/local/")
# model type (predictor(s))
predA <- "roads_" # useful for file naming
predB <- "rivers_" # useful for file naming

# Ancestry data file
Q_file <- "STRUCTURE_data_K4.rds"
# K cluster prefix
K_prefix <- "K4_"

# model specs
chains <- 4
iter <- 15000
warmup <- 7500
####################

########################################################
# DATA PREP
########################################################
# ======================== ANCESTRY COEFFICIENT DATA ============================
ancestry_coeff_data <- readRDS(paste0(Q_dir, Q_file))

# Normalize each row to sum to 1
Q <- sweep(ancestry_coeff_data, 1, rowSums(ancestry_coeff_data), "/") 

# ======================== PATCH DATA ================================
# Load and clean up data so that 'interstates' is defined in the dataset and 'Q' is a matrix or dataframe of observed probabilities
patch_data <- read.csv(paste0(patch_dir, "Puma_patches.csv")) # Read the data

########################################################
# RUNNING MODEL THROUGH `rstan`
########################################################
# ======================== ASSIGNING GEOGRAPHIC OUTLIERS ================================
# ------------------------------------------------------------------------------------------------------------
# manually assign pumas outisde of state boundary to river region 2
patch_data <- patch_data %>%
  mutate(rivers = ifelse(Puma_ID == "UA00031261", 2, rivers))
patch_data <- patch_data %>%
  mutate(rivers = ifelse(Puma_ID == "UA00047644", 2, rivers))
# ------------------------------------------------------------------------------------------------------------
# manually assign pumas outisde of state boundary to road region 1
patch_data <- patch_data %>%
  mutate(roads = ifelse(Puma_ID == "UA00031261", 1, roads))
patch_data <- patch_data %>%
  mutate(roads = ifelse(Puma_ID == "UA00047644", 1, roads))
# ------------------------------------------------------------------------------------------------------------

# ======================== COORDINATES FOR NNGP ================================
# Convert coordinates to spatial object
coords <- data.frame(Longitude = patch_data$Longitude, Latitude = patch_data$Latitude)
coordinates(coords) <- ~ Longitude + Latitude
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")

# Transform coordinates to UTM Zone 12N (Arizona)
utm_coords <- spTransform(coords, CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs"))
# Extract UTM coordinates as a matrix
utm_matrix <- as.matrix(coordinates(utm_coords))

# Add small jitter (if more that one puma was harvested at a location; helps NNGP)
dup_idx <- which(duplicated(utm_matrix) | duplicated(utm_matrix, fromLast = TRUE)) # find all duplicated rows
if(length(dup_idx) > 0){ # for each dup location, adds small jitter
  utm_matrix[dup_idx, ] <- utm_matrix[dup_idx, ] + matrix(rnorm(length(dup_idx)*2, 0, 0.1), ncol=2) # rnorm is *2 because need both x and y coords, then just mean 0 sd 0.1; also make sure ncol=2
}

# ======================== DUMMY VARIABLES ================================
# ROADS
# Convert predictor(s) to factor with all expected levels, drop dummy variable, choose base (e.g., base = 1 for first group)
## IMPORTANT!: levels=as.character(1:5) should match what is expected in from preds
patch_data$roads <- factor(patch_data$roads, levels=as.character(1:6)) # ensures correct dummy coding and all levels present
# Remove intercept (first group is base)
dummy_mat_predA <- model.matrix(~ roads, data = patch_data)[, -1, drop=FALSE]

# RIVERS
# Convert predictor(s) to factor with all expected levels, drop dummy variable, choose base (e.g., base = 1 for first group)
## IMPORTANT!: levels=as.character(1:5) should match what is expected in from preds
patch_data$rivers <- factor(patch_data$rivers, levels=as.character(1:5)) # ensures correct dummy coding and all levels present
# Remove intercept (first group is base)
dummy_mat_predB <- model.matrix(~ rivers, data = patch_data)[, -1, drop=FALSE]

# Note: it's good to double check dummy_mat before proceeding
head(dummy_mat_predA)
head(dummy_mat_predB)

# ======================== NNGP/VEcchia Neighbor Matrices ================================
N <- nrow(utm_matrix)
# M <- 10 # Assign number of NN
M <- 10 # Assign number of NN

# Code for Vecchia/NNGP neighbor matrices (following Zheng method)
NN_ind <- matrix(NA, nrow = N - 1, ncol = M)
NN_dist <- matrix(NA, nrow = N - 1, ncol = M)
NN_distM <- matrix(NA, nrow = N - 1, ncol = M * (M - 1) / 2)

for (i in 2:N) {
  candidate_idx <- 1:(i-1)
  this_M <- min(M, length(candidate_idx))
  dists <- sqrt(rowSums((utm_matrix[candidate_idx, , drop=FALSE] - utm_matrix[i, ])^2))
  nn <- order(dists)[1:this_M]
  NN_ind[i-1, 1:this_M] <- candidate_idx[nn]
  NN_dist[i-1, 1:this_M] <- dists[nn]
  # make pairwise distance matrix for the neighbors
  neighbor_coords <- utm_matrix[NN_ind[i-1, 1:this_M], , drop=FALSE]
  pairwise <- as.matrix(dist(neighbor_coords))
  tri_vec <- rep(NA, M * (M - 1) / 2)
  if(this_M > 1){
    tri_vec[1:(this_M * (this_M - 1) / 2)] <- pairwise[lower.tri(pairwise)]
  }
  NN_distM[i-1, ] <- tri_vec
}

# This block helps avoid errors when running stan model
# Replace any NA in NN_ind with 1
NN_ind[is.na(NN_ind)] <- 1
# Replace any NA in NN_dist and NN_distM with 0
NN_dist[is.na(NN_dist)] <- 0
NN_distM[is.na(NN_distM)] <- 0

# Note: It's good to check matrix in case anything looks weird

# ======================== RUN MODEL ================================
stan_data <- list(
  N = nrow(Q), # number of observations
  K = ncol(Q), # number of ancestral groups
  PA = ncol(dummy_mat_predA), # number of predictors (e.g., road dummies)
  PB = ncol(dummy_mat_predB), # number of predictors (e.g., river dummies)
  roads = dummy_mat_predA,
  rivers = dummy_mat_predB,
  Q = Q, # response variables for Dirichlet
  M = M, # number of NN
  NN_ind = NN_ind,
  NN_dist = NN_dist,
  NN_distM = NN_distM
)

# Run the model
## Sometimes you may need to restart R or it will have compilation issues
# fit <- stan(file = paste0(stan_dir, "roads_riv.stan"), data = stan_data, cores=detectCores(), chains = chains, iter = iter, warmup = warmup)
# Nonspatial
fit <- stan(file = paste0(stan_dir, "/nonspatial/roads_riv_NS.stan"), data = stan_data, cores=detectCores(), chains = chains, iter = iter, warmup = warmup)

# save fit
saveRDS(fit, file = paste0(output_dir, K_prefix, predA, predB, "fit.rds"))


###################################
##### DirichletReg Comparison #####
###################################
# DirichReg
library(DirichletReg)
ancestry_coeff_data <- readRDS(paste0(Q_dir, Q_file))   
Q <- sweep(ancestry_coeff_data, 1, rowSums(ancestry_coeff_data), "/") 
Q_DR <- DR_data(Q)
# mm <- model.matrix(~ patch_data$roads, data = patch_data)
# fitDir <- DirichReg(Q_DR ~ patch_data$roads, model="common")
fitDir <- DirichReg(Q_DR ~ patch_data$roads + patch_data$rivers, model="alternative")
summary(fitDir)
head(fitted(fitDir))

# DirichReg-like output of fitted meanss 
# Extract posterior draws for mu_out (dimensions: [iterations, N, K])
mu_draws <- rstan::extract(fit, "mu_pred")$mu_pred
# Compute posterior mean for each observation/component (mean over iterations)
stan_mu <- apply(mu_draws, c(2, 3), mean)  # [N, K]
# Assign column names to match DirichletReg output
colnames(stan_mu) <- c("Pop1", "Pop2", "Pop3", "Pop4")

# Show first few rows, just like DirichletReg
head(stan_mu)
```


## Stan Scripts Specifying Model and Parameters
The following stan script was used for the best supported model (roads and rivers).

### `roads_riv.stan`
```stan
functions {
  // NNGP prior on latent w (consistent with Lu Zhang's case study)
  real nngp_w_lpdf(vector w, real sigmasq, real phi,
                   matrix NN_dist, matrix NN_distM, int[,] NN_ind,
                   int N, int M) {
    vector[N] V;
    vector[N] I_Aw = w;
    int dim;
    int h;

    for (i in 2:N) {
      matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M ] iNNdistM;
      matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M ] iNNCholL;
      vector[ i < (M + 1) ? (i - 1) : M ] iNNcorr;
      vector[ i < (M + 1) ? (i - 1) : M ] v;
      row_vector[ i < (M + 1) ? (i - 1) : M ] v2;

      dim = (i < (M + 1)) ? (i - 1) : M;

      // Build neighbor covariance matrix for node i
      if (dim == 1) {
        iNNdistM[1, 1] = 1;
      } else {
        h = 0;
        for (j in 1:(dim - 1)) {
          for (k in (j + 1):dim) {
            h = h + 1;
            iNNdistM[j, k] = exp(-phi * NN_distM[(i - 1), h]);
            iNNdistM[k, j] = iNNdistM[j, k];
          }
        }
        for (j in 1:dim) {
          iNNdistM[j, j] = 1;
        }
      }
      iNNCholL = cholesky_decompose(iNNdistM);
      // Vectorized construction of iNNcorr
      iNNcorr = to_vector(exp(-phi * NN_dist[(i - 1), 1:dim]));

      v = mdivide_left_tri_low(iNNCholL, iNNcorr);
      V[i] = 1 - dot_self(v);
      v2 = mdivide_right_tri_low(v', iNNCholL);

      // Vectorized neighbor update for I_Aw[i]
      I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];
    }
    V[1] = 1;
    return -0.5 * (1 / sigmasq * dot_product(I_Aw, (I_Aw ./ V)) +
                   sum(log(V)) + N * log(sigmasq));
  }
}

data {
  int<lower=1> N; // number of observations/individuals
  int<lower=2> K; // number of ancestral groups
  int<lower=1> PA; // number of predictor A categories/regions
  int<lower=1> PB; // number of predictor B categories/regions)
  matrix[N, PA] roads; // roads category for each observation (dummy categories)
  matrix[N, PB] rivers; // rivers category for each observation (dummy categories)
  simplex[K] Q[N]; // observed probs ancestral groups
  int<lower=1> M; // # of neighbors

  int<lower=1, upper=N> NN_ind[N - 1, M];
  matrix[N - 1, M] NN_dist;
  matrix[N - 1, (M * (M - 1)) / 2] NN_distM;
}

parameters {
  vector[K-1] beta0; // intercept for K-1 groups (1st one becomes baseline)
  matrix[K-1, PA] beta_roads; // coefs for each predictor A and group (i.e. not base)
  matrix[K-1, PB] beta_rivers; // coefs for each predictor B and group (i.e. not base)
  real<lower=0> phi_Dirich; // precision parameter for Dirichlet
  real<lower=0> phi; // for NNGP 
  real<lower=0> sigmasq; // spatial process sd (used in NNGP)
  vector[N] w; // spatial random effect (used in NNGP)
}

model {
  beta0 ~ normal(0, 10);
  to_vector(beta_roads) ~ normal(0, 10);
  to_vector(beta_rivers) ~ normal(0, 10);
  phi ~ gamma(1, 1);  
  phi_Dirich ~ lognormal(log(10), 1);
  sigmasq ~ lognormal(log(10), 1);

  // NNGP joint prior for w
  w ~ nngp_w(sigmasq, phi, NN_dist, NN_distM, NN_ind, N, M);

  // Dirichlet regression likelihood
  for (n in 1:N) {
    vector[K] eta; // linear predictor (for each ancestral group)
    vector[K] mu; // expected proportions
    eta[1] = 0; // first group as baseline
    for (k in 2:K)
      eta[k] = beta0[k-1] + roads[n] * beta_roads[k-1]' + rivers[n] * beta_rivers[k-1]'+ w[n];
    mu = softmax(eta); // convert linear predictor to proportions
    Q[n] ~ dirichlet(mu * phi_Dirich); // likelihood for Dirichlet for observed proportions
  }
}

generated quantities {
  vector[K] mu_pred[N];
  simplex[K] Q_pred[N];
  vector[N] log_lik;
  for (n in 1:N) {
    vector[K] eta;
    eta[1] = 0;
    for (k in 2:K)
      eta[k] = beta0[k-1] + roads[n] * beta_roads[k-1]' + rivers[n] * beta_rivers[k-1]' + w[n];
    mu_pred[n] = softmax(eta);
    Q_pred[n] = dirichlet_rng(mu_pred[n] * phi_Dirich);
    log_lik[n] = dirichlet_lpdf(Q[n] | mu_pred[n] * phi_Dirich);
  }
}

```


## Model Summaries
After running models, the following script can be used to evaluate convergence metrics, Bayesian p-values, pairs plots, and other model summaries.
* Note that predAB refers to two parameters include in the model (Parameter A and Parameter B), though the script is generalizeable

### `2b_output_and_plots_predAB.R`
```R
########################################################
# MODEL OUTPUT
########################################################
# OUTPUT SUMMARY
library(rethinking)
sink(paste0(output_dir, K_prefix, predA, predB, "precis.txt")) # like 'pdf()' but for text files
precis(fit, depth=3, prob = 0.95)
sink()

# pairs(fit, pars = c("phi_Dirich", "phi", "sigmasq"))

########################################################
# MODEL EVALUATION
########################################################
# ======================== LOO VALIDATION CRITERIA ================================
library(loo)
loo_result <- loo(fit, pars = "log_lik")
sink(paste0(output_dir, K_prefix, predA, predB, "loo.txt"))  # like 'pdf()' but for text files
print(loo_result)
sink()

# ======================== PREDICTOR COLLINEARITY ================================
library(car)

cor_mat_predA <- cor(dummy_mat_predA) # make correlation matrix,
# Variance Inflation Factor (VIF)
dummy_df_predA <- as.data.frame(dummy_mat_predA) # convert matrix to df
vif_model_predA <- lm(Q[,1] ~ ., data = dummy_df_predA) # predA with Q[,1] as response and all columns in dummy_df as preds (necessary for VIF)
vif_vals_predA <- vif(vif_model_predA) # get VIFs
sink(paste0(output_dir, K_prefix, predA, predB, "VIF_results_", predA, ".txt"))
print(vif_vals_predA)
sink()

cor_mat_predB <- cor(dummy_mat_predB) # make correlation matrix,
# Variance Inflation Factor (VIF)
dummy_df_predB <- as.data.frame(dummy_mat_predB) # convert matrix to df
vif_model_predB <- lm(Q[,1] ~ ., data = dummy_df_predB) # predB with Q[,1] as response and all columns in dummy_df as preds (necessary for VIF)
vif_vals_predB <- vif(vif_model_predB) # get VIFs
sink(paste0(output_dir, K_prefix, predA, predB, "VIF_results_", predB, ".txt"))
print(vif_vals_predB)
sink()



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### # ONLY FILE NAMING IS DIFFERENT IN THIS BLOCK # ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

# ======================== RESIDUAL EVALUATION ================================
mu_pred <- rstan::extract(fit)$mu_pred # extract predicte mu
mu_pred_mean <- apply(mu_pred, c(2,3), mean) # expected matrix [N,K] (i.e. what c(2,3) is grabbing)

# grab vector of posterior phi draws
phi_mean <- mean(rstan::extract(fit)$phi)  # or median()
# Compute model based Dirichlet variance for each fitted value
V_dirich <- mu_pred_mean * (1 - mu_pred_mean) / (1 + phi_mean) # [N,K]

##### ##### #####
# "Standardized" residuals (trying to match DirichletReg R package method)
pop <- 4 ## Need to do for each pop
std_resids_pop <- (Q[,pop] - mu_pred_mean[,pop]) / sqrt(V_dirich[,pop])
# quantiles summary
sink(paste0(output_dir, K_prefix, predA, predB, "std_resids_quantiles_pop", pop, ".txt"))
quantile(std_resids_pop, probs = c(0, 0.25, 0.5, 0.75, 1)) # between -2 and 2 is generally ideal
sink()
# Q-Q plot
pdf(paste0(output_dir, K_prefix, predA, predB, "qqnorm_pop", pop, ".pdf"))
qqnorm(std_resids_pop, main=paste0("Q-Q plot (pop", pop, ")"))
qqline(std_resids_pop)
dev.off()
# ======================== RESIDUAL SPATIAL AUTOCORRELATION ================================
# Using same pop as specified above
# KNN appropriate 
library(spdep)
knn <- knearneigh(utm_matrix, k=10) # get 10 NN again
nb <- knn2nb(knn) # conver to neighbors list (class 'nb')
lw <- nb2listw(nb, style="W") # convert neighbot list to spatial weights, row standardized ("W")
moran_result <- moran.test(std_resids_pop, lw) 
sink(paste0(output_dir, K_prefix, predA, predB, "MoransI_pop", pop, ".txt"))
print(moran_result)
sink()
##### ##### #####

# Optional: Which individual observations are extreme? (e.g., |residual| > 3 is a often cutoff for outliers)
extreme_idx <- which(abs(std_resids_pop) > 3)
# make summary table of indices and values of individuals
test <- data.frame(
  Index = extreme_idx,
  Std_Residual = std_resids_pop[extreme_idx], # value
  Q_Observed = Q[extreme_idx, pop],         # observed
  Fitted = mu_pred_mean[extreme_idx, pop]   # fitted
)
head(test)



##########################################
# POSTERIOR PREDICTIVE CHECK
##########################################
library(bayesplot)

# Grab observed data from earlier
Q_obs <- Q

# # Add observation index and roads category
# Q_obs <- Q_obs %>% mutate(obs = row_number(), roads = factor(roads)) # not necessary

# Assuming Q_pred_draws is an array where:
# - The first dimension corresponds to predictive draws
# - The second dimension corresponds to observations (the same order as in Q_obs)
Q_pred_draws <- rstan::extract(fit)$Q_pred
thin_factor <- 1  # Keep every Nth draw
Q_pred_draws <- Q_pred_draws[seq(1, dim(Q_pred_draws)[1], by = thin_factor), , ]

# Get observed data for each pop and convert to numeric
Q_obs_pop1 <- as.numeric(Q_obs$Pop1)
Q_obs_pop2 <- as.numeric(Q_obs$Pop2)
Q_obs_pop3 <- as.numeric(Q_obs$Pop3)
Q_obs_pop4 <- as.numeric(Q_obs$Pop4)

# No need to transpose as it seems we need [draws, observations]
Q_pred_pop1_corrected <- Q_pred_draws[, , 1]  # Assuming the 1st ancestral group is what we're focusing on
Q_pred_pop2_corrected <- Q_pred_draws[, , 2]
Q_pred_pop3_corrected <- Q_pred_draws[, , 3]
Q_pred_pop4_corrected <- Q_pred_draws[, , 4]

# Print test statistic (e.g. mean) vs posterior hist
pdf(paste0(output_dir, K_prefix, predA, predB, "ppc_mean_Pop1.pdf"))
ppc_stat(y=Q_obs_pop1, yrep=Q_pred_pop1_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "ppc_mean_Pop2.pdf"))
ppc_stat(y=Q_obs_pop2, yrep=Q_pred_pop2_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "ppc_mean_Pop3.pdf"))
ppc_stat(y=Q_obs_pop3, yrep=Q_pred_pop3_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "ppc_mean_Pop4.pdf"))
ppc_stat(y=Q_obs_pop4, yrep=Q_pred_pop4_corrected, stat = "mean")
dev.off()

#########################################
# BAYESIAN P-VALUES
#########################################
# Extract log-likelihood
library(loo)
log_lik <- extract_log_lik(fit, parameter_name = "log_lik")

# Perform posterior predictive checks with loo
loo_object <- loo(log_lik)

# Simple function to calculate Bayesian p-values
bayesian_p_value <- function(y, yrep) {
  mean(yrep >= y)
}

# Calculate Bayesian p-values for each population
p_value_pop1 <- bayesian_p_value(Q_obs_pop1, Q_pred_pop1_corrected)
p_value_pop2 <- bayesian_p_value(Q_obs_pop2, Q_pred_pop2_corrected)
p_value_pop3 <- bayesian_p_value(Q_obs_pop3, Q_pred_pop3_corrected)
p_value_pop4 <- bayesian_p_value(Q_obs_pop4, Q_pred_pop4_corrected)
# Print 
sink(paste0(output_dir, K_prefix, predA, predB, "bayesian_p_values.txt"))  # like 'pdf()' but for text files
print(paste("Bayesian p-value for Pop1:", p_value_pop1))
print(paste("Bayesian p-value for Pop2:", p_value_pop2))
print(paste("Bayesian p-value for Pop3:", p_value_pop3))
print(paste("Bayesian p-value for Pop4:", p_value_pop4))
sink()

graphics.off()

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### END OF ONLY FILE NAMING ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 




########################################################
# PLOTTING
########################################################
library(ggplot2) 

########## MCMC trace plots ##########
pdf(paste0(output_dir, K_prefix, predA, predB, "trace_plot_beta0.pdf"))
stan_trace(fit, pars = c("beta0"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()



pdf(paste0(output_dir, K_prefix, predA, predB, "trace_plot_beta_roads.pdf"))
stan_trace(fit, pars = c("beta_roads"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "trace_plot_beta_rivers.pdf"))
stan_trace(fit, pars = c("beta_rivers"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "trace_plot_beta_ecoregs.pdf"))
stan_trace(fit, pars = c("beta_ecoregs"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()


pdf(paste0(output_dir, K_prefix, predA, predB, "trace_plot_phi.pdf"))
stan_trace(fit, pars = c("phi"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "trace_plot_phi_Dirich.pdf"))
stan_trace(fit, pars = c("phi_Dirich"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "trace_plot_sigmasq.pdf"))
stan_trace(fit, pars = c("sigmasq"))
dev.off()

graphics.off()

########## Fit plot ##########
pdf(paste0(output_dir, K_prefix, predA, predB, "fit_plot_beta0.pdf"))
plot(fit, pars=c("beta0"))
dev.off()



pdf(paste0(output_dir, K_prefix, predA, predB, "fit_plot_beta_roads.pdf"))
plot(fit, pars=c("beta_roads"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "fit_plot_beta_rivers.pdf"))
plot(fit, pars=c("beta_rivers"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "fit_plot_beta_ecoregs.pdf"))
plot(fit, pars=c("beta_ecoregs"))
dev.off()


pdf(paste0(output_dir, K_prefix, predA, predB, "fit_plot_phi.pdf"))
plot(fit, pars=c("phi"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "fit_plot_phi_Dirich.pdf"))
plot(fit, pars=c("phi_Dirich"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, "fit_plot_sigmasq.pdf"))
plot(fit, pars=c("sigmasq"))
dev.off()

graphics.off()

########## Pairs plots ##########
# for (i in 1:3) { # beta0 is length 3
#   pdf(paste0(output_dir, K_prefix, predA, predB, "pairs_plot_beta0_", i, ".pdf"), width = 10, height = 10)
#   pairs(fit, pars = paste0("beta0[", i, "]"), las = 1)
#   dev.off()
# }

for (i in 1:3) { # for each row of beta
  param_names <- paste0("beta_roads[", i, ",", 1:5, "]")
  pdf(paste0(output_dir, K_prefix, predA, predB, "pairs_plot_beta_roads_", i, ".pdf"), width = 10, height = 10)
  pairs(fit, pars = param_names, las = 1)
  dev.off()
}

for (i in 1:3) { # for each row of beta
  param_names <- paste0("beta_ecoregs[", i, ",", 1:4, "]")
  pdf(paste0(output_dir, K_prefix, predA, predB, "pairs_plot_beta_ecoregs_", i, ".pdf"), width = 10, height = 10)
  pairs(fit, pars = param_names, las = 1)
  dev.off()
}

for (i in 1:3) { # for each row of beta
  param_names <- paste0("beta_rivers[", i, ",", 1:4, "]")
  pdf(paste0(output_dir, K_prefix, predA, predB, "pairs_plot_beta_rivers_", i, ".pdf"), width = 10, height = 10)
  pairs(fit, pars = param_names, las = 1)
  dev.off()
}
```

## STRUCTURE-like Barplots
The following script is used to generate STRUCTURE-like barplots.

### `3a_output_barplots.R`
```R
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(dplyr)

# Assuming Q is your data frame
Q$Individual <- 1:nrow(Q)  # Add an identifier for each individual

############################
# GRAB DATA
############################
# predicted data
Q_pred <- rstan::extract(fit)$Q_pred
# means and CIs
Q_pred_means <- apply(Q_pred, c(2, 3), mean) # Resulting in [individuals, groups]
Q_pred_CI_lower <- apply(Q_pred, c(2, 3), function(x) quantile(x, probs = 0.025))
Q_pred_CI_upper <- apply(Q_pred, c(2, 3), function(x) quantile(x, probs = 0.975))
Q_obs <- Q

############################
# FORMAT DATA FOR PLOTTING
############################
colnames(Q_pred_means) <- c("Pop1", "Pop2", "Pop3", "Pop4")
colnames(Q_pred_CI_lower) <- c("Pop1", "Pop2", "Pop3", "Pop4")
colnames(Q_pred_CI_upper) <- c("Pop1", "Pop2", "Pop3", "Pop4")

Q_pred_means_df <- as.data.frame(Q_pred_means)
Q_pred_means_df$Individual <- 1:nrow(Q_pred_means_df)
Q_pred_means_long <- reshape2::melt(Q_pred_means_df, id.vars = "Individual",
                                    variable.name = "Response_Variable",
                                    value.name = "Mean_Value")

Q_pred_CI_lower_df <- as.data.frame(Q_pred_CI_lower)
Q_pred_CI_lower_df$Individual <- 1:nrow(Q_pred_CI_lower_df)
Q_pred_CI_lower_long <- reshape2::melt(Q_pred_CI_lower_df, id.vars = "Individual",
                                       variable.name = "Response_Variable",
                                       value.name = "CI_Lower")

Q_pred_CI_upper_df <- as.data.frame(Q_pred_CI_upper)
Q_pred_CI_upper_df$Individual <- 1:nrow(Q_pred_CI_upper_df)
Q_pred_CI_upper_long <- reshape2::melt(Q_pred_CI_upper_df, id.vars = "Individual",
                                       variable.name = "Response_Variable",
                                       value.name = "CI_Upper")

# Remove Individual and Type columns
Q_obs <- Q_obs %>%
  select(-Individual)
Q_obs$Individual <- 1:nrow(Q_pred_means_df)

Q_obs_long <- reshape2::melt(Q_obs, id.vars = "Individual",
                             variable.name = "Response_Variable",
                             value.name = "Observed")
###
                    
Q_pred_combined <- merge(Q_pred_means_long, Q_pred_CI_lower_long, 
                         by = c("Individual", "Response_Variable"))
Q_pred_combined <- merge(Q_pred_combined, Q_pred_CI_upper_long, 
                         by = c("Individual", "Response_Variable"))
Q_pred_combined <- merge(Q_pred_combined, Q_obs_long, 
                         by = c("Individual", "Response_Variable"))

#########################
# COMPARATIVE BAR GRAPH
#########################
# pdf(paste0(output_dir, K_prefix, predA, "bar_graph_comparative.pdf")) # 1 predictor
pdf(paste0(output_dir, K_prefix, predA, predB, "bar_graph_comparative.pdf")) # 2 predictor
# pdf(paste0(output_dir, K_prefix, predA, predB, predC, "bar_graph_comparative.pdf")) # 3 predictor

legend_labels <- c("Pop1" = "1", "Pop2" = "2", "Pop3" = "3", "Pop4" = "4")
ggplot(Q_pred_combined, aes(x = factor(Individual), y = Observed, fill = Response_Variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_point(aes(y = CI_Lower, color = "Lower CI"), position = position_dodge(width = 0.9), size = 0.1, shape = 20, show.legend = FALSE) +
  geom_point(aes(y = CI_Upper, color = "Upper CI"), position = position_dodge(width = 0.9), size = 0.1, shape = 20, show.legend = FALSE) +
  geom_point(aes(y = Mean_Value, color = "black"), position = position_dodge(width = 0.9), size = 0.1, shape = 20, show.legend = FALSE) +
  facet_wrap(~ Response_Variable, ncol = 1, scales = "free_y") +  # No need for labeller since we're removing labels
  labs(x = " ", y = "Observed Ancestry Proportions", title = " ", fill = "Population") +  # Add "Population" as legend title
  scale_fill_discrete(labels = legend_labels) +  # Keep default colors, just change labels
  scale_color_manual(values = c("Lower CI" = "blue", "Upper CI" = "red", "Observed" = "black")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        strip.text = element_blank())  # Hide facet labels

dev.off()


#########################
# STRUCTURE LIKE BAR PLOT
#########################
Q_obs_long$Type <- "Observed"
colnames(Q_obs_long) <- c("Individual", "Response_Variable", "Proportion", "Type")
Q_pred_means_long$Type <- "Predicted"
colnames(Q_pred_means_long) <- c("Individual", "Response_Variable", "Proportion", "Type")
obs_vs_pred <- bind_rows(Q_obs_long, Q_pred_means_long)


# pdf(paste0(output_dir, K_prefix, predA, "bar_graph_STRUCTURE_like.pdf"), width = 10, height = 7) # 1 predictor
pdf(paste0(output_dir, K_prefix, predA, predB, "bar_graph_STRUCTURE_like.pdf"), width = 10, height = 7) # 2 predictor
# pdf(paste0(output_dir, K_prefix, predA, predB, predC, "bar_graph_STRUCTURE_like.pdf"), width = 10, height = 7) # 3 predictor

legend_labels <- c("Pop1" = "1", "Pop2" = "2", "Pop3" = "3", "Pop4" = "4")
ggplot(obs_vs_pred, aes(x = Individual, y = Proportion, fill = Response_Variable)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Type, ncol = 1) +
  labs(title = " ",
       x = " ",
       y = " ",
       fill = "Population") +
  scale_fill_discrete(labels = legend_labels) +
  theme_minimal() +
  theme(panel.grid = element_blank(),    # Remove grid lines
        panel.background = element_blank(),  # Remove background
        axis.text.x = element_blank(),    # Hide individual labels
        axis.ticks.x = element_blank(),   # Hide x-axis ticks
        axis.text.y = element_blank(),    # Hide individual labels
        axis.ticks.y = element_blank(),
        strip.text = element_blank())  

dev.off()

```





