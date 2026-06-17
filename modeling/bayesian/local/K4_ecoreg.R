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
base_dir <- ""
# patch data directory
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
predA <- "ecoreg_" # useful for file naming

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
# Load and clean up data so that predictor is defined in the dataset and 'Q' is a matrix or dataframe of observed probabilities
patch_data <- read.csv(paste0(patch_dir, "Puma_patches.csv")) # Read the data

# ------------------------------------------------------------------------------------------------------------
# SNIPPET FOR ECOREGIONS (ADDS COLUMN WITH ECOREGION REPRESENTED BY NUMBER)
# Assign lone pumas in ecoregion to nearest ecoregion
patch_data$US_L3NAME[patch_data$Puma_ID == "UA00046759"] <- "Madrean Archipelago" # Chihuahuan Desert to Madrean Archipelago
patch_data$US_L3NAME[patch_data$Puma_ID == "UA00013730"] <- "Arizona/New Mexico Plateau" # Colorado Plateaus to Arizona/New Mexico Plateau
# Sort and get unique values from 'US_L3NAME'
US_L3NAME_num <- sort(unique(patch_data$US_L3NAME))
# Convert 'US_L3NAME' to a factor and then to numeric values
patch_data$US_L3NAME_num = as.numeric(factor(patch_data$US_L3NAME, levels = US_L3NAME_num))
unique(patch_data$US_L3NAME_num)
# ------------------------------------------------------------------------------------------------------------

########################################################
# RUNNING MODEL THROUGH `rstan`
########################################################
# ======================== ASSIGNING GEOGRAPHIC OUTLIERS ================================
# N/A
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
# Convert predictor(s) to factor with all expected levels, drop dummy variable, choose base (e.g., base = 1 for first group)
## IMPORTANT!: levels=as.character(1:5) should match what is expected in from preds
patch_data$US_L3NAME_num <- factor(patch_data$US_L3NAME_num, levels=as.character(1:5)) # ensures correct dummy coding and all levels present
# Remove intercept (first group is base)
dummy_mat <- model.matrix(~ US_L3NAME_num, data = patch_data)[, -1, drop=FALSE]

# Note: it's good to double check dummy_mat before proceeding
head(dummy_mat)

# ======================== NNGP/VEcchia Neighbor Matrices ================================
N <- nrow(utm_matrix)
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
  P = ncol(dummy_mat), # number of predictors (e.g., road dummies)
  ecoregs = dummy_mat,
  Q = Q, # response variables for Dirichlet
  M = M, # number of NN
  NN_ind = NN_ind,
  NN_dist = NN_dist,
  NN_distM = NN_distM
)

# Run the model
## Sometimes you may need to restart R or it will have compilation issues
fit <- stan(file = paste0(stan_dir, "ecoreg.stan"), data = stan_data, cores=detectCores(), chains = chains, iter = iter, warmup = warmup)
# Nonspatial:
# fit <- stan(file = paste0(stan_dir, "/nonspatial/ecoreg_NS.stan"), data = stan_data, cores=detectCores(), chains = chains, iter = iter, warmup = warmup)

# save fit
saveRDS(fit, file = paste0(output_dir, K_prefix, predA, "fit.rds"))
