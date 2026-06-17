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
