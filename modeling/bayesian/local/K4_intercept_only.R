########################################################
# RUNNING THROUGH `rstan`
########################################################
rm(list=ls())
library(rstan)
library(dplyr)
library(parallel)
# library(FNN)
# library(geosphere)
# library(sp)

########## PATH INFO ##########
# Base directory
base_dir <- ""
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
# model
predA <- "intercept_only_" # for naming

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
# rivers <- patch_data$rivers # categorical data on rivers patch/region

# ======================== RUN MODEL ================================
stan_data <- list(
  N = nrow(Q), # number of observations
  K = ncol(Q), # number of ancestral groups
  Q = Q # response variables for Dirichlet
)

# spatial only (testing)
fit <- stan(file = paste0(stan_dir, "/nonspatial/intercept_only.stan"), data = stan_data, cores=detectCores(), chains = chains, iter = iter, warmup = warmup)

# save fit
saveRDS(fit, file = paste0(output_dir, K_prefix, "intercept_only_", "fit.rds"))



