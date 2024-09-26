########################################################
# RUNNING THROUGH `rstan`
########################################################
# Need to:
## Prepare Stan Model: Write Stan model code in a text file with a .stan extension or as a string within R. For example, save model code in a file named model.stan
## Prepare Data in R: Organize data in a list that matches the data block in your Stan model
## Run the Model: Use the stan function to compile and fit the model. Need to also specify the path to the Stan model file, the data list, and any additional parameters like the number of chains and iterations.
rm(list=ls())
library(rstan)
library(dplyr)
library(parallel)
library(FNN)
library(geosphere)
########## PATH INFO ##########
# Base directory
base_dir <- "Puma-AZ/Github/puma-AZ/"
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
model <- "ecoreg_" # for barplots

# Ancestry data file
# Q_file <- "STRUCTURE_data_K4.rds"
Q_file <- "TESS_data_K4.rds"
# K cluster prefix
K_prefix <- "K4_"

# model specs
chains <- 4
iter <- 2000
warmup <- 1000
####################

########################################################
# DATA PREP
########################################################
# ======================== ANCESTRY COEFFICIENT DATA ============================
ancestry_coeff_data <- readRDS(paste0(Q_dir, Q_file))   

# Normalize each row to sum to 1
Q <- sweep(ancestry_coeff_data, 1, rowSums(ancestry_coeff_data), "/") 

# Rearrange for TESS data (if Q)file is "TESS_data_K4.rds" (Makes STRUCTURE-like plots nicer)
if (Q_file == "TESS_data_K4.rds") {
  colnames(Q) <- c("Pop3", "Pop1", "Pop4", "Pop2") # better
  Q <- Q[, c("Pop1", "Pop2", "Pop3", "Pop4")] # reorder
}

# ======================== PATCH DATA ================================
# Load and clean up data so that 'interstates' is defined in the dataset and 'Q' is a matrix or dataframe of observed probabilities
# setwd("Puma-AZ/Github/puma-AZ/modeling")
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
# Assign actual values to variables in your .stan script
N <- 536  # number of observations
K <- 3    # number of ancestral groups
E <- 5    # number of "patches"
Q <- Q    # response variables for Dirichlet

ecoreg <- patch_data$US_L3NAME_num # categorical data on ecoreg patch/region

# Convert your coordinates to a spatial object
library(sp)
coords <- data.frame(Longitude = patch_data$Longitude, Latitude = patch_data$Latitude)
coordinates(coords) <- ~ Longitude + Latitude
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")

# Transform the coordinates to UTM Zone 12N
utm_coords <- spTransform(coords, CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs"))

# Extract the UTM coordinates as a matrix
utm_matrix <- as.matrix(coordinates(utm_coords))

# Calculate the Euclidean distance matrix
dist_mat <- as.matrix(dist(utm_matrix))

# DISTANCE MATRIX (GAUSSIAN)
calculate_Z <- function(dist_mat, lambda) {
  N <- nrow(dist_mat)
  Z <- matrix(0, nrow = N, ncol = N)
  
  for (i in 1:N) {
    for (j in 1:N) {
      Z[i, j] <- exp(-(dist_mat[i, j]^2 / lambda^2))
    }
  }
  
  return(Z)
}

lambda <- 40000 # Adjust this value as necessary (often is range/3; range is estimated from semi-variogram)
Z <- calculate_Z(dist_mat, lambda)

# Prepare the Stan data list
stan_data <- list(
  N = nrow(patch_data),
  K = 4,
  E = 5,
  ecoreg = ecoreg,
  Z = Z, # covariance matrix calculated from observed data and dist_mat
  Q = as.matrix(Q)
)

# Run the model
## Sometimes you may need to restart R or it will have compilation issues
## Note that sigma_alpha and mu_alpha may need to be log transformed (to log_sigma_alpha and log_mu_alpha) in model to approach normality
fit <- stan(file = paste0(stan_dir, "ecoreg_RE.stan"), data = stan_data, cores=detectCores(), chains = chains, iter = iter, warmup = warmup)
# print(fit)

saveRDS(fit, file = paste0(output_dir, K_prefix, "ecoreg_RE.rds"))
# fit <- readRDS(file = paste0(output_dir, "riv_RE.rds"))

# OUTPUT SUMMARY
# package `rethinking` is used for the below
library(rethinking) 

sink(paste0(output_dir, K_prefix, "ecoreg_precis.txt")) # like 'pdf()' but for text files
precis(fit, depth=3)
sink()

#precisdf <- as.data.frame(precis(fit, depth=3)) # turning precis into a dataframe (useful if you want to adjust values based on log scale, etc.)
#precisdf
# from: https://mc-stan.org/rstanarm/reference/priors.html#:~:text=If%20concentration%20%3E%201%20%2C%20then%20the,this%20mode%20becomes%20more%20pronounced.
##  If concentration > 1, then the prior mode corresponds to all variables having the same (proportion of total) variance, which can be used to ensure the the posterior variances are not zero. As the concentration parameter approaches infinity, this mode becomes more pronounced. 
## In the unlikely case that concentration < 1, the variances are more polarized.

########################################################
# MODEL EVALUATION
########################################################
# MODEL EVALUATION
library(loo)
loo_result <- loo(fit, pars = "log_lik")

sink(paste0(output_dir, K_prefix, "ecoreg_loo.txt"))  # like 'pdf()' but for text files
print(loo_result)
sink()


##########################################
# POSTERIOR PREDICTIVE CHECK
##########################################
library(bayesplot)

# Grab observed data from earlier
Q_obs <- Q

# Add observation index and ecoreger category
Q_obs <- Q_obs %>% mutate(obs = row_number(), ecoreg = factor(ecoreg))

# Assuming Q_pred_draws is an array where:
# - The first dimension corresponds to predictive draws
# - The second dimension corresponds to observations (the same order as in Q_obs)
Q_pred_draws <- rstan::extract(fit)$Q_pred
thin_factor <- 1  # Keep every 100th draw
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
pdf(paste0(output_dir, K_prefix, "ecoreg_ppc_mean_Pop1.pdf"))
ppc_stat(y=Q_obs_pop1, yrep=Q_pred_pop1_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_ppc_mean_Pop2.pdf"))
ppc_stat(y=Q_obs_pop2, yrep=Q_pred_pop2_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_ppc_mean_Pop3.pdf"))
ppc_stat(y=Q_obs_pop3, yrep=Q_pred_pop3_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_ppc_mean_Pop4.pdf"))
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
sink(paste0(output_dir, K_prefix, "ecoreg_bayesian_p_values.txt"))  # like 'pdf()' but for text files
print(paste("Bayesian p-value for Pop1:", p_value_pop1))
print(paste("Bayesian p-value for Pop2:", p_value_pop2))
print(paste("Bayesian p-value for Pop3:", p_value_pop3))
print(paste("Bayesian p-value for Pop4:", p_value_pop4))
sink()

graphics.off()

########################################################
# PLOTTING
########################################################
library(ggplot2)

########## MCMC trace plots ##########
pdf(paste0(output_dir, K_prefix, "ecoreg_trace_plot_beta_ecoreg.pdf"))
# stan_trace(fit, pars = c("beta_river"))
stan_trace(fit, pars = c("beta_ecoreg"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_trace_plot_gamma.pdf"))
stan_trace(fit, pars = c("gamma"))
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_trace_plot_epsilon.pdf"))
stan_trace(fit, pars = c("epsilon"))
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_trace_plot_mu.pdf"))
stan_trace(fit, pars = c("mu")) # too many
dev.off()

graphics.off()

########## Fit plot ##########
pdf(paste0(output_dir, K_prefix, "ecoreg_fit_plot_beta_ecoreg.pdf"))
plot(fit, pars=c("beta_ecoreg"))
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_fit_plot_gamma.pdf"))
plot(fit, pars=c("gamma"))
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_fit_plot_epsilon.pdf"))
plot(fit, pars=c("epsilon"))
dev.off()

pdf(paste0(output_dir, K_prefix, "ecoreg_fit_plot_mu.pdf"))
plot(fit, pars=c("mu"))
dev.off()

graphics.off()

########## Pairs plots ##########
for (i in 1:K) {
  pdf(paste0(output_dir, K_prefix, "ecoreg_pairs_plot_beta_ecoreg", i, ".pdf"), width = 10, height = 10)
  pairs(fit, pars = paste0("beta_ecoreg[", i, ",", 1:E, "]"), las = 1)
  dev.off()
}

pdf(paste0(output_dir, K_prefix, "ecoreg_pairs_plot_mu.pdf"), width = 10, height = 10)
pairs(fit, pars = "mu", las = 1) # too many
dev.off()

###############################
# BAR PLOTS
###############################
source(paste0(R_dir, "output_barplots.R"))

