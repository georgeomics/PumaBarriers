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
library(sp)
library(rgdal)
library(gstat)
########## PATH INFO ##########
# Base directory
base_dir <- "puma-AZ/"
# Patch data directory
patch_dir <- paste0(base_dir, "GIS/DATA/Created/")
# Ancestry data directory
Q_dir <- paste0(base_dir, "modeling/input_data/")
# Stan directory
stan_dir <- paste0(base_dir, "modeling/bayesian/local/stan/")
# Output directory
output_dir <- paste0(base_dir, "modeling/bayesian/local/Routput/")

# Ancestry data file
Q_file <- "STRUCTURE_data_K4.rds"
# K cluster prefix
K_prefix <- "K4_"


########################################################
# DATA PREP
########################################################
# ======================== ANCESTRY COEFFICIENT DATA ============================
ancestry_coeff_data <- readRDS(paste0(Q_dir, Q_file))       
# TESS_data_K4 <- readRDS("modeling/input_data/TESS_data_K4.rds") 

# Normalize each row to sum to 1
Q <- sweep(ancestry_coeff_data, 1, rowSums(ancestry_coeff_data), "/") 
# Q <- sweep(TESS_data_K4, 1, rowSums(TESS_data_K4), "/") 

# ======================== PATCH DATA ================================
# Load and clean up data so that 'interstates' is defined in the dataset and 'Q' is a matrix or dataframe of observed probabilities
# setwd("/Users/georgezaragoza/Library/CloudStorage/OneDrive-Personal/UCF/Projects/Puma-AZ/Github/puma-AZ/modeling")
patch_data <- read.csv(paste0(patch_dir, "Puma_patches.csv")) # Read the data


# Convert your coordinates to a spatial object
coords <- data.frame(Longitude = patch_data$Longitude, Latitude = patch_data$Latitude)
coordinates(coords) <- ~ Longitude + Latitude
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")

# Transform the coordinates to UTM Zone 12N
utm_coords <- spTransform(coords, CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs"))

# Extract the UTM coordinates as a matrix
utm_matrix <- as.matrix(coordinates(utm_coords))


##########################
##### SEMI-VARIOGRAM #####
##########################
# Step 1: Prepare the data
# Extract the residuals for Pop4
residuals_Pop4 <- Q$Pop4 - mean(Q$Pop4)

# Step 3: Compute the empirical semi-variogram using the transformed coordinates
variogram_Pop4 <- variogram(residuals_Pop4 ~ 1, locations = utm_coords, width = 10000)

# Step 4: Plot the semi-variogram
plot(variogram_Pop4, main = "Empirical Semi-Variogram for Pop4", xlab = "Distance (meters)", ylab = "Semi-Variance")

# hist(as.vector(dist_mat), breaks = 50, main = "Histogram of Distances", xlab = "Distance (meters)")

# ##### Pop1
# initial_model <- vgm(psill = 0.02, model = "Gau", range = 150000, nugget = 0.03) # adjusting based on looking at plot from variogram_Pop1. Range refers to when the curve starts to generally flatten
# fitted_model_1 <- fit.variogram(variogram_Pop1, model = initial_model)
# print(fitted_model_1)

# ###### Pop2
# initial_model <- vgm(psill = 0.02, model = "Gau", range = 150000, nugget = 0.03) # adjusting based on looking at plot from variogram_Pop2\. Range refers to when the curve starts to generally flatten
# fitted_model_2 <- fit.variogram(variogram_Pop2, model = initial_model)
# print(fitted_model_2)

# ###### Pop3
# initial_model <- vgm(psill = 0.02, model = "Gau", range = 150000, nugget = 0.02) # adjusting based on looking at plot from variogram_Pop2\. Range refers to when the curve starts to generally flatten
# fitted_model_3 <- fit.variogram(variogram_Pop3, model = initial_model)
# print(fitted_model_3)

# ###### Pop4
initial_model <- vgm(psill = 0.02, model = "Gau", range = 150000, nugget = 0.02) # adjusting based on looking at plot from variogram_Pop2\. Range refers to when the curve starts to generally flatten
fitted_model_4 <- fit.variogram(variogram_Pop4, model = initial_model)
print(fitted_model_4)

# Print to see the estimated parameters
print(fitted_model)

# Plot empirical variogram and fitted model together
plot(variogram_Pop4, fitted_model, main = "Empirical Semi-Variogram with Fitted Model")
##### ############## #####
