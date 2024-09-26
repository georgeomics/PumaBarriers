rm(list=ls())
library(tidyr)
library(dplyr)

##########################################
# LEVEL III ECOREGIONS
##########################################
# RETRIEVE DATA
ecoreg_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-US_L3NAME_num.csv")
# rivers_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-rivers_05242024.csv")
# roads_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-roads_05242024.csv")

ecoreg_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-US_L3NAME_num.csv")
# rivers_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-rivers_05242024.csv")
# roads_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-roads_05242024.csv")


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

# Write results
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Pvalue-Zscore-US_L3NAME_num", row.names = FALSE)


##########################################
# RIVERS
##########################################
# RETRIEVE DATA
rm(list=ls())
rivers_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-rivers.csv")
rivers_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-rivers.csv")

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

# Write results
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Pvalue-Zscore-rivers", row.names = FALSE)


##########################################
# ROADS
##########################################
# RETRIEVE DATA
rm(list=ls())
roads_null <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Permutation-roads.csv")
roads_1Fst <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-roads.csv")

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

# Write results
write.csv(results, "Fst/Rscripts_local/Routput/DATA-Fst-Pvalue-Zscore-roads", row.names = FALSE)

















# ############################################################################################################
# # SANBOX: TESTING ONE AT A TIME
# ############################################################################################################
# null_distribution <- ecoreg_null %>% # normally distributed null (permutated)
#   filter(paircombo == 15) %>%
#   pull(Fst)

# obs_stat <- ecoreg_1Fst %>%
#   filter(paircombo == 15) %>% # filters rows
#   pull(Fst) %>% # grabs Fst column from filtered df
#   .[1] # selects first element of the vector (i.e. row 1 from previous filter step)

# ##########################################
# # GENERATE P-VALUES
# ##########################################
# # Calculate p-value (for a one-tailed test, assuming a positive effect is of interest)
# p_value <- sum(null_distribution >= obs_stat) / length(null_distribution)
# #sum(c(1,2,3,4,5) >= 3)

# # Calculate pseudo-confidence interval from the null distribution
# conf_interval <- quantile(null_distribution, probs = c(0.025, 0.975))

# # Plotting
# hist(null_distribution, breaks = 30, col = 'skyblue', main = "Null Distribution with Observed Statistic")
# abline(v = obs_stat, col = "red", lwd = 2, lty = 2) # Observed statistic
# abline(v = conf_interval, col = "darkblue", lwd = 2) # Confidence interval from null

# # Display results
# cat("P-value:", p_value, "\n")
# cat("95% Confidence Interval from Null:", conf_interval, "\n")

# ##########################################
# # GENERATE Z-SCORES
# ##########################################
# # Z-SCORES 
# # Z = (x-u)/s
# # X is observed test statistic (e.g. vanilla Fst)
# # u is mean of null (permutated) distribution
# # s is standard deviation of the null distribution
# # assumption about population (null) normality is OK because of central limit theorem (more subsampling will evntually result in a normal distribution)

# # NOTE: This might actually be Cohen's D version of a Z-score (since it's divided by SD, unlike a Z-score which is divided by SE)

# # Calculate the mean and standard deviation of the null distribution
# null_mean <- mean(null_distribution)
# null_sd <- sd(null_distribution)

# # Calculate Z-score
# z_score <- (obs_stat - null_mean) / null_sd
# print(z_score)
