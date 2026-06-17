########################################################
# MODEL OUTPUT
########################################################
# OUTPUT SUMMARY
library(rethinking)
sink(paste0(output_dir, K_prefix, predA, predB, predC,"precis.txt")) # like 'pdf()' but for text files
precis(fit, depth=3, prob = 0.95)
sink()

########################################################
# MODEL EVALUATION
########################################################
# ======================== LOO VALIDATION CRITERIA ================================
library(loo)
loo_result <- loo(fit, pars = "log_lik")
sink(paste0(output_dir, K_prefix, predA, predB, predC, "loo.txt"))  # like 'pdf()' but for text files
print(loo_result)
sink()

# ======================== PREDICTOR COLLINEARITY ================================
library(car)

cor_mat_predA <- cor(dummy_mat_predA) # make correlation matrix,
# Variance Inflation Factor (VIF)
dummy_df_predA <- as.data.frame(dummy_mat_predA) # convert matrix to df
vif_model_predA <- lm(Q[,1] ~ ., data = dummy_df_predA) # predA with Q[,1] as response and all columns in dummy_df as preds (necessary for VIF)
vif_vals_predA <- vif(vif_model_predA) # get VIFs
sink(paste0(output_dir, K_prefix, predA, predB, predC, "VIF_results_", predA, ".txt"))
print(vif_vals_predA)
sink()

cor_mat_predB <- cor(dummy_mat_predB) # make correlation matrix,
# Variance Inflation Factor (VIF)
dummy_df_predB <- as.data.frame(dummy_mat_predB) # convert matrix to df
vif_model_predB <- lm(Q[,1] ~ ., data = dummy_df_predB) # predB with Q[,1] as response and all columns in dummy_df as preds (necessary for VIF)
vif_vals_predB <- vif(vif_model_predB) # get VIFs
sink(paste0(output_dir, K_prefix, predA, predB, predC, "VIF_results_", predB, ".txt"))
print(vif_vals_predB)
sink()

cor_mat_predC <- cor(dummy_mat_predC) # make correlation matrix,
# Variance Inflation Factor (VIF)
dummy_df_predC<- as.data.frame(dummy_mat_predC) # convert matrix to df
vif_model_predC<- lm(Q[,1] ~ ., data = dummy_df_predC) # predCwith Q[,1] as response and all columns in dummy_df as preds (necessary for VIF)
vif_vals_predC<- vif(vif_model_predC) # get VIFs
sink(paste0(output_dir, K_prefix, predA, predB, predC, "VIF_results_", predC, ".txt"))
print(vif_vals_predC)
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
sink(paste0(output_dir, K_prefix, predA, predB,predC, "std_resids_quantiles_pop", pop, ".txt"))
quantile(std_resids_pop, probs = c(0, 0.25, 0.5, 0.75, 1)) # between -2 and 2 is generally ideal
sink()
# Q-Q plot
pdf(paste0(output_dir, K_prefix, predA, predB, predC, "qqnorm_pop", pop, ".pdf"))
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
sink(paste0(output_dir, K_prefix, predA, predB, predC, "MoransI_pop", pop, ".txt"))
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
pdf(paste0(output_dir, K_prefix, predA, predB, predC, "ppc_mean_Pop1.pdf"))
ppc_stat(y=Q_obs_pop1, yrep=Q_pred_pop1_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "ppc_mean_Pop2.pdf"))
ppc_stat(y=Q_obs_pop2, yrep=Q_pred_pop2_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "ppc_mean_Pop3.pdf"))
ppc_stat(y=Q_obs_pop3, yrep=Q_pred_pop3_corrected, stat = "mean")
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "ppc_mean_Pop4.pdf"))
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
sink(paste0(output_dir, K_prefix, predA, predB, predC, "bayesian_p_values.txt"))  # like 'pdf()' but for text files
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
pdf(paste0(output_dir, K_prefix, predA, predB, predC, "trace_plot_beta0.pdf"))
stan_trace(fit, pars = c("beta0"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()



pdf(paste0(output_dir, K_prefix, predA, predB, predC, "trace_plot_beta_roads.pdf"))
stan_trace(fit, pars = c("beta_roads"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "trace_plot_beta_rivers.pdf"))
stan_trace(fit, pars = c("beta_rivers"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "trace_plot_beta_ecoregs.pdf"))
stan_trace(fit, pars = c("beta_ecoregs"), nrow = 6, ncol = 4) +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) 
dev.off()


pdf(paste0(output_dir, K_prefix, predA, predB, predC, "trace_plot_phi.pdf"))
stan_trace(fit, pars = c("phi"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "trace_plot_phi_Dirich.pdf"))
stan_trace(fit, pars = c("phi_Dirich"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "trace_plot_sigmasq.pdf"))
stan_trace(fit, pars = c("sigmasq"))
dev.off()

graphics.off()

########## Fit plot ##########
pdf(paste0(output_dir, K_prefix, predA, predB, predC, "fit_plot_beta0.pdf"))
plot(fit, pars=c("beta0"))
dev.off()



pdf(paste0(output_dir, K_prefix, predA, predB, predC, "fit_plot_beta_roads.pdf"))
plot(fit, pars=c("beta_roads"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "fit_plot_beta_rivers.pdf"))
plot(fit, pars=c("beta_rivers"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "fit_plot_beta_ecoregs.pdf"))
plot(fit, pars=c("beta_ecoregs"))
dev.off()


pdf(paste0(output_dir, K_prefix, predA, predB, predC, "fit_plot_phi.pdf"))
plot(fit, pars=c("phi"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "fit_plot_phi_Dirich.pdf"))
plot(fit, pars=c("phi_Dirich"))
dev.off()

pdf(paste0(output_dir, K_prefix, predA, predB, predC, "fit_plot_sigmasq.pdf"))
plot(fit, pars=c("sigmasq"))
dev.off()

graphics.off()

########## Pairs plots ##########
# for (i in 1:3) { # beta0 is length 3
#   pdf(paste0(output_dir, K_prefix, predA, predB, predC, "pairs_plot_beta0_", i, ".pdf"), width = 10, height = 10)
#   pairs(fit, pars = paste0("beta0[", i, "]"), las = 1)
#   dev.off()
# }

for (i in 1:3) { # for each row of beta
  param_names <- paste0("beta_roads[", i, ",", 1:5, "]")
  pdf(paste0(output_dir, K_prefix, predA, predB, predC, "pairs_plot_beta_roads_", i, ".pdf"), width = 10, height = 10)
  pairs(fit, pars = param_names, las = 1)
  dev.off()
}

for (i in 1:3) { # for each row of beta
  param_names <- paste0("beta_rivers[", i, ",", 1:4, "]")
  pdf(paste0(output_dir, K_prefix, predA, predB, predC, "pairs_plot_beta_rivers_", i, ".pdf"), width = 10, height = 10)
  pairs(fit, pars = param_names, las = 1)
  dev.off()
}

for (i in 1:3) { # for each row of beta
  param_names <- paste0("beta_ecoregs[", i, ",", 1:4, "]")
  pdf(paste0(output_dir, K_prefix, predA, predB, predC, "pairs_plot_beta_ecoregs_", i, ".pdf"), width = 10, height = 10)
  pairs(fit, pars = param_names, las = 1)
  dev.off()
}

###############################
# BAR PLOTS
###############################
# source(paste0(R_dir, "output_barplots.R"))
