library(ggplot2)
library(reshape2)
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
pdf(paste0(output_dir, K_prefix, model, "bar_graph_comparative.pdf"))
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
# Q_pred_df <- as.data.frame(Q_pred_means)
# colnames(Q_pred_df) <- c("Pop1", "Pop2", "Pop3", "Pop4")
# Q_pred_df$Individual <- 1:nrow(Q_pred_df)
# Q_pred_df <- Q_pred_df[, c("Pop1", "Pop2", "Pop3", "Pop4", "Individual")]
# str(Q_pred_df)

library(ggplot2)
library(tidyr)
library(dplyr)

Q_obs_long$Type <- "Observed"
colnames(Q_obs_long) <- c("Individual", "Response_Variable", "Proportion", "Type")
Q_pred_means_long$Type <- "Predicted"
colnames(Q_pred_means_long) <- c("Individual", "Response_Variable", "Proportion", "Type")
obs_vs_pred <- bind_rows(Q_obs_long, Q_pred_means_long)


pdf(paste0(output_dir, K_prefix, model, "bar_graph_STRUCTURE_like.pdf"), width = 10, height = 7)
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










