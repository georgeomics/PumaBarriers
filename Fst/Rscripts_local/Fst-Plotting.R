rm(list=ls())
library(ggplot2)
library(ggthemes)

# GRAB DATA
Fstvan_roads <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-roads.csv", header=TRUE)
Fstvan_rivers <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-rivers.csv", header=TRUE)
Fstvan_ecoreg <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Vanilla-US_L3NAME_num.csv", header=TRUE)

FstPZ_roads <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Pvalue-Zscore-roads", header=TRUE)
FstPZ_rivers <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Pvalue-Zscore-rivers", header=TRUE)
FstPZ_ecoreg <- read.csv("Fst/Rscripts_local/Routput/DATA-Fst-Pvalue-Zscore-US_L3NAME_num", header=TRUE)

# ORGANIZE DATA
Fst_roads <- cbind(Fstvan_roads, FstPZ_roads)
Fst_roads <- Fst_roads[,-6] # remove extra paircombon column
Fst_roads$patch <- "roads"

Fst_rivers <- cbind(Fstvan_rivers, FstPZ_rivers)
Fst_rivers <- Fst_rivers[,-6] # remove extra paircombon column
Fst_rivers$patch <- "rivers"

Fst_ecoreg <- cbind(Fstvan_ecoreg, FstPZ_ecoreg)
Fst_ecoreg <- Fst_ecoreg[,-6] # remove extra paircombon column
Fst_ecoreg$patch <- "ecoreg"

Fst <- rbind(Fst_roads, Fst_rivers, Fst_ecoreg)



####################
# PLOTTING
####################
# Index for plotting
Fst$index <- seq_along(Fst$z_score)

# Create a new column 'pair_label' for customized x-axis labels
Fst$pair_label <- paste(Fst$pair.a, "v", Fst$pair.b)

# Begin the plot
png("Fst/Rscripts_local/Routput/test_plot.png", width = 8, height = 6, units = 'in', res = 300)
p <- ggplot(Fst, aes(x = factor(index), y = z_score, fill = factor(patch))) +
  geom_col(show.legend = TRUE) +  # Display actual Z-score values as bars
  scale_x_discrete(labels = Fst$pair_label) +  # Set custom labels for x-axis based on 'pair_label'
  labs(x = "Pair Combination", y = "Z-Score") +
  scale_fill_brewer(palette = "Dark2", name = "Patch Type", labels = c("Ecoregions", "Rivers", "Roads")) +  # Customize legend labels
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Adjust label positioning closer to the bars
    legend.title = element_blank(),  # Ensures no text is displayed for any legend title
    plot.title = element_text(hjust = 0.5, margin = margin(t = 20, b = 10)),  # Center the plot title
    panel.grid.major.y = element_line(color = "gray", size = 0.5),  # Add vertical major grid lines
    panel.grid.minor.y = element_line(color = "lightgray", size = 0.25),  # Add vertical minor grid lines
    panel.grid.major.x = element_blank(),  # Ensure no horizontal major grid lines
    panel.grid.minor.x = element_blank(),  # Ensure no horizontal minor grid lines
    panel.border = element_blank(),  # Remove the box around the plot area
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),  # Add more white space around the plot boundaries
    axis.title.x = element_text(margin = margin(t = 20, b = 10)),  # Add space between the plot and x-axis label
    axis.title.y = element_text(margin = margin(r = 20, l = 10))  # Add space between the plot and y-axis label
  )
print(p)
dev.off()




# Begin the plot
png("Fst/Rscripts_local/Routput/test_plot.png", width = 8, height = 6, units = 'in', res = 300)
p <- ggplot(Fst, aes(x = factor(index), y = z_score, color = factor(patch))) +
  geom_point(size = 3, show.legend = TRUE) +  # Display actual Z-score values as points
  geom_segment(aes(xend = factor(index), yend = 0), linetype = "dotted", alpha = 0.5) +  # Draw vertical lines from points to x-axis
  scale_x_discrete(labels = Fst$pair_label) +  # Set custom labels for x-axis based on 'pair_label'
  scale_y_continuous(limits = c(-5, 55), breaks = seq(0, 55, 10), name = "Z-Score") +  # Set y-axis limits and breaks
  labs(x = "Pair Combination", y = "Z-Score") +
  scale_color_brewer(palette = "Dark2", name = "Patch Type", labels = c("Ecoregions", "Rivers", "Roads")) +  # Customize legend labels using color aesthetic
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Adjust label positioning
    legend.title = element_blank(),  # Ensures no text is displayed for any legend title
    plot.title = element_text(hjust = 0.5, margin = margin(t = 20, b = 10)),  # Center the plot title
    panel.grid.major.y = element_line(color = "gray", size = 0.5),  # Add vertical major grid lines
    panel.grid.minor.y = element_line(color = "lightgray", size = 0.25),  # Add vertical minor grid lines
    panel.grid.major.x = element_blank(),  # Ensure no horizontal major grid lines
    panel.grid.minor.x = element_blank(),  # Ensure no horizontal minor grid lines
    panel.border = element_blank(),  # Remove the box around the plot area
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),  # Add more white space around the plot boundaries
    axis.title.x = element_text(margin = margin(t = 20, b = 10)),  # Add space between the plot and x-axis label
    axis.title.y = element_text(margin = margin(r = 20, l = 10))  # Add space between the plot and y-axis label
  )
print(p)
dev.off()
