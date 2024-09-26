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

####################
# GETTING K=5 DATA
# Extract data and calculate means
cluster1 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K5_%d_f", i)]][["Cluster1"]]
})
cluster1_means <- rowMeans(cluster1, na.rm = TRUE)

cluster2 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K5_%d_f", i)]][["Cluster2"]]
})
cluster2_means <- rowMeans(cluster2, na.rm = TRUE)

cluster3 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K5_%d_f", i)]][["Cluster3"]]
})
cluster3_means <- rowMeans(cluster3, na.rm = TRUE)

cluster4 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K5_%d_f", i)]][["Cluster4"]]
})
cluster4_means <- rowMeans(cluster4, na.rm = TRUE)

cluster5 <- sapply(1:10, function(i) {
  slistK[[sprintf("outputfile_K5_%d_f", i)]][["Cluster5"]]
})
cluster5_means <- rowMeans(cluster5, na.rm = TRUE)

K5_data <- cbind(cluster1_means, cluster2_means, cluster3_means, cluster4_means, cluster5_means)
STRUCTURE_data_K5 <- as.data.frame(K5_data)
colnames(STRUCTURE_data_K5) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")

# Save to file for downstream analyses
saveRDS(STRUCTURE_data_K5, file = "modeling/input_data/STRUCTURE_data_K5.rds")


###################################################################################################
# ADMIXTURE
###################################################################################################
# Load qlist
rm(list=ls())
alistK <- readRDS("pop-structure/ADMIXTURE/Routput/alistK.rds")

####################
# GETTING K=3 DATA
# Extract data and calculate means
cluster1 <- sapply(1:10, function(i) {
  alistK[[sprintf("pumaplink.3.%d.Q", i)]][["Cluster1"]]
})
cluster1_means <- rowMeans(cluster1, na.rm = TRUE)

cluster2 <- sapply(1:10, function(i) {
  alistK[[sprintf("pumaplink.3.%d.Q", i)]][["Cluster2"]]
})
cluster2_means <- rowMeans(cluster2, na.rm = TRUE)

cluster3 <- sapply(1:10, function(i) {
  alistK[[sprintf("pumaplink.3.%d.Q", i)]][["Cluster3"]]
})
cluster3_means <- rowMeans(cluster3, na.rm = TRUE)

K3_data <- cbind(cluster1_means, cluster2_means, cluster3_means)
ADMIXTURE_data_K3 <- as.data.frame(K3_data)
colnames(ADMIXTURE_data_K3) <- c("Pop1", "Pop2", "Pop3")

# Save to file for downstream analyses
saveRDS(ADMIXTURE_data_K3, file = "modeling/input_data/ADMIXTURE_data_K3.rds")

####################
# GETTING K=5 DATA
# Extract data and calculate means
cluster1 <- sapply(1:10, function(i) {
  alistK[[sprintf("pumaplink.5.%d.Q", i)]][["Cluster1"]]
})
cluster1_means <- rowMeans(cluster1, na.rm = TRUE)

cluster2 <- sapply(1:10, function(i) {
  alistK[[sprintf("pumaplink.5.%d.Q", i)]][["Cluster2"]]
})
cluster2_means <- rowMeans(cluster2, na.rm = TRUE)

cluster3 <- sapply(1:10, function(i) {
  alistK[[sprintf("pumaplink.5.%d.Q", i)]][["Cluster3"]]
})
cluster3_means <- rowMeans(cluster3, na.rm = TRUE)

cluster4 <- sapply(1:10, function(i) {
  alistK[[sprintf("pumaplink.5.%d.Q", i)]][["Cluster4"]]
})
cluster4_means <- rowMeans(cluster4, na.rm = TRUE)

cluster5 <- sapply(1:10, function(i) {
  alistK[[sprintf("pumaplink.5.%d.Q", i)]][["Cluster5"]]
})
cluster5_means <- rowMeans(cluster5, na.rm = TRUE)

K5_data <- cbind(cluster1_means, cluster2_means, cluster3_means, cluster4_means, cluster5_means)
ADMIXTURE_data_K5 <- as.data.frame(K5_data)
colnames(ADMIXTURE_data_K5) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")

# Save to file for downstream analyses
saveRDS(ADMIXTURE_data_K5, file = "modeling/input_data/ADMIXTURE_data_K5.rds")

###################################################################################################
# TESS
###################################################################################################
# Load qlist
rm(list=ls())
tlistK <- readRDS("pop-structure/TESS/SCRIPTS/Routput/tlistK.rds")

tlistK[["sample3"]][["Cluster1"]]

####################
# GETTING K=3 DATA
# Extract data and calculate means
cluster1 <- tlistK[[sprintf("sample3")]][["Cluster1"]]

cluster2 <- tlistK[[sprintf("sample3")]][["Cluster2"]]

cluster3 <- tlistK[[sprintf("sample3")]][["Cluster3"]]

K3_data <- cbind(cluster1, cluster2, cluster3)
TESS_data_K3 <- as.data.frame(K3_data)
colnames(TESS_data_K3) <- c("Pop1", "Pop2", "Pop3")

# Save to file for downstream analyses
saveRDS(TESS_data_K3, file = "modeling/input_data/TESS_data_K3.rds")


####################
# GETTING K=4 DATA
# Extract data and calculate means
cluster1 <- tlistK[[sprintf("sample4")]][["Cluster1"]]

cluster2 <- tlistK[[sprintf("sample4")]][["Cluster2"]]

cluster3 <- tlistK[[sprintf("sample4")]][["Cluster3"]]

cluster4 <- tlistK[[sprintf("sample4")]][["Cluster4"]]

K4_data <- cbind(cluster1, cluster2, cluster3, cluster4)
TESS_data_K4 <- as.data.frame(K4_data)
colnames(TESS_data_K4) <- c("Pop1", "Pop2", "Pop3", "Pop4")

# Save to file for downstream analyses
saveRDS(TESS_data_K4, file = "modeling/input_data/TESS_data_K4.rds")

####################
# GETTING K=5 DATA
# Extract data and calculate means
cluster1 <- tlistK[[sprintf("sample5")]][["Cluster1"]]

cluster2 <- tlistK[[sprintf("sample5")]][["Cluster2"]]

cluster3 <- tlistK[[sprintf("sample5")]][["Cluster3"]]

cluster4 <- tlistK[[sprintf("sample5")]][["Cluster4"]]

cluster5 <- tlistK[[sprintf("sample5")]][["Cluster5"]]

K5_data <- cbind(cluster1, cluster2, cluster3, cluster4, cluster5)
TESS_data_K5 <- as.data.frame(K5_data)
colnames(TESS_data_K5) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")

# Save to file for downstream analyses
saveRDS(TESS_data_K5, file = "modeling/input_data/TESS_data_K5.rds")

