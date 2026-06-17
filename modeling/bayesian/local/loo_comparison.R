rm(list=ls())
library(loo)

filepath <- "modeling/bayesian/local/Routput/"
nonspatial <- "02172026_nonspatial/"
spatial <- "02172026_spatial/"

##########################################################
# Generate loo object from fit
##########################################################
# Note: it's better to save because the fit objects are too large if loading at once
## alternatively can rm(), but saving loo objs speeds things up

#------------------------------
# NON-SPATIAL
#------------------------------
# # SINGLE MODELS
# fit_rivers_NS <- readRDS(file = paste0(filepath, nonspatial, "K4_river_fit.rds"))
# loo_fit_rivers_NS <- loo(fit_rivers_NS)
# saveRDS(loo_fit_rivers_NS, file = paste0(filepath, nonspatial, "loo_fit_rivers_NS.rds"))

# fit_roads_NS <- readRDS(file = paste0(filepath, nonspatial, "K4_roads_fit.rds"))
# loo_fit_roads_NS <- loo(fit_roads_NS)
# saveRDS(loo_fit_roads_NS, file = paste0(filepath, nonspatial, "loo_fit_roads_NS.rds"))

# fit_ecoreg_NS <- readRDS(file = paste0(filepath, nonspatial, "K4_ecoreg_fit.rds"))
# loo_fit_ecoreg_NS <- loo(fit_ecoreg_NS)
# saveRDS(loo_fit_ecoreg_NS, file = paste0(filepath, nonspatial, "loo_fit_ecoreg_NS.rds"))

# # DOUBLE MODELS
# fit_rivers_ecoregs_NS <- readRDS(file = paste0(filepath, nonspatial, "K4_rivers_ecoregs_fit.rds"))
# loo_fit_rivers_ecoregs_NS <- loo(fit_rivers_ecoregs_NS)
# saveRDS(loo_fit_rivers_ecoregs_NS, file = paste0(filepath, nonspatial, "loo_fit_rivers_ecoregs_NS.rds"))

# fit_roads_ecoregs_NS <- readRDS(file = paste0(filepath, nonspatial, "K4_roads_ecoregs_fit.rds"))
# loo_fit_roads_ecoregs_NS <- loo(fit_roads_ecoregs_NS)
# saveRDS(loo_fit_roads_ecoregs_NS, file = paste0(filepath, nonspatial, "loo_fit_roads_ecoregs_NS.rds"))

# fit_roads_rivers_NS <- readRDS(file = paste0(filepath, nonspatial, "K4_roads_rivers_fit.rds"))
# loo_fit_roads_rivers_NS <- loo(fit_roads_rivers_NS)
# saveRDS(loo_fit_roads_rivers_NS, file = paste0(filepath, nonspatial, "loo_fit_roads_rivers_NS.rds"))

# # FULL MODEL
# fit_roads_rivers_ecoregs_NS <- readRDS(file = paste0(filepath, nonspatial, "K4_roads_rivers_ecoregs_fit.rds"))
# loo_fit_roads_rivers_ecoregs_NS <- loo(fit_roads_rivers_ecoregs_NS)
# saveRDS(loo_fit_roads_rivers_ecoregs_NS, file = paste0(filepath, nonspatial, "loo_fit_roads_rivers_ecoregs_NS.rds"))

#------------------------------
# SPATIAL
#------------------------------
# # SINGLE MODELS
# fit_rivers <- readRDS(file = paste0(filepath, spatial, "K4_river_fit.rds"))
# loo_fit_rivers <- loo(fit_rivers)
# saveRDS(loo_fit_rivers, file = paste0(filepath, spatial, "loo_fit_rivers.rds"))

# fit_roads <- readRDS(file = paste0(filepath, spatial, "K4_roads_fit.rds"))
# loo_fit_roads <- loo(fit_roads)
# saveRDS(loo_fit_roads, file = paste0(filepath, spatial, "loo_fit_roads.rds"))

# fit_ecoreg <- readRDS(file = paste0(filepath, spatial, "K4_ecoreg_fit.rds"))
# loo_fit_ecoreg <- loo(fit_ecoreg)
# saveRDS(loo_fit_ecoreg, file = paste0(filepath, spatial, "loo_fit_ecoreg.rds"))

# # DOUBLE MODELS
# fit_rivers_ecoregs <- readRDS(file = paste0(filepath, spatial, "K4_rivers_ecoregs_fit.rds"))
# loo_fit_rivers_ecoregs <- loo(fit_rivers_ecoregs)
# saveRDS(loo_fit_rivers_ecoregs, file = paste0(filepath, spatial, "loo_fit_rivers_ecoregs.rds"))

# fit_roads_ecoregs <- readRDS(file = paste0(filepath, spatial, "K4_roads_ecoregs_fit.rds"))
# loo_fit_roads_ecoregs <- loo(fit_roads_ecoregs)
# saveRDS(loo_fit_roads_ecoregs, file = paste0(filepath, spatial, "loo_fit_roads_ecoregs.rds"))

# fit_roads_rivers <- readRDS(file = paste0(filepath, spatial, "K4_roads_rivers_fit.rds"))
# loo_fit_roads_rivers <- loo(fit_roads_rivers)
# saveRDS(loo_fit_roads_rivers, file = paste0(filepath, spatial, "loo_fit_roads_rivers.rds"))

# # FULL MODEL
# fit_roads_rivers_ecoregs <- readRDS(file = paste0(filepath, spatial, "K4_roads_rivers_ecoregs_fit.rds"))
# loo_fit_roads_rivers_ecoregs <- loo(fit_roads_rivers_ecoregs)
# saveRDS(loo_fit_roads_rivers_ecoregs, file = paste0(filepath, spatial, "loo_fit_roads_rivers_ecoregs.rds"))

#------------------------------
# NULLS
#------------------------------
# fit_spatial_only <- readRDS(file = paste0(filepath, "K4_spatial_only_fit.rds"))
# loo_fit_spatial_only <- loo(fit_spatial_only)
# saveRDS(loo_fit_spatial_only, file = paste0(filepath, "loo_fit_spatial_only.rds"))

# fit_intercept_only <- readRDS(file = paste0(filepath, "K4_intercept_only_fit.rds"))
# loo_fit_intercept_only <- loo(fit_intercept_only)
# saveRDS(loo_fit_intercept_only, file = paste0(filepath, "loo_fit_intercept_only.rds"))


##########################################################
# Reload loo objects
##########################################################
rm(list=ls())
library(loo)

filepath <- "modeling/bayesian/local/Routput/"
nonspatial <- "02172026_nonspatial/"
spatial <- "02172026_spatial/"

#------------------------------
# NON-SPATIAL
#------------------------------
loo_fit_rivers_NS <- readRDS(file = paste0(filepath, nonspatial, "loo_fit_rivers_NS.rds"))
loo_fit_roads_NS <- readRDS(file = paste0(filepath, nonspatial, "loo_fit_roads_NS.rds"))
loo_fit_ecoreg_NS <- readRDS(file = paste0(filepath, nonspatial, "loo_fit_ecoreg_NS.rds"))

loo_fit_rivers_ecoregs_NS <- readRDS(file = paste0(filepath, nonspatial, "loo_fit_rivers_ecoregs_NS.rds"))
loo_fit_roads_ecoregs_NS <- readRDS(file = paste0(filepath, nonspatial, "loo_fit_roads_ecoregs_NS.rds"))
loo_fit_roads_rivers_NS <- readRDS(file = paste0(filepath, nonspatial, "loo_fit_roads_rivers_NS.rds"))

loo_fit_roads_rivers_ecoregs_NS <- readRDS(file = paste0(filepath, nonspatial, "loo_fit_roads_rivers_ecoregs_NS.rds"))

#------------------------------
# SPATIAL
#------------------------------
loo_fit_rivers <- readRDS(file = paste0(filepath, spatial, "loo_fit_rivers.rds"))
loo_fit_roads <- readRDS(file = paste0(filepath, spatial, "loo_fit_roads.rds"))
loo_fit_ecoreg <- readRDS(file = paste0(filepath, spatial, "loo_fit_ecoreg.rds"))

loo_fit_rivers_ecoregs <- readRDS(file = paste0(filepath, spatial, "loo_fit_rivers_ecoregs.rds"))
loo_fit_roads_ecoregs <- readRDS(file = paste0(filepath, spatial, "loo_fit_roads_ecoregs.rds"))
loo_fit_roads_rivers <- readRDS(file = paste0(filepath, spatial, "loo_fit_roads_rivers.rds"))

loo_fit_roads_rivers_ecoregs <- readRDS(file = paste0(filepath, spatial, "loo_fit_roads_rivers_ecoregs.rds"))

#------------------------------
# NULL
#------------------------------
loo_fit_spatial_only <- readRDS(file = paste0(filepath, "loo_fit_spatial_only.rds"))
loo_fit_intercept_only <- readRDS(file = paste0(filepath, "loo_fit_intercept_only.rds"))


##########################################################
# Model comparisons
##########################################################
model_comp <- loo_compare(
    loo_fit_rivers_NS, # 1
    loo_fit_roads_NS, # 2
    loo_fit_ecoreg_NS, # 3
    loo_fit_rivers_ecoregs_NS, # 4
    loo_fit_roads_ecoregs_NS, # 5
    loo_fit_roads_rivers_NS, # 6
    loo_fit_roads_rivers_ecoregs_NS, # 7
    loo_fit_rivers, # 8
    loo_fit_roads, # 9
    loo_fit_ecoreg, # 10
    loo_fit_rivers_ecoregs, # 11
    loo_fit_roads_ecoregs, # 12
    loo_fit_roads_rivers, # 13
    loo_fit_roads_rivers_ecoregs, # 14
    loo_fit_spatial_only, # 15
    loo_fit_intercept_only #16
    )

model_comp