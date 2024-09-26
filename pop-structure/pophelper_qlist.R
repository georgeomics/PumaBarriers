##############################################################################################################
# STRUCTURE
##############################################################################################################
# This script is for processing and visualizing output of STRUCTURE, as well as running evanno tests
rm(list=ls())
library(pophelper)
library(ggplot2)
library(gridExtra)

##################################################
# list STRUCTURE output files in character vector
##################################################
sfiles <- list.files("pop-structure/STRUCTURE/coombs/DATA/output", full.names=TRUE) # need 'full.names=TRUE' to get correct path to folder with files
sfiles <- sfiles[sfiles != "pop-structure/STRUCTURE/coombs/DATA/output/seed.txt"] # remove seed file (otherwise can't convert to qlist)

# convert run files (q-matrices) to qlist
slist <- readQ(sfiles, filetype="structure", indlabfromfile = TRUE)
slistK <- alignK(slist)
slist1 <- slist[c(31:33, 51:53)] # test line
slist2 <- alignK(slist[c(31:33, 51:53)]) # test line

# Save slistK for downstream analyses
saveRDS(slistK, file = "pop-structure/STRUCTURE/SCRIPTS/Routput/slistK.rds")
# slistK <- readRDS("slistK.rds") # read file later
  
# some basic summary stats
tabulateQ(qlist=readQ(sfiles))
head(summariseQ(tabulateQ(slistK)))

##################################################
# EVANNO METHOD
##################################################
# Used to estimate the optimal number of K. The summarised runs table output fromsummariseQ() function can be input to evannoMethodStructure(). 
sr1 <- summariseQ(tabulateQ(slistK))
evannoMethodStructure(data=sr1)

# Evanno graphs
# See 'https://www.rdocumentation.org/packages/pophelper/versions/2.3.1/topics/evannoMethodStructure' for options
p <- evannoMethodStructure(data=sr1,
                           exportplot=F,
                           returnplot=T,
                           returndata=F,
                           basesize=12,
                           linesize=0.7,
                           xaxisbreaks=c(1,2,3,4,5,6,7,8,9,10),
                           xaxislabels=c(1,2,3,4,5,6,7,8,9,10))

pdf("pop-structure/STRUCTURE/SCRIPTS/Routput/Rplot_STRUCTURE-evanno-graphs.pdf")
# Actual plots (i.e. arrange the grid of plots)
grid.arrange(p)
# Close the PDF device
dev.off()

# Barplots
# 'alignK' aligns/orders 'slist' names for easy grabbing of files for plotting
slist1 <- alignK(slistK[c(31:33, 51:53)]) # from list, choose runs/reps of same K

# see 'https://www.rdocumentation.org/packages/pophelper/versions/2.3.1/topics/plotQ' for options
# use `splab` to label strip panels, e.g. `splab=c("test1","test2","test3","test4","test5","test6")`
p1 <- plotQ(slist1, imgoutput="join", returnplot=T, exportplot=F, basesize=11)
pdf("pop-structure/STRUCTURE/SCRIPTS/Routput/Rplot_STRUCTURE-barplots.pdf")
grid.arrange(p1$plot[[1]])
dev.off()

##############################################################################################################
# TESS3
##############################################################################################################
# Once the tess3r.obj output of analysis is saved as a 'tess3r.RData' file, you can use the 'tess3r.RData' file to run the following code
load("pop-structure/TESS/coombs/SCRIPTS/Routput/tess3r.RData") # loads tess3.obj into environment
str(tess3.obj)
# str(tess3.obj, max.level = 2)
# Accessing the Q matrix when K=9
# tess3.obj[[9]]$tess3.run[[3]]$Q

# Convert TESS3 R object to pophelper qlist
## Note: seems to aggregate runs which is OK
tlist <- readQTess3(t3list = tess3.obj, progressbar = FALSE)

# align groups
tlistK <- alignK(tlist) 

# Save tlistK for downstream analyses
saveRDS(tlistK, file = "pop-structure/TESS/SCRIPTS/Routput/tlistK.rds")
# tlistK <- readRDS("tlistK.rds") # read file later

# tabulated data
# some basic summary stats
tabQ <- tabulateQ((tlistK))
summariseQ((tabQ))

# tabulateQ((qlist_tess3))
# sr1_tess3 <- summariseQ(tabulateQ(qlist_tess3))

# bar plots (outputs a plot for each K run)
plotQ(tlistK, exportpath=getwd())