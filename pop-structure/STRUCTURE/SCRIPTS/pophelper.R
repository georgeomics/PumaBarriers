##############################################################################################################
# This script is for processing and visualizing output of STRUCTURE, as well as running evanno tests
##############################################################################################################
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
  
# some basic summary stats
tabulateQ(qlist=readQ(sfiles))
head(summariseQ(tabulateQ(slist)))

##################################################
# EVANNO METHOD
##################################################
# Used to estimate the optimal number of K. The summarised runs table output fromsummariseQ() function can be input to evannoMethodStructure(). 
sr1 <- summariseQ(tabulateQ(slist))
evannoMethodStructure(data=sr1)

# Evanno graphs
# See 'https://www.rdocumentation.org/packages/pophelper/versions/2.3.1/topics/evannoMethodStructure' for options
p <- evannoMethodStructure(data=sr1,
                           exportplot=F,
                           returnplot=T,
                           returndata=F,
                           basesize=9,
                           linesize=0.7,
                           xaxisbreaks=c(1,2,3,4,5,6,7,8,9,10),
                           xaxislabels=c(1,2,3,4,5,6,7,8,9,10))

png("pop-structure/STRUCTURE/SCRIPTS/Routput/Rplot_STRUCTURE-evanno-graphs_HWEonly.png", width = 8, height = 6, units = "in", res = 300)
# Actual plots (i.e. arrange the grid of plots)
grid.arrange(p)
dev.off()

# Barplots
# 'alignK' aligns/orders 'slist' names for easy grabbing of files for plotting
slist1 <- alignK(slist[c(31:33, 51:53)]) # from list, choose runs/reps of same K

# see 'https://www.rdocumentation.org/packages/pophelper/versions/2.3.1/topics/plotQ' for options
# use `splab` to label strip panels, e.g. `splab=c("test1","test2","test3","test4","test5","test6")`
p1 <- plotQ(slist1, imgoutput="join", returnplot=T, exportplot=F, basesize=11)
pdf("pop-structure/STRUCTURE/SCRIPTS/Routput/Rplot_STRUCTURE-barplots_HWEonlu.pdf")
grid.arrange(p1$plot[[1]])
dev.off()

