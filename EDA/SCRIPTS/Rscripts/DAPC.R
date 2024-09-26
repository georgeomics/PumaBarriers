# DAPC
library(LEA)

library(adegenet)
# data(dapcIllus) # example datasheet
# class(dapcIllus) # example datasheet
# names(dapcIllus) # example datasheet
# x <- dapcIllus$a # example datasheet
# x # example datasheet

# Using own dataset
data1 <- read.csv("FINAL-DATA.csv", header=T)
data2 <- data1[,-(c(1:7,33,34,35,36,37))]

########### GENIND ##########
# convert df to genind
datagen <- df2genind(data2, ploidy=2, ncode=1, ind.names=data1$Puma_ID, pop=data1$Region)
datagen

########### DAPC PREP ##########
# from: https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
grp <- find.clusters(datagen, max.n.clust = 40) # then, retain 200 PCs as in tutorial..., then choose clusters based on that

# look at 'grp' list for various info

# too see have known original groups compare to inferred groups
table(pop(datagen), grp$grp) # Note: 'pop' comes from 'df2genind' function
table.value(table(pop(datagen), grp$grp), col.lab=paste("inf", 1:6),
            row.lab=paste("ori", 1:6)) # inferred vs original groups
  # Note: Inferred groups are absent geographic knowledge (ori. group 1 may not equal inf. group 1), but this comparison could still be useful to determine where there might be substructuring. For example, original group 6 was inferred to one group (6). However, original group 2 was split into several groups, meaning it could be one place additional substructing is occuring.

########### DAPC ANALYSIS ########### 
dapc1 <- dapc(datagen, grp$grp) # best to choose PCs where retaining more leads to little increase in capturing meaningful information. For discriminant functions, if only a small number of clusters, all eigenvalues can be retained (100 is used because it seems if only have, say 5, it will just take max available)

View(dapc1)
# Essentially, the slots ind.coord and grp.coord may contain the "coordinates" (NOT geographic) of the individuals and of the groups used in scatterplots.
# Contributions of the alleles to each discriminant function are stored in the slot var.contr

# Plots (Note: plotting inferred groups, not original)
scatter(dapc1) # Note: DAPC scatterplots are the main result of DAPC
# scatter can also represent a single discriminant function, which is useful when only one of these has been retained (e.g. in the case k = 2). This is achieved by plotting the densities of individuals on a given discriminant function with different colors for different groups:
# From StatsQuest (https://www.youtube.com/watch?v=azXCzI57Yfc), variation in LDA is "scatter". LDA wants to maximize distance between means and minimize scatter/variation.
  # PCA looks at the genes with the most variation
  # LDA tries to maximize the separation of known categories
myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(dapc1,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

##### Plot of PCs #####
myPal <- colorRampPalette(c("blue","gold","red"))
scatter(dapc1, col=transp(myPal(6)), scree.da=FALSE,
        cell=1.5, cex=2, bg="white",cstar=0) # remember, these are inferred reasons that don't necessarily correspond to original regions

# Alleles highlighting groups
set.seed(4)
contrib <- loadingplot(dapc1$var.contr, axis=2,
                       thres=.07, lab.jitter=1)

# temp is a list invisibly returned by loadingplot which contains the most contributing alleles(i.e., contributions above a given threshold – argument threshold).
  # We can look into the highlighted allele's frequencies across groups (e.g. region or year)
freq09 <- tab(genind2genpop(datagen[loc=c("PP09")]),freq=TRUE)
freq25 <- tab(genind2genpop(datagen[loc=c("PP25")]),freq=TRUE)

par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,.1),las=3)
matplot(freq09, pch=c("A","G"), type="b",
        xlab="cluster",ylab="allele frequency", xaxt="n",
        cex=1.5, main="SNP PP09")
axis(side=1, at=1:6, lab=NULL)
matplot(freq25, pch=c("C","T"), type="b", xlab="cluster",
        ylab="allele frequency", xaxt="n", cex=1.5,
        main="SNP PP25")
axis(side=1, at=1:6, lab=NULL)

##### INTERPRETING GROUP MEMEBERSHIPS #####
# Caution: If membership is based on many PCs, there is a risk of overfitting
# Membership probabilities are based on retained discriminant functions
class(dapc1$posterior)
dim(dapc1$posterior)
round(head(dapc1$posterior), 3) # rows correspond to individuals, columns to groups

summary(dapc1)
# The slot assign.per.pop indicates the proportions of successful reassignment (based on the discriminant functions) of individuals to their original clusters. Large values indicate clear-cut clusters, while low values suggest admixed groups.
par(mfrow=c(1,1))
assignplot(dapc1, subset=1:50) # visualize the above info as a heatmap
# blue crosses represent the prior cluster provided to DAPC

# Can also show this similar to STRUCTURE
compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=1, xlab="individuals", col=funky(6))
# We can also have a closer look at a subset of individuals; for instance, for the first 50 individuals:
compoplot(dapc1, subset=1:50, posi="bottomright",
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=2, xlab="individuals", col=funky(6))

# Can also use R to determine most "admixed" individuals (i.e., having no more than 90% probability of membership in a single cluster):
temp <- which(apply(dapc1$posterior, 1, function(e) all(e<0.9)))
temp
compoplot(dapc1, subset=temp, posi="bottomright",
          txt.leg=paste("Cluster", 1:6),
          ncol=2, col=funky(6))












#################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
####### LANDSCAPE GENOMICS WORKSHOP (UNABLE TO REPEAT) ##########
# convert structure file to genotyping matrix
struct2geno("pumadata.str", ploidy = 2, FORMAT = 2, extra.row = 1, extra.column = 2)

# try getting data into DAPC
# library(data.table)
# example_qced <- fread(input = "example_qced.raw", h=T) # not really "qc'd" yet, but using names from workshop
# write.table(example_qced, "example_qced.raw", col.names=T, row.names=F, 
#             quote=FALSE, sep=" ")


library(adegenet)
dapc_input

#K-means analysis on the principal coponents:
# Max K to test in K-means analysis
maxk < 10
grp <- find.cluster(dapc_input, pca.select  )
