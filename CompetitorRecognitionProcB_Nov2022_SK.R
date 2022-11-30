# SK 20/06/2022

########################################
########################################

## COMPETITOR RECOGNITION
## Test hypotheses from Peiman & Robinson (2010)

########################################
########################################

# 3 main hypotheses from P&R (2010):
# 1. individuals are less willing to approach a heterospecific (proximity)
# 2. signalling (0.1 chase = signalling) is more common amongst conspecifics
# 3. attacks are more variable between heterospecifics (chase distance)

# Additional questions:
# 4.	Is the intensity of competitor recognition predictable 
#     by phylogenetic relatedness?
# 5.	Is competitor recognition maintained after mass coral mortality?


####################################################################
####################################################################

## 1. UPLOAD EXISTING DATA & PREP WORKSPACE

####################################################################
####################################################################

rm(list=ls())
setwd("/Users/sallykeith/Documents/Dropbox/DataAnalysis/ButterflyfishBiogeography/UncertaintyHypothesis/")
library(EnvStats)
library(vegan)
library(ggplot2)
library(cvequality)
library(gmodels)
library(dplyr)
library(tidyr)
library(DescTools)
library(effects)
library(performance)
library(tidyverse)
library(AICcmodavg)
library(lme4)
library(coin)
library(pscl)


# data on pairwise interactions and species traits including phylogenetic relatedness
# includes conspecific and heterospecific, before and after bleaching
load("DataForCompetitorRecognition2.RData")
# load phylogenetic data
load("PairwisePhylogeneticDistance.RData")
bodysize <- read.csv("SppBodysize.csv")

par(mfcol=c(2,2))

##############

###################################
# PREPARE PHYLOGENETIC DATA
###################################
 
# phylogenetic relatedness data
phylo.dist$spi <- rownames(phylo.dist)
phylo.dist <- phylo.dist[,c(ncol(phylo.dist),2:(ncol(phylo.dist)-1))]
phylo <- data.frame(pivot_longer(phylo.dist,colnames(phylo.dist[2:ncol(phylo.dist)])))
phylo$ij <- paste(phylo[,1],phylo[,2],sep=".")
phylo$ji <- paste(phylo[,2],phylo[,1],sep=".")
phylo <- phylo[,3:5]
colnames(phylo)[1] <- "phylo"


####################################################################
####################################################################

# 2. DO CONSPECIFICS APPROACH EACH OTHER MORE CLOSELY
#    THAN HETEROSPECIFICS?

####################################################################
####################################################################

##############################################################################

# 2.1. DATA PREP FOR
#      (1) MINIMUM OF 5 SAMPLES (encounters) ACROSS 5 INDIVIDUALS
#      (2) MATCHED PAIRS - SPECIES PRESENT IN BOTH CONSPECIFIC AND HETEROSPECIFIC DATA
#      N.B. this iterative process that needs some MANUAL input (eurgh!)

##############################################################################

###################################################
# substitute proximity for ordinal numeric category
# NOTE: requires full dataset i.e., not only aggressive encounters

# before
pre.intra$minimum.distance.nearest25cm. <- gsub("0-25","1",pre.intra$minimum.distance.nearest25cm.)
pre.intra$minimum.distance.nearest25cm. <- gsub("25-50","2",pre.intra$minimum.distance.nearest25cm.)
pre.intra$minimum.distance.nearest25cm. <- gsub("50-75","3",pre.intra$minimum.distance.nearest25cm.)
pre.intra$minimum.distance.nearest25cm. <- gsub("75-100","4",pre.intra$minimum.distance.nearest25cm.)
pre.intra$minimum.distance.nearest25cm. <- as.numeric(pre.intra$minimum.distance.nearest25cm.)

pre.inter$minimum.distance.nearest25cm. <- gsub("0-25","1",pre.inter$minimum.distance.nearest25cm.)
pre.inter$minimum.distance.nearest25cm. <- gsub("25-50","2",pre.inter$minimum.distance.nearest25cm.)
pre.inter$minimum.distance.nearest25cm. <- gsub("50-75","3",pre.inter$minimum.distance.nearest25cm.)
pre.inter$minimum.distance.nearest25cm. <- gsub("75-100","4",pre.inter$minimum.distance.nearest25cm.)
pre.inter$minimum.distance.nearest25cm. <- as.numeric(pre.inter$minimum.distance.nearest25cm.)

# after
post.intra$minimum.distance.nearest25cm. <- gsub("0-25","1",post.intra$minimum.distance.nearest25cm.)
post.intra$minimum.distance.nearest25cm. <- gsub("25-50","2",post.intra$minimum.distance.nearest25cm.)
post.intra$minimum.distance.nearest25cm. <- gsub("50-75","3",post.intra$minimum.distance.nearest25cm.)
post.intra$minimum.distance.nearest25cm. <- gsub("75-100","4",post.intra$minimum.distance.nearest25cm.)
post.intra$minimum.distance.nearest25cm. <- as.numeric(post.intra$minimum.distance.nearest25cm.)

post.inter$minimum.distance.nearest25cm. <- gsub("0-25","1",post.inter$minimum.distance.nearest25cm.)
post.inter$minimum.distance.nearest25cm. <- gsub("25-50","2",post.inter$minimum.distance.nearest25cm.)
post.inter$minimum.distance.nearest25cm. <- gsub("50-75","3",post.inter$minimum.distance.nearest25cm.)
post.inter$minimum.distance.nearest25cm. <- gsub("75-100","4",post.inter$minimum.distance.nearest25cm.)
post.inter$minimum.distance.nearest25cm. <- as.numeric(post.inter$minimum.distance.nearest25cm.)

# view data as barplots
par(mfcol=c(2,2))
barplot(table(pre.intra$minimum.distance.nearest25cm.),main="Conspecific")
barplot(table(pre.inter$minimum.distance.nearest25cm.),main="Heterospecific")
barplot(table(post.intra$minimum.distance.nearest25cm.),main="Conspecific")
barplot(table(post.inter$minimum.distance.nearest25cm.),main="Heterospecific")
####

########################
# (1) MINIMUM 5 SAMPLES
########################

######################
# CONSPECIFICS
######################
csprox <- pre.intra
sort(unique(as.character(csprox$focal.species)))

# BEFORE
subset(as.data.frame(table(csprox$focal.species)),Freq>4)
# MANUAL CHECK - select only those species with >4 encounters
csprox.pairs5 <- c("auriga","baronessa","citrinellus","guttatissimus","kleinii",
                   "lunula","lunulatus","plebeius","rafflesii","trifascialis","unimaculatus")
csprox.sample5 <- csprox[csprox$focal.species %in% csprox.pairs5,]
csprox.sample5$focal.species <- factor(csprox.sample5$focal.species)
table(csprox.sample5$focal.species) # check

# MANUAL CHECK there are at least 5 focal INDIVIDUALS too
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="auriga",]$fish.id)),Freq>0))         # no
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="baronessa",]$fish.id)),Freq>0))      # yes
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="citrinellus",]$fish.id)),Freq>0))    # yes
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="guttatissimus",]$fish.id)),Freq>0))  # yes
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="kleinii",]$fish.id)),Freq>0))        # yes
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="lunula",]$fish.id)),Freq>0))         # yes
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="lunulatus",]$fish.id)),Freq>0))      # yes
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="plebeius",]$fish.id)),Freq>0))       # no
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="rafflesii",]$fish.id)),Freq>0))      # yes
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="trifascialis",]$fish.id)),Freq>0))   # yes
nrow(subset(as.data.frame(table(csprox.sample5[csprox.sample5$focal.species=="unimaculatus",]$fish.id)),Freq>0))   # yes
# MANUAL CHECK - Remove auriga and plebeius
csprox.sample5 <- csprox.sample5[!csprox.sample5$focal.species=="auriga",]
csprox.sample5 <- csprox.sample5[!csprox.sample5$focal.species=="plebeius",]
csprox.sample5$focal.species <- factor(csprox.sample5$focal.species)

# AFTER
csprox.post <- post.intra
# remove Philippines because there was no coral mortality
csprox.post <- csprox.post[!csprox.post$region=="Philippines",]
sort(unique(as.character(csprox.post$focal.species)))
subset(as.data.frame(table(csprox.post$focal.species)),Freq>4)

# MANUAL CHECK - select only those species with >4 encounters
csprox.post.pairs5 <- c("citrinellus","collare","lunulatus","meyeri",
                  "trifasciatus","trifascialis","unimaculatus")
csprox.post.sample5 <- csprox.post[csprox.post$focal.species %in% csprox.post.pairs5,]
csprox.post.sample5$focal.species <- factor(csprox.post.sample5$focal.species)
table(csprox.post.sample5$focal.species) # check

# MANUAL CHECK there are at least 5 focal INDIVIDUALS too
nrow(subset(as.data.frame(table(csprox.post.sample5[csprox.post.sample5$focal.species=="citrinellus",]$fish.id)),Freq>0))     # yes
nrow(subset(as.data.frame(table(csprox.post.sample5[csprox.post.sample5$focal.species=="collare",]$fish.id)),Freq>0))         # yes
nrow(subset(as.data.frame(table(csprox.post.sample5[csprox.post.sample5$focal.species=="lunulatus",]$fish.id)),Freq>0))       # yes
nrow(subset(as.data.frame(table(csprox.post.sample5[csprox.post.sample5$focal.species=="meyeri",]$fish.id)),Freq>0))          # yes
nrow(subset(as.data.frame(table(csprox.post.sample5[csprox.post.sample5$focal.species=="trifascialis",]$fish.id)),Freq>0))    # yes
nrow(subset(as.data.frame(table(csprox.post.sample5[csprox.post.sample5$focal.species=="trifasciatus",]$fish.id)),Freq>0))    # yes
nrow(subset(as.data.frame(table(csprox.post.sample5[csprox.post.sample5$focal.species=="unimaculatus",]$fish.id)),Freq>0))    # yes

# MANUAL CHECK - all species have >4 focal individuals
csprox.post.sample5$focal.species <- factor(csprox.post.sample5$focal.species)

####


######################
# HETEROSPECIFICS
######################
hsprox <- pre.inter
hsprox.post <- post.inter
sort(unique(c(as.character(hsprox$focal.species),as.character(hsprox$encountered.species))))
sort(unique(c(as.character(hsprox.post$focal.species),as.character(hsprox.post$encountered.species))))


# before and after data prepared together (legacy!)
hsprox.pairs <- paste(hsprox$focal.species,hsprox$encountered.species,sep=".")
hsprox.pairs.id <- data.frame(hsprox.pairs,rep(NA,length(hsprox.pairs)))
# remove Philippines because there was no coral mortality
hsprox.post <- hsprox.post[!hsprox.post$region=="Philippines",]
hsprox.post.pairs <- paste(hsprox.post$focal.species,hsprox.post$encountered.species,sep=".")
hsprox.post.pairs.id <- data.frame(hsprox.post.pairs,rep(NA,length(hsprox.post.pairs)))

# need to recognise that species pair can be flipped
# get unique species names
sp.names <- sort(unique(c(as.character(hsprox$focal.species),as.character(hsprox$encountered.species),
               as.character(hsprox.post$focal.species),as.character(hsprox.post$encountered.species))))

# get unique species pair combinations (so not duplicated by being written two different
# ways around)
for(spi in 1:length(sp.names)){
   for(spj in 1:length(sp.names)){
      comb1 <- paste(sp.names[spi],sp.names[spj],sep=".")
      comb2 <- paste(sp.names[spj],sp.names[spi],sep=".")
      pair.id <- which(hsprox.pairs==comb1 | hsprox.pairs==comb2)
      pair.id.post <- which(hsprox.post.pairs==comb1 | hsprox.post.pairs==comb2)
      hsprox.pairs.id[pair.id,2] <- comb1
      hsprox.post.pairs.id[pair.id.post,2] <- comb1
   }
}

# bind column with non-duplicated pair names to full heterospecific data
hsprox.pairs2 <- cbind(hsprox,hsprox.pairs.id[,2])
colnames(hsprox.pairs2)[23] <- "pair.id"
hsprox.post.pairs2 <- cbind(hsprox.post,hsprox.post.pairs.id[,2])
colnames(hsprox.post.pairs2)[23] <- "pair.id"

# find how many samples exist for each species pair.
# note that this does not account for repeated measures within focal individuals
pair.id.samples <- subset(as.data.frame(table(hsprox.pairs2$pair.id)),Freq>4)
post.pair.id.samples <- subset(as.data.frame(table(hsprox.post.pairs2$pair.id)),Freq>4)

# reduce heterospecific aggression datasets to species pairs
# with at least 5 samples
hsprox.sample5 <- hsprox.pairs2[which(hsprox.pairs2$pair.id %in% pair.id.samples[,1]),]
hsprox.post.sample5 <- hsprox.post.pairs2[which(hsprox.post.pairs2$pair.id %in% post.pair.id.samples[,1]),]
# MANUAL CHECK to ensure all species pairs with <5 samples have been removed
as.data.frame(table(as.character(hsprox.sample5$pair.id)))
as.data.frame(table(as.character(hsprox.post.sample5$pair.id)))
# All have at least 5 samples

#####


########################################################################################
# FINAL SAMPLE SIZES USED IN PROXIMITY ANALYSIS FOR MINIMUM 5 SAMPLES, 5 INDIVIDUALS
########################################################################################
FIGURE1.cs.prox.min5.sample.size <- as.data.frame(table(csprox.sample5$focal.species))
FIGURE1.cs.prox.post.min5.sample.size <- as.data.frame(table(csprox.post.sample5$focal.species))
FIGURE1.hs.prox.min5.sample.size <- as.data.frame(table(as.character(hsprox.sample5$pair.id)))
FIGURE1.hs.prox.post.min5.sample.size <- as.data.frame(table(as.character(hsprox.post.sample5$pair.id)))

#####

########################
# (2) MATCHED SPECIES
########################

# BEFORE
csmatch <- as.character(unique(csprox.sample5$focal.species))
hsmatch <- unique(c(as.character(hsprox.sample5$focal.species),as.character(hsprox.sample5$encountered.species)))
csprox.match <- csprox.sample5[csprox.sample5$focal.species %in% hsmatch,]
hsprox.match <- hsprox.sample5[hsprox.sample5$focal.species %in% csprox.match$focal.species,]
hsprox.match <- hsprox.match[hsprox.match$encountered.species %in% csprox.match$focal.species,]
csprox.match$focal.species <- factor(csprox.match$focal.species)
hsprox.match$focal.species <- factor(hsprox.match$focal.species)
hsprox.match$encountered.species <- factor(hsprox.match$encountered.species)
# check species are matched
sort(unique(as.character(csprox.match$focal.species)))
sort(unique(c(as.character(hsprox.match$focal.species),as.character(hsprox.match$encountered.species))))

# AFTER
cspostmatch <- as.character(unique(csprox.post.sample5$focal.species))
hspostmatch <- unique(c(as.character(hsprox.post.sample5$focal.species),as.character(hsprox.post.sample5$encountered.species)))
csprox.post.match <- csprox.post.sample5[csprox.post.sample5$focal.species %in% hspostmatch,]
hsprox.post.match <- hsprox.post.sample5[hsprox.post.sample5$focal.species %in% csprox.post.match$focal.species,]
hsprox.post.match <- hsprox.post.match[hsprox.post.match$encountered.species %in% csprox.post.match$focal.species,]
csprox.post.match <- csprox.post.sample5[csprox.post.sample5$focal.species %in% c(as.character(hsprox.post.match$focal.species),
                                                                                 as.character(hsprox.post.match$encountered.species)),]
csprox.post.match$focal.species <- factor(csprox.post.match$focal.species)
hsprox.post.match$focal.species <- factor(hsprox.post.match$focal.species)
hsprox.post.match$encountered.species <- factor(hsprox.post.match$encountered.species)
# check species are matched
sort(unique(as.character(csprox.post.match$focal.species)))
sort(unique(c(as.character(hsprox.post.match$focal.species),as.character(hsprox.post.match$encountered.species))))



# keep only species pairs in heterospecific datasets (before/after)
# that match species found in conspecific datasets

# BEFORE
hsprox.match5 <- hsprox.sample5[hsprox.sample5$focal.species %in% csprox.match$focal.species & 
                                   hsprox.sample5$encountered.species %in% csprox.match$focal.species,]
# MANUAL check only matched species are included
sort(unique(as.character(csprox.match$focal.species)))
hsprox.match5$pair.id <- factor(hsprox.match5$pair.id)
unique(hsprox.match5$pair.id)

# AFTER
hsprox.post.match5 <- hsprox.post.sample5[hsprox.post.sample5$focal.species %in% csprox.post.match$focal.species & 
                                   hsprox.post.sample5$encountered.species %in% csprox.post.match$focal.species,]
# MANUAL check only matched species are included
sort(unique(as.character(csprox.post.match$focal.species)))
hsprox.post.match5$pair.id <- factor(hsprox.post.match5$pair.id)
unique(hsprox.post.match5$pair.id)
#####

########################################################################################
# FINAL SAMPLE SIZES USED IN PROXIMITY ANALYSIS FOR MATCHED SPECIES
# AND COMBINED SAMPLE SIZES FOR PROXIMITY
########################################################################################
FIGURE1.cs.prox.match.sample.size <- as.data.frame(table(csprox.match$focal.species))
FIGURE1.cs.prox.post.match.sample.size <- as.data.frame(table(csprox.post.match$focal.species))
FIGURE1.hs.prox.match.sample.size <- as.data.frame(table(as.character(hsprox.match$pair.id)))
FIGURE1.hs.prox.post.match.sample.size <- as.data.frame(table(as.character(hsprox.post.match$pair.id)))

proximity.sample.sizes <- list(FIGURE1.cs.prox.min5.sample.size,FIGURE1.cs.prox.post.min5.sample.size,
                                 FIGURE1.hs.prox.min5.sample.size,FIGURE1.hs.prox.post.min5.sample.size,
                          FIGURE1.cs.prox.match.sample.size,FIGURE1.cs.prox.post.match.sample.size,
                          FIGURE1.hs.prox.match.sample.size,FIGURE1.hs.prox.post.match.sample.size)
names(proximity.sample.sizes) <- c("cs.prox.min5.sample.size","cs.prox.post.min5.sample.size",
                                 "hs.prox.min5.sample.size","hs.prox.post.min5.sample.size",
                                 "cs.match.sample.size","cs.post.match.sample.size",
                          "hs.match.sample.size","hs.post.match.sample.size")

# number of species 5 samples
length(unique(c(as.character(csprox.sample5$focal.species),
         as.character(hsprox.sample5$focal.species),
         as.character(hsprox.sample5$encountered.species))))
length(unique(c(as.character(csprox.post.sample5$focal.species),
         as.character(hsprox.post.sample5$focal.species),
         as.character(hsprox.post.sample5$encountered.species))))
# number of species matched
length(unique(c(as.character(csprox.match$focal.species),
         as.character(hsprox.match$focal.species),
         as.character(hsprox.match$encountered.species))))
length(unique(c(as.character(csprox.post.match$focal.species),
         as.character(hsprox.post.match$focal.species),
         as.character(hsprox.post.match$encountered.species))))

# summary of sample size data
lapply(proximity.sample.sizes,summary)

summary(rbind(FIGURE1.cs.prox.min5.sample.size,FIGURE1.hs.prox.min5.sample.size))
summary(rbind(FIGURE1.cs.prox.post.min5.sample.size,FIGURE1.hs.prox.post.min5.sample.size))
summary(rbind(FIGURE1.cs.prox.match.sample.size,FIGURE1.hs.prox.match.sample.size))
summary(rbind(FIGURE1.cs.prox.post.match.sample.size,FIGURE1.hs.prox.post.match.sample.size))

sum(rbind(FIGURE1.cs.prox.min5.sample.size,FIGURE1.hs.prox.min5.sample.size)[,2])
sum(rbind(FIGURE1.cs.prox.post.min5.sample.size,FIGURE1.hs.prox.post.min5.sample.size)[,2])
sum(rbind(FIGURE1.cs.prox.match.sample.size,FIGURE1.hs.prox.match.sample.size)[,2])
sum(rbind(FIGURE1.cs.prox.post.match.sample.size,FIGURE1.hs.prox.post.match.sample.size)[,2])

#####

#############################################################################
# view data as barplots for final check there is nothing strange going on
#############################################################################
par(mfcol=c(2,4))
barplot(table(csprox.sample5$minimum.distance.nearest25cm.),main="Conspecific before 5")
barplot(table(hsprox.sample5$minimum.distance.nearest25cm.),main="Heterospecific before 5")
barplot(table(csprox.post.sample5$minimum.distance.nearest25cm.),main="Conspecific after 5")
barplot(table(hsprox.post.sample5$minimum.distance.nearest25cm.),main="Heterospecific after 5")

barplot(table(csprox.match$minimum.distance.nearest25cm.),main="Conspecific before match")
barplot(table(hsprox.match$minimum.distance.nearest25cm.),main="Heterospecific before match")
barplot(table(csprox.post.match$minimum.distance.nearest25cm.),main="Conspecific after match")
barplot(table(hsprox.post.match$minimum.distance.nearest25cm.),main="Heterospecific after match")
# looks good, no obvious strange things :)
#####


# ADD BODY SIZE
bodysize <- cbind(bodysize,bodysize[,1])
colnames(bodysize) <- c("focal.species","size","encountered.species")
csprox.sample5 <- merge(csprox.sample5,bodysize[,1:2],by='focal.species')
csprox.post.sample5 <- merge(csprox.post.sample5,bodysize[,1:2],by='focal.species')
hsprox.sample5 <- merge(hsprox.sample5,bodysize[,1:2],by="focal.species")
colnames(hsprox.sample5)[ncol(hsprox.sample5)] <- "size.focal"
hsprox.sample5 <- merge(hsprox.sample5,bodysize[,2:3],by="encountered.species")
colnames(hsprox.sample5)[ncol(hsprox.sample5)] <- "size.enc"
hsprox.post.sample5 <- merge(hsprox.post.sample5,bodysize[,1:2],by="focal.species")
colnames(hsprox.post.sample5)[ncol(hsprox.post.sample5)] <- "size.focal"
hsprox.post.sample5 <- merge(hsprox.post.sample5,bodysize[,2:3],by="encountered.species")
colnames(hsprox.post.sample5)[ncol(hsprox.post.sample5)] <- "size.enc"

#################################################

# 2.2. FUNCTION TO TEST PROXIMITY DURING ENCOUNTER

#################################################

Proximity.Diff <- function(cs.prox.before,hs.prox.before,
                            cs.prox.after,hs.prox.after,title){
   
   # Calculate mean proxmity for each individual
   csprox.pre.mean <- aggregate(minimum.distance.nearest25cm.~fish.id+focal.species,data=cs.prox.before,FUN=mean)
   hsprox.pre.mean <- aggregate(minimum.distance.nearest25cm.~fish.id+pair.id,data=hs.prox.before,FUN=mean)
   csprox.post.mean <- aggregate(minimum.distance.nearest25cm.~fish.id+focal.species,data=cs.prox.after,FUN=mean)
   hsprox.post.mean <- aggregate(minimum.distance.nearest25cm.~fish.id+pair.id,data=hs.prox.after,FUN=mean)
   colnames(csprox.pre.mean) <- c("fish.id","pair.id","mean.proximity")
   colnames(hsprox.pre.mean) <- c("fish.id","pair.id","mean.proximity")
   colnames(csprox.post.mean) <- c("fish.id","pair.id","mean.proximity")
   colnames(hsprox.post.mean) <- c("fish.id","pair.id","mean.proximity")
   
   boxplot(csprox.pre.mean$mean.proximity,hsprox.pre.mean$mean.proximity,
           csprox.post.mean$mean.proximity,hsprox.post.mean$mean.proximity,
           names = c("conspecific","heterospecific","conspecific","heterospecific"),
           col = c("white","white","grey","grey"),ylab="mean proximity (cm)",
           xlab="encounter",cex.lab=1.5,main=title)
   
   legend("topleft",c("before","after"),cex=1.5,fill=c("white","grey"))
   
   # statistical test to see if proximity is significantly less in x than y (one-tailed)
   # Mann Whitney U test from coin package - permutation tests more robust
   
   csprox.pre.df <- data.frame(csprox.pre.mean,rep("conspecific",nrow(csprox.pre.mean)),rep("before",nrow(csprox.pre.mean)))
   hsprox.pre.df <- data.frame(hsprox.pre.mean,rep("heterospecific",nrow(hsprox.pre.mean)),rep("before",nrow(hsprox.pre.mean)))
   csprox.post.df <- data.frame(csprox.post.mean,rep("conspecific",nrow(csprox.post.mean)),rep("after",nrow(csprox.post.mean)))
   hsprox.post.df <- data.frame(hsprox.post.mean,rep("heterospecific",nrow(hsprox.post.mean)),rep("after",nrow(hsprox.post.mean)))
   colnames(csprox.pre.df) <- c("fish.id","pair.id","mean.proximity","encounter","time")
   colnames(hsprox.pre.df) <- c("fish.id","pair.id","mean.proximity","encounter","time")
   colnames(csprox.post.df) <- c("fish.id","pair.id","mean.proximity","encounter","time")
   colnames(hsprox.post.df) <- c("fish.id","pair.id","mean.proximity","encounter","time")
   allprox <- rbind(csprox.pre.df,hsprox.pre.df,csprox.post.df,hsprox.post.df)
   
   # compare conspecific against heterospecific
   print("cs vs hs before Wilcoxon permutation test")
   print(wilcox_test(mean.proximity~encounter,data=allprox[allprox$time=="before",], 
                     alternative="less",distribution="approximate"))
   print("cs vs hs after Wilcoxon permutation test")
   print(wilcox_test(mean.proximity~encounter,data=allprox[allprox$time=="after",], 
                     alternative="less",distribution="approximate"))
   
   # compare different times
   print("before vs after conspecifics Wilcoxon permutation test")
   print(wilcox_test(mean.proximity~time,data=allprox[allprox$encounter=="conspecific",], 
                     alternative="less",distribution="approximate"))
   print("before vs after heterospecifics Wilcoxon permutation test")
   print(wilcox_test(mean.proximity~time,data=allprox[allprox$encounter=="heterospecific",], 
                     alternative="less",distribution="approximate"))

   return("allprox" = allprox)
   
   }
#################################################


#################################################
# Run proximity function.
# Generates results of permutation based 1-tailed 
# Mann Whitney-U 
#################################################

# ONLY SPECIES WITH 5 SAMPLES OR MORE
prox.sample5.out <- Proximity.Diff(csprox.sample5,hsprox.sample5,csprox.post.sample5,hsprox.post.sample5,"Minimum 5 samples")
proxsig.sample5.out <- Proximity.Diff(cssig.sample5,hssig.sample5,cssig.post.sample5,hssig.post.sample5,"Minimum 5 samples")

# plot stacked bar for proportions of encounters with each proximity
prox.sample5.out$round.prox <- round(prox.sample5.out$mean.proximity,0)
prox.plot <- as.data.frame(table(prox.sample5.out$round.prox,prox.sample5.out$encounter,prox.sample5.out$time))
colnames(prox.plot) <- c("proximity","encounter","time","freq")
prox.plot$freq <- prox.plot$freq
prox.plot$label <- paste(prox.plot$time,prox.plot$encounter,sep=" ")
prox.plot$label <- ordered(prox.plot$label, levels=c("before conspecific","before heterospecific","after conspecific","after heterospecific"))
# use log proportion of encounters so it is visually clear
# NB. this still fits the Mann Whitney U stats because that is done on rank values
ggplot(prox.plot, aes(x = label, y = log(freq), fill = proximity)) +
  geom_col(colour = "black", position = "fill") +
  theme(axis.text.x = element_text(size = 16, angle = 45,hjust=1),
        axis.text.y = element_text(size = 16, angle = 0),
        axis.title.x=element_blank(),axis.title.y=element_text(size=18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
        scale_fill_grey(name = "proximity (cm)", labels = c("0-24","25-49","50-74","75-100"),start=1,end=0.25) +
        labs(y="log(proportion of encounters)")


# MATCHED - ONLY SPECIES WITH 5 SAMPLES OR MORE AND FOUND IN CS AND HS
prox.match.out <- Proximity.Diff(csprox.match,hsprox.match,csprox.post.match,hsprox.post.match,"Matched species pairs")
prox.match.out$round.prox <- round(prox.match.out$mean.proximity,0)
prox.match.plot <- as.data.frame(table(prox.match.out$round.prox,prox.match.out$encounter,prox.match.out$time))
colnames(prox.match.plot) <- c("proximity","encounter","time","freq")
prox.match.plot$freq <- prox.match.plot$freq
prox.match.plot$label <- paste(prox.match.plot$time,prox.match.plot$encounter,sep=" ")
prox.match.plot$label <- ordered(prox.match.plot$label, levels=c("before conspecific","before heterospecific","after conspecific","after heterospecific"))
# use log proportion of encounters so it is visually clear
# NB. this still fits the Mann Whitney U stats because that is done on rank values
ggplot(prox.match.plot, aes(x = label, y = log(freq), fill = proximity)) +
  geom_col(colour = "black", position = "fill") +
  theme(axis.text.x = element_text(size = 16, angle = 45,hjust=1),
        axis.text.y = element_text(size = 16, angle = 0),
        axis.title.x=element_blank(),axis.title.y=element_text(size=18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
        scale_fill_grey(name = "proximity (cm)", labels = c("0-24","25-49","50-74","75-100"),start=1,end=0.25) +
        labs(y="log(proportion of encounters)")

# SAVE THESE PLOTS FOR USE IN MANUSCRIPT (SEE END OF SCRIPT)

####



#################################################

# 2.3. REGRESSION ANALYSIS

#      SEPARATE REGRESSION MODELS REQUIRED FOR DIET AND FOR 
#      PHYLOGENETIC RELATEDNESS BECAUSE THERE ARE NOT ENOUGH
#      SAMLPES:PREDICTORS FOR A MODEL THAT INCLUDES ALL AS 
#      PREDICTORS AT ONCE
#      Use minimum 5 samples dataset before bleaching

#################################################

prox.match.phylo <- phylo[(phylo$ij %in% as.character(prox.sample5.out$pair.id)) |
                                    (phylo$ji %in% as.character(prox.sample5.out$pair.id)),]

prox.match.phylo1 <- prox.match.phylo[,c(2,1)]
prox.match.phylo2 <- prox.match.phylo[,c(3,1)]  
colnames(prox.match.phylo1) <- c("pair.id","phylo")
colnames(prox.match.phylo2) <- c("pair.id","phylo")
prox.merge1 <- merge(prox.sample5.out,prox.match.phylo1,all.x=T)
prox.merge2 <- merge(prox.merge1,prox.match.phylo2,all.x=T)
prox.pre.merged <- prox.merge2
prox.pre.merged <- aggregate(mean.proximity~pair.id+encounter+time+phylo,data=prox.pre.merged,FUN=mean)
prox.phylo <- prox.pre.merged[prox.pre.merged$time=="before",]

# add body size difference
hsprox.sample5$size.diff <- abs(hsprox.sample5$size.focal-hsprox.sample5$size.enc)
prox.phylo.size <- merge(prox.phylo,hsprox.sample5[,c(23,26)],by="pair.id")
prox.phylo.size <- aggregate(mean.proximity~pair.id+size.diff+phylo,data=prox.phylo.size,mean)
#####


#################################################
# REGRESSION WITH PHYLOGENETIC RELATEDNESS AS PREDICTOR
############################

###################
# Predictor variable relationships
###################################

plot(prox.phylo.size$size.diff,prox.phylo.size$phylo,xlab="size difference",
     ylab="phylogenetic relatedness",cex.lab=1.5)
lines(lowess(prox.phylo.size$size.diff,prox.phylo.size$phylo,f=2/3))
cor.test(prox.phylo.size$size.diff,prox.phylo.size$phylo,method="spearman")

plot(prox.phylo.size$size.diff,prox.phylo.size$mean.proximity,xlab="size difference",
     ylab="mean proximity",cex.lab=1.5)
lines(lowess(prox.phylo.size$size.diff,prox.phylo.size$mean.proximity,f=2/3))

plot(prox.phylo$phylo.size,prox.phylo$mean.proximity,xlab="phylogenetic distance",
     ylab="mean proximity",cex.lab=1.5)
lines(lowess(prox.phylo.size$phylo,prox.phylo.size$mean.proximity,f=2/3))

pairs(prox.phylo.size[,2:4])

####################

# linear regression
lm.prox.phylo <- lm(mean.proximity~phylo+size.diff,data=prox.phylo.size)
# add in quadratic
prox.phylo$phylo2 <- prox.phylo$phylo^2
prox.phylo2 <- prox.phylo[order(prox.phylo$phylo2), ]
lm.prox.phylo.quad <- lm(mean.proximity~phylo+phylo2,data=prox.phylo2)
# plot QQ plot etc to explore model fit
par(mfcol=c(4,2))
plot(lm.prox.phylo)
plot(lm.prox.phylo.quad)
# What is R-squared value for each?
summary(lm.prox.phylo)
summary(lm.prox.phylo.quad)
# Use AICc to see which model is best fit
AICc(lm.prox.phylo)       # best one
AICc(lm.prox.phylo.quad) 
####


####################################################################
####################################################################

# 3. IS SIGNALLING (chase distance 0.1) MORE COMMON IN CONSPECIFICS? 

####################################################################
####################################################################

# Aggression = cause actual or potential physical harm
# Signalling does not have the potential to cause physical harm

#########################################################

# 3.1. DATA PREP FOR MIN. 5 SAMPLES / MATCHED PAIRS

#########################################################

#######################################
# (1) MINIMUM 5 SAMPLES, 5 INDIVIDUALS
#######################################

######################
# CONSPECIFICS
######################

cssig <- csprox.sample5[csprox.sample5$outcome==1,]
sort(unique(as.character(cssig$focal.species)))
 
# BEFORE
subset(as.data.frame(table(cssig$focal.species)),Freq>4)
# MANUAL CHECK - select only those species with >4 encounters
cssig.pairs5 <- c("baronessa","citrinellus","kleinii",
                   "lunulatus","trifascialis","unimaculatus")
cssig.sample5 <- cssig[cssig$focal.species %in% cssig.pairs5,]
cssig.sample5$focal.species <- factor(cssig.sample5$focal.species)
table(cssig.sample5$focal.species) # check

# MANUAL CHECK there are at least 5 focal INDIVIDUALS too
nrow(subset(as.data.frame(table(cssig.sample5[cssig.sample5$focal.species=="baronessa",]$fish.id)),Freq>0))      # yes
nrow(subset(as.data.frame(table(cssig.sample5[cssig.sample5$focal.species=="citrinellus",]$fish.id)),Freq>0))    # yes
nrow(subset(as.data.frame(table(cssig.sample5[cssig.sample5$focal.species=="kleinii",]$fish.id)),Freq>0))        # yes
nrow(subset(as.data.frame(table(cssig.sample5[cssig.sample5$focal.species=="lunulatus",]$fish.id)),Freq>0))      # yes
nrow(subset(as.data.frame(table(cssig.sample5[cssig.sample5$focal.species=="trifascialis",]$fish.id)),Freq>0))   # yes
nrow(subset(as.data.frame(table(cssig.sample5[cssig.sample5$focal.species=="unimaculatus",]$fish.id)),Freq>0))   # yes
# All species have at least 5 focal individuals
cssig.sample5$focal.species <- factor(cssig.sample5$focal.species)


# AFTER
cssig.post <- csprox.post.sample5[csprox.post.sample5$outcome==1,]
# remove Philippines because there was no coral mortality
cssig.post <- cssig.post[!cssig.post$region=="Philippines",]

subset(as.data.frame(table(cssig.post$focal.species)),Freq>4)
# MANUAL CHECK - select only those species with >4 encounters
cssig.post.pairs5 <- c("lunulatus","meyeri","trifasciatus","trifascialis")
cssig.post.sample5 <- cssig.post[cssig.post$focal.species %in% cssig.post.pairs5,]
cssig.post.sample5$focal.species <- factor(cssig.post.sample5$focal.species)
table(cssig.post.sample5$focal.species) # check

# MANUAL CHECK there are at least 5 focal INDIVIDUALS too
nrow(subset(as.data.frame(table(cssig.post.sample5[cssig.post.sample5$focal.species=="lunulatus",]$fish.id)),Freq>0))       # yes
nrow(subset(as.data.frame(table(cssig.post.sample5[cssig.post.sample5$focal.species=="meyeri",]$fish.id)),Freq>0))          # yes
nrow(subset(as.data.frame(table(cssig.post.sample5[cssig.post.sample5$focal.species=="trifascialis",]$fish.id)),Freq>0))    # yes
nrow(subset(as.data.frame(table(cssig.post.sample5[cssig.post.sample5$focal.species=="trifasciatus",]$fish.id)),Freq>0))    # yes
# MANUAL CHECK - all species have >4 focal individuals
cssig.post.sample5$focal.species <- factor(cssig.post.sample5$focal.species)
####

######################
# HETEROSPECIFICS
######################

hssig.sample5 <- hsprox.sample5[hsprox.sample5$outcome==1,]
sort(unique(c(as.character(hssig.sample5$focal.species),as.character(hssig.sample5$encountered.species))))
hssig.post.sample5 <- hsprox.post.sample5[hsprox.post.sample5$outcome==1,]
sort(unique(c(as.character(hssig.post.sample5$focal.species),as.character(hssig.post.sample5$encountered.species))))

# find how many samples exist for each species pair.
# note that this does not account for repeated measures within focal individuals
pair.id.samples.sig <- subset(as.data.frame(table(hssig.sample5$pair.id)),Freq>4)
post.pair.id.samples.sig <- subset(as.data.frame(table(hssig.post.sample5$pair.id)),Freq>4)

# reduce heterospecific aggression datasets to species pairs
# with at least 5 samples
hssig.sample5 <- hssig.sample5[which(hssig.sample5$pair.id %in% pair.id.samples.sig[,1]),]
hssig.post.sample5 <- hssig.post.sample5[which(hssig.post.sample5$pair.id %in% post.pair.id.samples.sig[,1]),]
# MANUAL CHECK to ensure all species pairs with <5 samples have been removed
as.data.frame(table(as.character(hssig.sample5$pair.id)))
as.data.frame(table(as.character(hssig.post.sample5$pair.id)))

####

########################################################################################
# FINAL SAMPLE SIZES USED IN SIGNALLING ANALYSIS FOR MINIMUM 5 SAMPLES, 5 INDIVIDUALS
########################################################################################
FIGURE1.cs.sig.min5.sample.size <- as.data.frame(table(cssig.sample5$focal.species))
FIGURE1.cs.sig.post.min5.sample.size <- as.data.frame(table(cssig.post.sample5$focal.species))
FIGURE1.hs.sig.min5.sample.size <- as.data.frame(table(as.character(hssig.sample5$pair.id)))
FIGURE1.hs.sig.post.min5.sample.size <- as.data.frame(table(as.character(hssig.post.sample5$pair.id)))
#####

########################
# (2) MATCHED SPECIES
########################

# BEFORE
cssig.match1 <- as.character(unique(cssig.sample5$focal.species))
hssig.match1 <- unique(c(as.character(hssig.sample5$focal.species),as.character(hssig.sample5$encountered.species)))
cssig.match <- cssig.sample5[cssig.sample5$focal.species %in% hssig.match1,]
hssig.match <- hssig.sample5[hssig.sample5$focal.species %in% cssig.match$focal.species,]
hssig.match <- hssig.match[hssig.match$encountered.species %in% cssig.match$focal.species,]
cssig.match$focal.species <- factor(cssig.match$focal.species)
hssig.match$focal.species <- factor(hssig.match$focal.species)
hssig.match$encountered.species <- factor(hssig.match$encountered.species)
# check species are matched
sort(unique(as.character(cssig.match$focal.species)))
sort(unique(c(as.character(hssig.match$focal.species),as.character(hssig.match$encountered.species))))

# AFTER
cssigpost.match1 <- as.character(unique(cssig.post.sample5$focal.species))
hssigpost.match1 <- unique(c(as.character(hssig.post.sample5$focal.species),as.character(hssig.post.sample5$encountered.species)))
cssig.post.match <- cssig.post.sample5[cssig.post.sample5$focal.species %in% hssigpost.match1,]
hssig.post.match <- hssig.post.sample5[hssig.post.sample5$focal.species %in% cssig.post.match$focal.species,]
hssig.post.match <- hssig.post.match[hssig.post.match$encountered.species %in% cssig.post.match$focal.species,]
cssig.post.match <- cssig.post.sample5[cssig.post.sample5$focal.species %in% c(as.character(hssig.post.match$focal.species),
                                                                                 as.character(hssig.post.match$encountered.species)),]
cssig.post.match$focal.species <- factor(cssig.post.match$focal.species)
hssig.post.match$focal.species <- factor(hssig.post.match$focal.species)
hssig.post.match$encountered.species <- factor(hssig.post.match$encountered.species)
# check species are matched
sort(unique(as.character(cssig.post.match$focal.species)))
sort(unique(c(as.character(hssig.post.match$focal.species),as.character(hssig.post.match$encountered.species))))

# keep only species pairs in heterospecific datasets (before/after)
# that match species found in conspecific datasets

# BEFORE
hssig.match5 <- hssig.sample5[hssig.sample5$focal.species %in% cssig.match$focal.species & 
                                   hssig.sample5$encountered.species %in% cssig.match$focal.species,]
# MANUAL check only matched species are included
sort(unique(as.character(cssig.match$focal.species)))
hssig.match5$pair.id <- factor(hssig.match5$pair.id)
unique(hssig.match5$pair.id)

# AFTER
hssig.post.match5 <- hssig.post.sample5[hssig.post.sample5$focal.species %in% cssig.post.match$focal.species & 
                                   hssig.post.sample5$encountered.species %in% cssig.post.match$focal.species,]
# MANUAL check only matched species are included
sort(unique(as.character(cssig.post.match$focal.species)))
hssig.post.match5$pair.id <- factor(hssig.post.match5$pair.id)
unique(hssig.post.match5$pair.id)
#####

########################################################################################
# FINAL SAMPLE SIZES USED IN SIGNALLING ANALYSIS FOR MATCHED SPECIES
# AND COMBINED SAMPLE SIZES FOR SIGNALLING
########################################################################################
FIGURE1.cs.sig.match.sample.size <- as.data.frame(table(cssig.match$focal.species))
FIGURE1.cs.sig.post.match.sample.size <- as.data.frame(table(cssig.post.match$focal.species))
FIGURE1.hs.sig.match.sample.size <- as.data.frame(table(as.character(hssig.match$pair.id)))
FIGURE1.hs.sig.post.match.sample.size <- as.data.frame(table(as.character(hssig.post.match$pair.id)))

signalling.sample.sizes <- list(FIGURE1.cs.sig.min5.sample.size,FIGURE1.cs.sig.post.min5.sample.size,
                                 FIGURE1.hs.sig.min5.sample.size,FIGURE1.hs.sig.post.min5.sample.size,
                          FIGURE1.cs.sig.match.sample.size,FIGURE1.cs.sig.post.match.sample.size,
                          FIGURE1.hs.sig.match.sample.size,FIGURE1.hs.sig.post.match.sample.size)
names(signalling.sample.sizes) <- c("cs.sig.min5.sample.size","cs.sig.post.min5.sample.size",
                                 "hs.sig.min5.sample.size","hs.sig.post.min5.sample.size",
                                 "cs.sig.match.sample.size","cs.sig.post.match.sample.size",
                          "hs.sig.match.sample.size","hs.sig.post.match.sample.size")

# number of species 5 samples
length(unique(c(as.character(cssig.sample5$focal.species),
         as.character(hssig.sample5$focal.species),
         as.character(hssig.sample5$encountered.species))))
length(unique(c(as.character(cssig.post.sample5$focal.species),
         as.character(hssig.post.sample5$focal.species),
         as.character(hssig.post.sample5$encountered.species))))
# number of species matched
length(unique(c(as.character(cssig.match$focal.species),
         as.character(hssig.match$focal.species),
         as.character(hssig.match$encountered.species))))
length(unique(c(as.character(cssig.post.match$focal.species),
         as.character(hssig.post.match$focal.species),
         as.character(hssig.post.match$encountered.species))))

# summary of sample size data
lapply(signalling.sample.sizes,summary)
summary(rbind(FIGURE1.cs.sig.min5.sample.size,FIGURE1.hs.sig.min5.sample.size))
summary(rbind(FIGURE1.cs.sig.post.min5.sample.size,FIGURE1.hs.sig.post.min5.sample.size))
summary(rbind(FIGURE1.cs.sig.match.sample.size,FIGURE1.hs.sig.match.sample.size))
summary(rbind(FIGURE1.cs.sig.post.match.sample.size,FIGURE1.hs.sig.post.match.sample.size))

sum(rbind(FIGURE1.cs.sig.min5.sample.size,FIGURE1.hs.sig.min5.sample.size)[,2])
sum(rbind(FIGURE1.cs.sig.post.min5.sample.size,FIGURE1.hs.sig.post.min5.sample.size)[,2])
sum(rbind(FIGURE1.cs.sig.match.sample.size,FIGURE1.hs.sig.match.sample.size)[,2])
sum(rbind(FIGURE1.cs.sig.post.match.sample.size,FIGURE1.hs.sig.post.match.sample.size)[,2])

#####



######################################

# 3.2. SIGNALLING PROPORTION FUNCTION

######################################

# Need to ensure that each dataset contains a signal event or code doesn't work properly.
# Deal with problem of non-independence of samples because some
# individuals have multiple encounters.
# Can't take mean for this one, instead use resampling that chooses 
# one sample per individual each time.
# For many would be the same every time but sampled from individuals with multiple encounters.

######################################
## FUNCTION:
#  1. takes 1 sample encounter per individual
#  2. records number of signalling and chase events
#  3. if no signalling events, record as 0
#  4. on bootstrapped dataset, performs chi-square & g-test

## OUTPUT: 
#  1. p values and 95% confidence intervals for chi-squared test & g-test 
#           corrected for small sample size
#  2. odds ratio for likelihood of signalling (rather than chase) per encounter
#           for conspecifics relative to heterospecifics
#  3. proportion of aggressive encounters that are signalling for cs/hs
#  4. raw frequencies data for each bootstrap

Signal.Proportion <- function(cssig.data,hssig.data,bootstrap.n,title) {
   
   bootstrap.n <- bootstrap.n
   nsig.cs.bs <- c()
   nchase.cs.bs <- c()
   nsig.hs.bs <- c()
   nchase.hs.bs <- c()
   chisq.bs <- matrix(data=NA,nrow=bootstrap.n,ncol=2)
   g.bs <- matrix(data=NA,nrow=bootstrap.n,ncol=2)
   or.bs <- matrix(data=NA,nrow=bootstrap.n,ncol=5)
   pairs.sig.prop.bs <- data.frame(NA,NA,NA,NA)
   colnames(pairs.sig.prop.bs) <- c("bs","pair.id","encounter","chase.distance.m")
   
   cssig.data <- cssig.data[,c(23,1,3,19)]
   hssig.data <- hssig.data[,c(22,23,4,19)]
   colnames(cssig.data) <- c("fish.id","pair.id","region","chase.distance.m")
   colnames(hssig.data) <- c("fish.id","pair.id","region","chase.distance.m")
   
   for(i in 1:bootstrap.n){
      
      cssig.label <- cbind(cssig.data,rep("conspecific",nrow(cssig.data)))
      colnames(cssig.label)[5] <- "encounter"
      hssig.label <- cbind(hssig.data,rep("heterospecific",nrow(hssig.data)))
      colnames(hssig.label)[5] <- "encounter"
      full <- rbind(cssig.label,hssig.label)
      full$pair.id <- as.character(full$pair.id)
      
      # sample dataframe so only 1 encounter per individual
      # to avoid non-independence from repeated measures. 
      # Next sum number of signalling/chase encounters
      
      # conspecifics
      cssig.label.bs <- aggregate(chase.distance.m~fish.id+pair.id+region+encounter,data=cssig.label,
                             FUN=function(x) sample(x,size=1))
      nchase.cs.bs[i] <- sum(cssig.label.bs$chase.distance.m > 0.1)
      nsig.cs.bs[i] <- nrow(cssig.label.bs) - sum(cssig.label.bs$chase.distance.m > 0.1)
      
      # heterospecifics
      hssig.label.bs <- aggregate(chase.distance.m~fish.id+pair.id+region+encounter,data=hssig.label,
                             FUN=function(x) sample(x,size=1))
      nchase.hs.bs[i] <- sum(hssig.label.bs$chase.distance.m > 0.1)
      nsig.hs.bs[i] <- nrow(hssig.label.bs) - sum(hssig.label.bs$chase.distance.m > 0.1)
      
      # join hs and cs
      sig.sample.bs <- rbind(cssig.label.bs,hssig.label.bs)
      pairs.sig <- aggregate(chase.distance.m~pair.id+encounter,data=sig.sample.bs,
                             FUN=function(x) sum(x == 0.1) / length(x))
      pairs.sig2 <- cbind(rep(i,nrow(pairs.sig)),pairs.sig)
      colnames(pairs.sig2)[1] <- "bs"
      pairs.sig.prop.bs <- rbind(pairs.sig.prop.bs,pairs.sig2) 
      
      # test for significant difference between signal vs chase frequency
      # in cs vs hs (2-way contingency table) with chi squared permutation test
      # from coin package - robust to small sample sizes
      chisq.table <- matrix(data=c(nchase.cs.bs[i],nsig.cs.bs[i],
                                   nchase.hs.bs[i],nsig.hs.bs[i]),nrow=2,ncol=2)
      chisq.table <- as.table(chisq.table)
      chisq.res <- chisq_test(chisq.table,simulate.p.value=T,distribution="approximate",nresample=1000)
      chisq.bs[i,1] <- statistic(chisq.res)
      chisq.bs[i,2] <- pvalue(chisq.res)   
   
   }
   
   # combine raw frequencies into a dataframe for export
   raw.freq <- data.frame(nchase.cs.bs,nsig.cs.bs,nchase.hs.bs,nsig.hs.bs)
   colnames(raw.freq) <- c("conspecific.chase","conspecific.signal",
                           "heterospecific.chase","heterospecific.signal")
   colnames(pairs.sig.prop.bs)[4] <- "sig.prop"
   # plot the proportion of agressive encounters with signalling to
   # identify direction of any significant results (chisq is one-tailed)
   # NULL: whether encounter is conspecifics and heterospecifics is independent of
   #       whether aggressive encounter results in signal or chase
   sig.prop.bs <- data.frame(nsig.cs.bs/(nsig.cs.bs+nchase.cs.bs),
                             nsig.hs.bs/(nsig.hs.bs+nchase.hs.bs))
   colnames(sig.prop.bs) <- c("conspecific","heterospecific")
   pairs.sig.prop.bs <- pairs.sig.prop.bs[-1,]

   chisq.res <- apply(chisq.bs,2,function(x) round(ci(x,confidence=0.95),3))
   colnames(chisq.res) <- c("statistic","p.value")
   
   # return mean and 95% confidence intervals for chi squared test,
   # significant p value means non-independence of variables
   return(list("chisq.res" = chisq.res,"signal.proportion" = sig.prop.bs,
               "raw.frequencies" = raw.freq,"pairs.sig.prop.bs" = pairs.sig.prop.bs))
}
####


######################################

# 3.3. RUN SIGNALLING PROPORTION FUNCTION

######################################

par(mfcol=c(2,2))

#########
# BEFORE
#########
   
# ONLY SPECIES WITH 5 SAMPLES OR MORE
 sig.before5 <- Signal.Proportion(cssig.sample5,hssig.sample5,1000,"Min 5 samples species before")
# convert to dataframe
sig.before5.bs.mean <- aggregate(sig.prop~pair.id+encounter,data=sig.before5$pairs.sig.prop.bs,mean)
boxplot(sig.prop~encounter,data=sig.before5$pairs.sig.prop.bs,main="Min 5 samples species before",
        ylab="signalling proportion",cex.lab=2,cex.axis=1.5,cex.main=2,notch=T)
boxplot(sig.prop~encounter,data=sig.before5.bs.mean,main="Min 5 samples species before",
        ylab="signalling proportion",cex.lab=2,cex.axis=1.5,cex.main=2,ylim=c(0,0.8),notch=T)
# see statistical test results
sig.before5$chisq.res

# ONLY SPECIES WITH 5 SAMPLES OR MORE AND FOUND IN CS AND HS
sig.before.match <- Signal.Proportion(cssig.sample5,hssig.match5,1000,"Matched species before")
# convert to dataframe
sig.before.match.bs.mean <- aggregate(sig.prop~pair.id+encounter,data=sig.before.match$pairs.sig.prop.bs,mean)
boxplot(sig.prop~encounter,data=sig.before.match$pairs.sig.prop.bs,main="Matched species before",
        ylab="signalling proportion",cex.lab=2,cex.axis=1.5,cex.main=2)
boxplot(sig.prop~encounter,data=sig.before.match.bs.mean,main="Matched species before",
        ylab="signalling proportion",cex.lab=2,cex.axis=1.5,cex.main=2,ylim=c(0,0.8))
# see statistical test results
sig.before.match$chisq.res

# What is the median proportion of signalling for cs and hs?
paste("conspecific")
summary(sig.before5.bs.mean[sig.before5.bs.mean$encounter=="conspecific",])
paste("heterospecific")
summary(sig.before5.bs.mean[sig.before5.bs.mean$encounter=="heterospecific",])

boxp.before <- boxplot(sig.prop~encounter,data=sig.before5.bs.mean,main="Min 5 samples species before",
        ylab="signalling proportion",cex.lab=2,cex.axis=1.5,cex.main=2,ylim=c(0,0.8),notch=T)
boxp.before


     
####

#########
# AFTER
#########
   
# ONLY SPECIES WITH 5 SAMPLES OR MORE
sig.after5 <- Signal.Proportion(cssig.post.sample5,hssig.post.sample5,1000,"Min 5 samples species after")
sig.after5.bs.mean <- aggregate(sig.prop~pair.id+encounter,data=sig.after5$pairs.sig.prop.bs,mean)
boxplot(sig.prop~encounter,data=sig.after5$pairs.sig.prop.bs,main="Min 5 samples species after")
boxplot(sig.prop~encounter,data=sig.after5.bs.mean,main="Min 5 samples species after",
        ylab="signalling proportion",cex.lab=2,cex.axis=1.5,cex.main=2,ylim=c(0,0.8))
# see statistical test results
sig.after5$chisq.res
   
# ONLY SPECIES WITH 5 SAMPLES OR MORE AND FOUND IN CS AND HS
sig.after.match <- Signal.Proportion(cssig.post.match,hssig.post.match,1000,"Matched species after")
sig.after.match.bs.mean <- aggregate(sig.prop~pair.id+encounter,data=sig.after.match$pairs.sig.prop.bs,mean)
boxplot(sig.prop~encounter,data=sig.after.match$pairs.sig.prop.bs,main="Matched species after",
        cex.lab=2,cex.axis=1.5,cex.main=2,ylim=c(0,0.8),notch=T)
# see statistical test results
sig.after.match$chisq.res

# What is the median proportion of signalling for cs and hs?
paste("conspecific")
summary(sig.after5.bs.mean[sig.after5.bs.mean$encounter=="conspecific",])
paste("heterospecific")
summary(sig.after5.bs.mean[sig.after5.bs.mean$encounter=="heterospecific",])

boxp.after <- boxplot(sig.prop~encounter,data=sig.after5.bs.mean,main="Min 5 samples species after",
        ylab="signalling proportion",cex.lab=2,cex.axis=1.5,cex.main=2,ylim=c(0,0.8))
boxp.after

######################################

# 3.4. TEST FOR DIFFERENCE IN SIGNALLING 
#      PROPORTION BEFORE/AFTER CORAL MORTALITY

######################################

cs.ba <- cbind(sig.before5$raw.frequencies[,1:2],sig.after5$raw.frequencies[,1:2])
colnames(cs.ba) <- c("bcs.chase","bcs.signal","acs.chase","acs.signal")
hs.ba <- cbind(sig.before5$raw.frequencies[,3:4],sig.after5$raw.frequencies[,3:4])
colnames(hs.ba) <- c("bhs.chase","bhs.signal","ahs.chase","ahs.signal")

# conspecific
chisq.bs.cs.ba <- matrix(data=NA,nrow=nrow(cs.ba),ncol=2)

for(irow in 1:nrow(cs.ba)){
   
   chisq.mat <- matrix(data=c(cs.ba[irow,1],cs.ba[irow,2],cs.ba[irow,3],cs.ba[irow,4]),nrow=2,ncol=2)
   chisq.table <- as.table(chisq.mat)
   chisq.res.ba <- chisq_test(chisq.table,simulate.p.value=T,distribution="approximate",nresample=1000)
   chisq.bs.cs.ba[irow,1] <- statistic(chisq.res.ba)
   chisq.bs.cs.ba[irow,2] <- pvalue(chisq.res.ba)  
}
 
colMeans(chisq.bs.cs.ba)  
ba.cs.out <- apply(chisq.bs.cs.ba,2,function(x) round(ci(x,confidence=0.95),3))
colnames(ba.cs.out) <- c("statistic","p.value")
ba.cs.out

#h eterospecific
chisq.bs.hs.ba <- matrix(data=NA,nrow=nrow(hs.ba),ncol=2)

for(irow in 1:nrow(hs.ba)){
   
   chisq.mat <- matrix(data=c(hs.ba[irow,1],hs.ba[irow,2],hs.ba[irow,3],hs.ba[irow,4]),nrow=2,ncol=2)
   chisq.table <- as.table(chisq.mat)
   chisq.res.ba <- chisq_test(chisq.table,simulate.p.value=T,distribution="approximate",nresample=1000)
   chisq.bs.hs.ba[irow,1] <- statistic(chisq.res.ba)
   chisq.bs.hs.ba[irow,2] <- pvalue(chisq.res.ba)  
}
 
colMeans(chisq.bs.hs.ba)  
ba.hs.out <- apply(chisq.bs.hs.ba,2,function(x) round(ci(x,confidence=0.95),3))
colnames(ba.hs.out) <- c("statistic","p.value")
ba.hs.out

###############################


#################################################

# 3.5. REGRESSION ANALYSIS

#################################################

sig.match.phylo <- phylo[(phylo$ij %in% as.character(sig.before5.bs.mean$pair.id)) |
                           (phylo$ji %in% as.character(sig.before5.bs.mean$pair.id)),]

sig.match.phylo1 <- sig.match.phylo[,c(2,1)]
sig.match.phylo2 <- sig.match.phylo[,c(3,1)]  
colnames(sig.match.phylo1) <- c("pair.id","phylo")
colnames(sig.match.phylo2) <- c("pair.id","phylo")
sig.merge1 <- merge(sig.before5.bs.mean,sig.match.phylo1,all.x=T)
sig.merge2 <- merge(sig.merge1,sig.match.phylo2,all.x=T)
sig.pre.merged <- sig.merge2
# remove NAs
sig.pre.merged <- sig.pre.merged[!is.na(sig.pre.merged$phylo),]
sig.phylo <- aggregate(sig.prop~pair.id+phylo,data=sig.pre.merged,FUN=mean)

# add body size difference
sig.phylo.size <- merge(sig.phylo,hssig.sample5[,c(23,26)],by="pair.id")
sig.phylo.size <- aggregate(sig.prop~pair.id+size.diff+phylo,data=sig.phylo.size,mean)


#################################################
# REGRESSION WITH PHYLOGENETIC RELATEDNESS AS PREDICTOR
################################################

plot(sig.phylo.size$size.diff,sig.phylo.size$phylo,xlab="size difference",
     ylab="phylogenetic relatedness",cex.lab=1.5)
lines(lowess(sig.phylo$phylo,sig.phylo$size.diff,f=2/3))
cor.test(sig.phylo.size$size.diff,sig.phylo.size$phylo,method="spearman")

plot(sig.phylo$phylo,sig.phylo$sig.prop,xlab="phylogenetic relatedness",
     ylab="signalling proportion",cex.lab=1.5)
lines(lowess(sig.phylo$phylo,sig.phylo$sig.prop,f=2/3))

plot(sig.phylo.size$size.diff,sig.phylo.size$sig.prop,xlab="size difference",
     ylab="signalling proportion",cex.lab=1.5)
lines(lowess(sig.phylo.size$size.diff,sig.phylo.size$sig.prop,f=2/3))

# linear regression
glm.sig.phylo <- glm(sig.prop~phylo+size.diff,family=binomial,data=sig.phylo.size)
anova(glm.sig.phylo,test="Chisq")
# add in quadratic
sig.phylo.size$phylo2 <- sig.phylo.size$phylo^2
sig.phylo2.size <- sig.phylo.size[order(sig.phylo.size$phylo2), ]
glm.sig.phylo.quad <- glm(sig.prop~phylo+phylo2+size.diff,family=binomial,data=sig.phylo2.size)
# What is summary for each?
summary(glm.sig.phylo)
summary(glm.sig.phylo.quad)
# Use AICc to see which model is best fit
AICc(glm.sig.phylo)       # LOWEST AICc
AICc(glm.sig.phylo.quad) 
# McFadden R2 index (R2 equivalent for GLM)
pR2(glm.sig.phylo)

#####



####################################################################
####################################################################

# 4. HOW VARIABLE ARE CHASE DISTANCES?

####################################################################
####################################################################

###############################################

# 4.1. DATA PREPARATION

###############################################

# remove extreme outliers from csagg (30 m and 50 m chase)
# group all chases >10 m into one bin
cssig.sample5$chase.distance.m[cssig.sample5$chase.distance.m>10] <- 10
hssig.sample5$chase.distance.m[hssig.sample5$chase.distance.m>10] <- 10


#######################################
# (1) MINIMUM 5 SAMPLES, 5 INDIVIDUALS
#######################################

######################
# CONSPECIFICS
######################

cschase.sample5 <- cssig.sample5[which(cssig.sample5$chase.distance.m>0.1),]
sort(unique(as.character(cschase.sample5$focal.species)))
 
# BEFORE
subset(as.data.frame(table(cschase.sample5$focal.species)),Freq>4)
# MANUAL CHECK - same species as signalling with >4 encounters

# MANUAL CHECK there are at least 5 focal INDIVIDUALS too
nrow(subset(as.data.frame(table(cschase.sample5[cschase.sample5$focal.species=="baronessa",]$fish.id)),Freq>0))      # yes
nrow(subset(as.data.frame(table(cschase.sample5[cschase.sample5$focal.species=="citrinellus",]$fish.id)),Freq>0))    # yes
nrow(subset(as.data.frame(table(cschase.sample5[cschase.sample5$focal.species=="kleinii",]$fish.id)),Freq>0))        # yes
nrow(subset(as.data.frame(table(cschase.sample5[cschase.sample5$focal.species=="lunulatus",]$fish.id)),Freq>0))      # yes
nrow(subset(as.data.frame(table(cschase.sample5[cschase.sample5$focal.species=="trifascialis",]$fish.id)),Freq>0))   # yes
nrow(subset(as.data.frame(table(cschase.sample5[cschase.sample5$focal.species=="unimaculatus",]$fish.id)),Freq>0))   # yes
# All species have at least 5 focal individuals
cschase.sample5$focal.species <- factor(cschase.sample5$focal.species)


# AFTER
cschase.post.sample5 <- cssig.post.sample5[cssig.post.sample5$chase.distance.m>0.1,]
# remove Philippines because there was no coral mortality
cschase.post.sample5 <- cschase.post.sample5[!cschase.post.sample5$region=="Philippines",]

subset(as.data.frame(table(cschase.post.sample5$focal.species)),Freq>4)
# MANUAL CHECK - same species as signalling with >4 encounters

# MANUAL CHECK there are at least 5 focal INDIVIDUALS too
nrow(subset(as.data.frame(table(cschase.post.sample5[cschase.post.sample5$focal.species=="lunulatus",]$fish.id)),Freq>0))       # yes
nrow(subset(as.data.frame(table(cschase.post.sample5[cschase.post.sample5$focal.species=="meyeri",]$fish.id)),Freq>0))          # yes
nrow(subset(as.data.frame(table(cschase.post.sample5[cschase.post.sample5$focal.species=="trifascialis",]$fish.id)),Freq>0))    # yes
nrow(subset(as.data.frame(table(cschase.post.sample5[cschase.post.sample5$focal.species=="trifasciatus",]$fish.id)),Freq>0))    # yes
# MANUAL CHECK - all species have >4 focal individuals

#####

######################
# HETEROSPECIFICS
######################

hschase.sample5 <- hssig.sample5[hssig.sample5$chase.distance>0.1,]
sort(unique(c(as.character(hschase.sample5$focal.species),as.character(hschase.sample5$encountered.species))))
hschase.post.sample5 <- hssig.post.sample5[hssig.post.sample5$chase.distance>0.1,]
sort(unique(c(as.character(hschase.post.sample5$focal.species),as.character(hschase.post.sample5$encountered.species))))

# find how many samples exist for each species pair.
# note that this does not account for repeated measures within focal individuals
pair.id.samples.chase <- subset(as.data.frame(table(hschase.sample5$pair.id)),Freq>4)
post.pair.id.samples.chase <- subset(as.data.frame(table(hschase.post.sample5$pair.id)),Freq>4)

# reduce heterospecific aggression datasets to species pairs
# with at least 5 samples
hschase.sample5 <- hschase.sample5[which(hschase.sample5$pair.id %in% pair.id.samples.chase[,1]),]
hschase.post.sample5 <- hschase.post.sample5[which(hschase.post.sample5$pair.id %in% post.pair.id.samples.chase[,1]),]
# MANUAL CHECK to ensure all species pairs with <5 samples have been removed
as.data.frame(table(as.character(hschase.sample5$pair.id)))
as.data.frame(table(as.character(hschase.post.sample5$pair.id)))

####


########################################################################################
# FINAL SAMPLE SIZES USED IN CHASE ANALYSIS FOR MINIMUM 5 SAMPLES, 5 INDIVIDUALS
########################################################################################
FIGURE1.cs.chase.min5.sample.size <- as.data.frame(table(cschase.sample5$focal.species))
FIGURE1.cs.chase.post.min5.sample.size <- as.data.frame(table(cschase.post.sample5$focal.species))
FIGURE1.hs.chase.min5.sample.size <- as.data.frame(table(as.character(hschase.sample5$pair.id)))
FIGURE1.hs.chase.post.min5.sample.size <- as.data.frame(table(as.character(hschase.post.sample5$pair.id)))
#####

########################
# (2) MATCHED SPECIES
########################

# BEFORE
cschase.match1 <- as.character(unique(cschase.sample5$focal.species))
hschase.match1 <- unique(c(as.character(hschase.sample5$focal.species),as.character(hschase.sample5$encountered.species)))
cschase.match <- cschase.sample5[cschase.sample5$focal.species %in% hschase.match1,]
hschase.match <- hschase.sample5[hschase.sample5$focal.species %in% cschase.match$focal.species,]
hschase.match <- hschase.match[hschase.match$encountered.species %in% cschase.match$focal.species,]
cschase.match$focal.species <- factor(cschase.match$focal.species)
hschase.match$focal.species <- factor(hschase.match$focal.species)
hschase.match$encountered.species <- factor(hschase.match$encountered.species)
# check species are matched
sort(unique(as.character(cschase.match$focal.species)))
sort(unique(c(as.character(hschase.match$focal.species),as.character(hschase.match$encountered.species))))

# AFTER
cschasepost.match1 <- as.character(unique(cschase.post.sample5$focal.species))
hschasepost.match1 <- unique(c(as.character(hschase.post.sample5$focal.species),as.character(hssig.post.sample5$encountered.species)))
cschase.post.match <- cschase.post.sample5[cschase.post.sample5$focal.species %in% hschasepost.match1,]
hschase.post.match <- hschase.post.sample5[hschase.post.sample5$focal.species %in% cschase.post.match$focal.species,]
hschase.post.match <- hschase.post.match[hschase.post.match$encountered.species %in% cschase.post.match$focal.species,]
cschase.post.match <- cschase.post.sample5[cschase.post.sample5$focal.species %in% c(as.character(hschase.post.match$focal.species),
                                                                                 as.character(hschase.post.match$encountered.species)),]
cschase.post.match$focal.species <- factor(cschase.post.match$focal.species)
hschase.post.match$focal.species <- factor(hschase.post.match$focal.species)
hschase.post.match$encountered.species <- factor(hschase.post.match$encountered.species)
# check species are matched
sort(unique(as.character(cschase.post.match$focal.species)))
sort(unique(c(as.character(hschase.post.match$focal.species),as.character(hschase.post.match$encountered.species))))

# keep only species pairs in heterospecific datasets (before/after)
# that match species found in conspecific datasets

# BEFORE
hschase.match5 <- hschase.sample5[hschase.sample5$focal.species %in% cschase.match$focal.species & 
                                   hschase.sample5$encountered.species %in% cschase.match$focal.species,]
# MANUAL check only matched species are included
sort(unique(as.character(cschase.match$focal.species)))
hschase.match5$pair.id <- factor(hschase.match5$pair.id)
unique(hschase.match5$pair.id)

# AFTER
hschase.post.match5 <- hschase.post.sample5[hschase.post.sample5$focal.species %in% cschase.post.match$focal.species & 
                                   hschase.post.sample5$encountered.species %in% cschase.post.match$focal.species,]
# MANUAL check only matched species are included
sort(unique(as.character(cschase.post.match$focal.species)))
hschase.post.match5$pair.id <- factor(hschase.post.match5$pair.id)
unique(hschase.post.match5$pair.id)
####

########################################################################################
# FINAL SAMPLE SIZES USED IN CHASE ANALYSIS FOR MATCHED SPECIES
# AND COMBINED SAMPLE SIZES FOR CHASES
########################################################################################
FIGURE1.cs.chase.match.sample.size <- as.data.frame(table(cschase.match$focal.species))
FIGURE1.cs.chase.post.match.sample.size <- as.data.frame(table(cschase.post.match$focal.species))
FIGURE1.hs.chase.match.sample.size <- as.data.frame(table(as.character(hschase.match$pair.id)))
FIGURE1.hs.chase.post.match.sample.size <- as.data.frame(table(as.character(hschase.post.match$pair.id)))

chase.sample.sizes <- list(FIGURE1.cs.chase.min5.sample.size,FIGURE1.cs.chase.post.min5.sample.size,
                                 FIGURE1.hs.chase.min5.sample.size,FIGURE1.hs.chase.post.min5.sample.size,
                          FIGURE1.cs.chase.match.sample.size,FIGURE1.cs.chase.post.match.sample.size,
                          FIGURE1.hs.chase.match.sample.size,FIGURE1.hs.chase.post.match.sample.size)
names(chase.sample.sizes) <- c("cs.chase.min5.sample.size","cs.chase.post.min5.sample.size",
                                 "hs.chase.min5.sample.size","hs.chase.post.min5.sample.size",
                                 "cs.chase.match.sample.size","cs.chase.post.match.sample.size",
                          "hs.chase.match.sample.size","hs.chase.post.match.sample.size")

# number of species 5 samples
length(unique(c(as.character(cschase.sample5$focal.species),
         as.character(hschase.sample5$focal.species),
         as.character(hschase.sample5$encountered.species))))
length(unique(c(as.character(cschase.post.sample5$focal.species),
         as.character(hschase.post.sample5$focal.species),
         as.character(hschase.post.sample5$encountered.species))))
# number of species matched
length(unique(c(as.character(cschase.match$focal.species),
         as.character(hschase.match$focal.species),
         as.character(hschase.match$encountered.species))))
length(unique(c(as.character(cschase.post.match$focal.species),
         as.character(hschase.post.match$focal.species),
         as.character(hschase.post.match$encountered.species))))

# summary of sample size data
lapply(chase.sample.sizes,summary)
summary(rbind(FIGURE1.cs.chase.min5.sample.size,FIGURE1.hs.chase.min5.sample.size))
summary(rbind(FIGURE1.cs.chase.post.min5.sample.size,FIGURE1.hs.chase.post.min5.sample.size))
summary(rbind(FIGURE1.cs.chase.match.sample.size,FIGURE1.hs.chase.match.sample.size))
summary(rbind(FIGURE1.cs.chase.post.match.sample.size,FIGURE1.hs.chase.post.match.sample.size))

sum(rbind(FIGURE1.cs.chase.min5.sample.size,FIGURE1.hs.chase.min5.sample.size)[,2])
sum(rbind(FIGURE1.cs.chase.post.min5.sample.size,FIGURE1.hs.chase.post.min5.sample.size)[,2])
sum(rbind(FIGURE1.cs.chase.match.sample.size,FIGURE1.hs.chase.match.sample.size)[,2])
sum(rbind(FIGURE1.cs.chase.post.match.sample.size,FIGURE1.hs.chase.post.match.sample.size)[,2])

#####




##########################################################

# 4.2. CALCULATE MEANS FOR EACH SPECIES PAIR WITHIN EACH
#      INDIVIDUAL OBSERVED TO DEAL WITH REPEATED MEASURES

##########################################################

cschase.sample5.mean <- aggregate(chase.distance.m~fish.id+focal.species+region,
                                data=cschase.sample5,FUN=mean)
cschase.post.sample5.mean <- aggregate(chase.distance.m~fish.id+focal.species+region,
                                    data=cschase.post.sample5,FUN=mean)
hschase.sample5.mean <- aggregate(chase.distance.m~fish.id+pair.id+region,
                                data=hschase.sample5,FUN=mean)
hschase.post.sample5.mean <- aggregate(chase.distance.m~fish.id+pair.id+region,
                                     data=hschase.post.sample5,FUN=mean)
# rename columns to ensure they match across data subsets
colnames(cschase.sample5.mean) <- c("fish.id","pair.id","region","chase.distance.m")
colnames(cschase.post.sample5.mean) <- c("fish.id","pair.id","region","chase.distance.m")
colnames(hschase.sample5.mean) <- c("fish.id","pair.id","region","chase.distance.m")
colnames(hschase.post.sample5.mean) <- c("fish.id","pair.id","region","chase.distance.m")

cschase.match.mean <- aggregate(chase.distance.m~fish.id+focal.species+region,
                                data=cschase.match,FUN=mean)
cschase.post.match.mean <- aggregate(chase.distance.m~fish.id+focal.species+region,
                                    data=cschase.post.match,FUN=mean)
hschase.match.mean <- aggregate(chase.distance.m~fish.id+pair.id+region,
                                data=hschase.match,FUN=mean)
hschase.post.match.mean <- aggregate(chase.distance.m~fish.id+pair.id+region,
                                     data=hschase.post.match,FUN=mean)
# rename columns to ensure they match across data subsets
colnames(cschase.match.mean) <- c("fish.id","pair.id","region","chase.distance.m")
colnames(cschase.post.match.mean) <- c("fish.id","pair.id","region","chase.distance.m")
colnames(hschase.match.mean) <- c("fish.id","pair.id","region","chase.distance.m")
colnames(hschase.post.match.mean) <- c("fish.id","pair.id","region","chase.distance.m")

######


#################################################################

# 4.3 MODIFIED FUNCTION TO TEST CHASE DISTANCE 
#     VARIATION FOR EACH SPECIES PAIR SEPARATELY.
#     Idea is to level the playing field on expected
#     variation due to number of different species included

################################################################

Chase.Variation.Pairs <- function(cs.chase,hs.chase){

   
   # boxplot of all species pairs combined 
   boxplot(cs.chase$chase.distance.m,hs.chase$chase.distance.m,
           names=c("conspecific","heterospecific"),
           ylab="chase distance (m)")
  
   # do boxplot with each species pair plotted separately
   cs.label <- cbind(cs.chase,rep("conspecific",nrow(cs.chase)))
   colnames(cs.label)[5] <- "encounter"
   hs.label <- cbind(hs.chase,rep("heterospecific",nrow(hs.chase)))
   colnames(hs.label)[5] <- "encounter"
   full <- rbind(cs.label,hs.label)
   full$pair.id <- as.character(full$pair.id)
   full$boxplot.label <- as.factor(paste(full$encounter,full$pair.id,sep="."))
   
   
   #############################################################
   # calculate coefficient of variation for each species pair
   #############################################################

   encounter.label <- c(rep("conspecific",length(unique(cs.chase$pair.id))),
                        rep("heterospecific",length(unique(hs.chase$pair.id))))
   full.pair.id <- unique(full$pair.id)
   pair.cv.res <- data.frame(encounter.label,full.pair.id,rep(NA,length(full.pair.id)))
   pair.list <- list()
   
   for(pairi in 1:length(full.pair.id)){
   
       pair.subset <- full[full$pair.id==full.pair.id[pairi],]
       pair.cv.res[pairi,3] <- cv(pair.subset$chase.distance.m)
       pair.list[[pairi]] <- pair.subset$chase.distance.m
       names(pair.list)[[pairi]] <- as.character(full.pair.id[pairi])
   
   }

   par(mfcol=c(1,2))
   colnames(pair.cv.res) <- c("encounter","pair.id","CV")
   boxplot(CV~encounter,data=pair.cv.res,ylab="CV across species pairs")
   
   # set up loop to compare each species pair to each.
   # Deals with problem that there might be a higher number of different species pairs 
   # included in heterospecific than conspecific (fewer possible combinations).
   
   mslr.pair.res <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA)
   colnames(mslr.pair.res) <- c("pairi","pairj","n.pairi","n.pairj","cv.pairi","cv.pairj","mslr..stat","mslr.p")
   
   for(pairi in 1:length(full.pair.id)){
      
      mslr.pair.resj <- data.frame(rep(NA,length(full[full$pair.id==full.pair.id[pairi],])),
                        rep(NA,length(full[full$pair.id==full.pair.id[pairi],])),
                        rep(NA,length(full[full$pair.id==full.pair.id[pairi],])),
                        rep(NA,length(full[full$pair.id==full.pair.id[pairi],])),
                        rep(NA,length(full[full$pair.id==full.pair.id[pairi],])),
                        rep(NA,length(full[full$pair.id==full.pair.id[pairi],])),
                        rep(NA,length(full[full$pair.id==full.pair.id[pairi],])),
                        rep(NA,length(full[full$pair.id==full.pair.id[pairi],])))
      colnames(mslr.pair.resj) <- c("pairi","pairj","n.pairi","n.pairj","cv.pairi","cv.pairj","mslr..stat","mslr.p")
      
         for(pairj in 1:length(full.pair.id)){
            
            pair.subset <- full[full$pair.id==full.pair.id[pairi] | full$pair.id==full.pair.id[pairj],]
            pair.subset$pair.id <- as.character(pair.subset$pair.id)
            mslr.pair.resj[pairj,1:2] <- c(as.character(full.pair.id[pairi]),as.character(full.pair.id[pairj]))
            mslr.pair.resj[pairj,3] <- nrow(pair.subset[pair.subset$pair.id == full.pair.id[pairi],])
            mslr.pair.resj[pairj,4] <- nrow(pair.subset[pair.subset$pair.id == full.pair.id[pairj],])
            mslr.pair.resj[pairj,5] <- cv(pair.subset[pair.subset$pair.id == full.pair.id[pairi],]$chase.distance.m) 
            mslr.pair.resj[pairj,6] <- cv(pair.subset[pair.subset$pair.id == full.pair.id[pairj],]$chase.distance.m) 
            # MSLR doesn't work if sample size too small (use 5 samples as cut off)
            if(mslr.pair.resj$n.pairi>4 && mslr.pair.resj$n.pairj[pairj]>4){
               mslr.pair.resj[pairj,7] <- mslr_test(x=pair.subset$chase.distance.m,y=pair.subset$pair.id)$MSLRT
               mslr.pair.resj[pairj,8] <- mslr_test(x=pair.subset$chase.distance.m,y=pair.subset$pair.id)$p_value
            } else {
               mslr.pair.resj[pairj,7] <- NA
               mslr.pair.resj[pairj,8] <- NA 
            }
         }
      
         mslr.pair.res <- rbind(mslr.pair.res,mslr.pair.resj)
      
   }
   
   return(list("mslr.out"=mslr.pair.res, "full.df"=full,"pair.cv.res"=pair.cv.res))
}   
#####  

 

###############################################

# 4.4. RUN CHASE VARIATION FUNCTION 

###############################################

#########
# BEFORE
#########
   
   # ONLY SPECIES WITH 5 SAMPLES OR MORE
   chase.sample5.out <- Chase.Variation.Pairs(cschase.sample5.mean,hschase.sample5.mean)
   boxplot(CV~encounter,data=chase.sample5.out$pair.cv.res,ylab="CV across species pairs",main="before 5 samples")
   boxplot(chase.distance.m~encounter,data=chase.sample5.out$full.df,ylab="chase distance (m)",main="before 5 samples",ylim=c(0,10))
   # boxplot that shows each species pair separately
   boxplot(chase.distance.m~boxplot.label,data=chase.sample5.out$full.df,ylab="chase distance (m)",main="before 5 samples",las=2)
   hist(chase.sample5.out$full.df[chase.sample5.out$full.df$encounter=="heterospecific",4])
   chase.sample5.out$mslr.out[chase.sample5.out$mslr.out$mslr.p<0.051,]
   mslr_test(x=chase.sample5.out$full.df$chase.distance.m,y=chase.sample5.out$full.df$encounter)
   # get coefficient of variation
   cschase.sample5.out.cd <- chase.sample5.out$full.df[chase.sample5.out$full.df$encounter=="conspecific",]
   hschase.sample5.out.cd <- chase.sample5.out$full.df[chase.sample5.out$full.df$encounter=="heterospecific",]
   mean(cschase.sample5.out.cd$chase.distance.m)
   mean(hschase.sample5.out.cd$chase.distance.m)
   cv(cschase.sample5.out.cd$chase.distance.m)
   cv(hschase.sample5.out.cd$chase.distance.m)
   
   # ONLY SPECIES WITH 5 SAMPLES OR MORE AND FOUND IN CS AND HS
   chase.match.out <- Chase.Variation.Pairs(cschase.match.mean,hschase.match.mean)
   boxplot(CV~encounter,data=chase.match.out$pair.cv.res,ylab="CV across species pairs",main="before match")
   boxplot(chase.distance.m~encounter,data=chase.match.out$full.df,ylab="chase distance (m)",main="before match")
   boxplot(CV~pair.id,data=chase.match.out$pair.cv.res[chase.match.out$pair.cv.res$encounter=="heterospecific",],
           ylab="CV within species pairs",main="before match",las=2)
   chase.match.out$mslr.out[chase.match.out$mslr.out$mslr.p<0.051,]
   mslr_test(x=chase.match.out$full.df$chase.distance.m,y=chase.match.out$full.df$encounter)
   # get coefficient of variation
   cschase.match.out.cd <- chase.match.out$full.df[chase.match.out$full.df$encounter=="conspecific",]
   hschase.match.out.cd <- chase.match.out$full.df[chase.match.out$full.df$encounter=="heterospecific",]
   mean(cschase.match.out.cd$chase.distance.m)
   mean(hschase.match.out.cd$chase.distance.m)
   cv(cschase.match.out.cd$chase.distance.m)
   cv(hschase.match.out.cd$chase.distance.m)
####
   
###############################################################
# AFTER
###############################################################
   
   # ONLY SPECIES WITH 5 SAMPLES OR MORE
   chase.post.sample5.out <- Chase.Variation.Pairs(cschase.post.sample5.mean,hschase.post.sample5.mean)
   boxplot(CV~encounter,data=chase.post.sample5.out$pair.cv.res,ylab="CV across species pairs",main="after 5 samples")
   boxplot(chase.distance.m~encounter,data=chase.post.sample5.out$full.df,ylab="chase distance (m)",main="after 5 samples",ylim=c(0,10))
   mslr_test(x=chase.post.sample5.out$full.df$chase.distance.m,y=chase.post.sample5.out$full.df$encounter)
   # get coefficient of variation
   cschase.post.sample5.out.cd <- chase.post.sample5.out$full.df[chase.post.sample5.out$full.df$encounter=="conspecific",]
   hschase.post.sample5.out.cd <- chase.post.sample5.out$full.df[chase.post.sample5.out$full.df$encounter=="heterospecific",]
   mean(cschase.post.sample5.out.cd$chase.distance.m)
   mean(hschase.post.sample5.out.cd$chase.distance.m)
   cv(cschase.post.sample5.out.cd$chase.distance.m)
   cv(hschase.post.sample5.out.cd$chase.distance.m)
   boxplot(cschase.post.sample5.out.cd$chase.distance.m)
   boxplot(hschase.post.sample5.out.cd$chase.distance.m)
   
   # ONLY SPECIES WITH 5 SAMPLES OR MORE AND FOUND IN CS AND HS
   chase.post.match.out <- Chase.Variation.Pairs(cschase.post.match.mean,hschase.post.match.mean)
   boxplot(CV~encounter,data=chase.post.match.out$pair.cv.res,ylab="CV across species pairs",main="after match")
   boxplot(chase.distance.m~encounter,data=chase.post.match.out$full.df,ylab="chase distance (m)",main="after match")
   mslr_test(x=chase.post.match.out$full.df$chase.distance.m,y=chase.post.match.out$full.df$encounter)
   # get coefficient of variation
   cschase.post.match.out.cd <- chase.post.match.out$full.df[chase.post.match.out$full.df$encounter=="conspecific",]
   hschase.post.match.out.cd <- chase.post.match.out$full.df[chase.post.match.out$full.df$encounter=="heterospecific",]
   mean(cschase.post.match.out.cd$chase.distance.m)
   mean(hschase.post.match.out.cd$chase.distance.m)
   cv(cschase.post.match.out.cd$chase.distance.m)
   cv(hschase.post.match.out.cd$chase.distance.m)
   
#####
   
   boxplot(cschase.sample5.out.cd$chase.distance.m,cschase.post.sample5.out.cd$chase.distance.m,
           main="conspecifics",ylab="chase distance")
   boxplot(hschase.sample5.out.cd$chase.distance.m,hschase.post.sample5.out.cd$chase.distance.m,
           main="heterospecifics",ylab="chase distance")
   # how have chase distances changed for heterospecifics?
   plot(density(chase.sample5.out$full.df[chase.sample5.out$full.df$encounter=="heterospecific",4]))
   lines(density(chase.post.sample5.out$full.df[chase.post.sample5.out$full.df$encounter=="heterospecific",4]),lty=2,col=2)
   # significant result in MSLR is for matched species only
   plot(density(chase.match.out$full.df[chase.match.out$full.df$encounter=="heterospecific",4]))
   lines(density(chase.post.match.out$full.df[chase.post.match.out$full.df$encounter=="heterospecific",4]),lty=2,col=2)
   
####################################
# COMPARE BEFORE AND AFTER WITH MSLR
####################################
   
chase.sample5.time <- rbind(chase.sample5.out$full.df,chase.post.sample5.out$full.df)
chase.sample5.time$time <- c(rep("before",nrow(chase.sample5.out$full.df)),rep("after",nrow(chase.post.sample5.out$full.df)))
chase.sample5.time.cs <- chase.sample5.time[chase.sample5.time$encounter=="conspecific",]
chase.sample5.time.hs <- chase.sample5.time[chase.sample5.time$encounter=="heterospecific",]
mslr_test(x=chase.sample5.time$chase.distance.m,y=chase.sample5.time$time)
mslr_test(x=chase.sample5.time.cs$chase.distance.m,y=chase.sample5.time.cs$time)
mslr_test(x=chase.sample5.time.hs$chase.distance.m,y=chase.sample5.time.hs$time)
wilcox_test(chase.sample5.time.hs$chase.distance.m~as.factor(chase.sample5.time.hs$time))


####
# matched data
chase.match.time <- rbind(chase.match.out$full.df,chase.post.match.out$full.df)
chase.match.time$time <- c(rep("before",nrow(chase.match.out$full.df)),rep("after",nrow(chase.post.match.out$full.df)))
chase.match.time.cs <- chase.match.time[chase.match.time$encounter=="conspecific",]
chase.match.time.hs <- chase.match.time[chase.match.time$encounter=="heterospecific",]
mslr_test(x=chase.match.time$chase.distance.m,y=chase.match.time$time)
mslr_test(x=chase.match.time.cs$chase.distance.m,y=chase.match.time.cs$time)
mslr_test(x=chase.match.time.hs$chase.distance.m,y=chase.match.time.hs$time)
wilcox_test(chase.match.time.hs$chase.distance.m~as.factor(chase.match.time.hs$time))
 
chase.sample5.out$pair.cv.res$time <- rep("before",nrow(chase.sample5.out$pair.cv.res))
chase.post.sample5.out$pair.cv.res$time <- rep("after",nrow(chase.post.sample5.out$pair.cv.res))
cv.time <- rbind(chase.sample5.out$pair.cv.res,chase.post.sample5.out$pair.cv.res)
summary(aov(CV~encounter+time+encounter*time,data=cv.time))
boxplot(CV~encounter+time,data=cv.time)
   
############################################################

# 4.5. REGRESSION ANALYSIS

############################################################ 
   
#####################
# PREPARE DATAFRAME WITH PHYLOGENETIC DATA
#####################

# join phylogenetic data to dataframe output of coefficient of variation
chase.match.phylo <- phylo[(phylo$ij %in% as.character(chase.sample5.out$pair.cv.res$pair.id)) |
                                    (phylo$ji %in% as.character(chase.sample5.out$pair.cv.res$pair.id)),]

chase.match.phylo1 <- chase.match.phylo[,c(2,1)]
chase.match.phylo2 <- chase.match.phylo[,c(3,1)]  
colnames(chase.match.phylo1) <- c("pair.id","phylo")
colnames(chase.match.phylo2) <- c("pair.id","phylo")
chase.merge1 <- merge(chase.sample5.out$pair.cv.res,prox.match.phylo1,all.x=T)
chase.merge2 <- merge(chase.merge1,chase.match.phylo2,all.x=T)
chase.pre.merged <- chase.merge2
chase.phylo <- aggregate(CV~pair.id+phylo,data=chase.pre.merged,FUN=mean)

# add body size difference
chase.phylo.size <- merge(chase.phylo,hschase.sample5[,c(23,26)],by="pair.id")
chase.phylo.size <- aggregate(CV~pair.id+size.diff+phylo,data=chase.phylo.size,mean)
#####


#################################################
# REGRESSION WITH PHYLOGENETIC RELATEDNESS AS PREDICTOR
################################################

plot(chase.phylo.size$size.diff,chase.phylo.size$phylo,xlab="size difference",
     ylab="phylogenetic relatedness",cex.lab=1.5)
lines(lowess(chase.phylo$phylo,chase.phylo$size.diff,f=2/3))
cor.test(chase.phylo.size$size.diff,chase.phylo.size$phylo,method="spearman")

plot(chase.phylo.size$size.diff,chase.phylo.size$CV,xlab="size difference",
     ylab="chase variation",cex.lab=1.5)
lines(lowess(chase.phylo.size$size.diff,chase.phylo.size$CV,f=2/3))

# phylogenetic relatedness
plot(chase.phylo$phylo,chase.phylo$CV,xlab="phylogenetic distance",
     ylab="coefficient of variation (chase distance)",cex.lab=1.5)
lines(lowess(chase.phylo$phylo,chase.phylo$CV,f=2/3))

lm.chase.phylo <- lm(CV~phylo+size.diff,data=chase.phylo.size)
# add in quadratic
chase.phylo$phylo2 <- chase.phylo$phylo^2
chase.phylo2 <- chase.phylo[order(chase.phylo$phylo2),]
lm.chase.phylo.quad <- lm(CV~phylo+phylo2,data=chase.phylo2)
# plot QQ plot etc to explore model fit
par(mfcol=c(4,2))
plot(lm.chase.phylo)
plot(lm.chase.phylo.quad)
# What is R-squared value for each?
summary(lm.chase.phylo)
summary(lm.chase.phylo.quad)
# Use AICc to see which model is best fit
AICc(lm.chase.phylo) # best one
AICc(lm.chase.phylo.quad) 

###########



######################################
#####################################

# 5. SAVE DATA FOR FIGURES

#####################################
#####################################

###################
## Run proximity function for data reduced to signalling/chase only
###################

# SIGNALLING
proxS.sample5.out <- Proximity.Diff(cssig.sample5,hssig.sample5,cssig.post.sample5,hssig.post.sample5,"Minimum 5 samples")
# plot stacked bar for proportions of encounters with each proximity
proxS.sample5.out$round.prox <- round(proxS.sample5.out$mean.proximity,0)
proxS.plot <- as.data.frame(table(proxS.sample5.out$round.prox,proxS.sample5.out$encounter,proxS.sample5.out$time))
colnames(proxS.plot) <- c("proximity","encounter","time","freq")
proxS.plot$freq <- proxS.plot$freq
proxS.plot$label <- paste(proxS.plot$time,proxS.plot$encounter,sep=" ")
proxS.plot$label <- ordered(proxS.plot$label, levels=c("before conspecific","before heterospecific","after conspecific","after heterospecific"))
# use log proportion of encounters so it is visually clear
# NB. this still fits the Mann Whitney U stats because that is done on rank values
ggplot(proxS.plot, aes(x = label, y = log(freq), fill = proximity)) +
  geom_col(colour = "black", position = "fill") +
  theme(axis.text.x = element_text(size = 16, angle = 45,hjust=1),
        axis.text.y = element_text(size = 16, angle = 0),
        axis.title.x=element_blank(),axis.title.y=element_text(size=18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
        scale_fill_grey(name = "proximity of aggression initiation (cm)", labels = c("0-24","25-49","50-74","75-100"),start=1,end=0.25) +
        labs(y="log(proportion of encounters)")

# CHASE
proxC.sample5.out <- Proximity.Diff(cschase.sample5,hschase.sample5,cschase.post.sample5,hschase.post.sample5,"Minimum 5 samples")
# plot stacked bar for proportions of encounters with each proximity
proxC.sample5.out$round.prox <- round(proxC.sample5.out$mean.proximity,0)
proxC.plot <- as.data.frame(table(proxC.sample5.out$round.prox,proxC.sample5.out$encounter,proxC.sample5.out$time))
colnames(proxC.plot) <- c("proximity","encounter","time","freq")
proxC.plot$freq <- proxC.plot$freq
proxC.plot$label <- paste(proxC.plot$time,proxC.plot$encounter,sep=" ")
proxC.plot$label <- ordered(proxC.plot$label, levels=c("before conspecific","before heterospecific","after conspecific","after heterospecific"))
# use log proportion of encounters so it is visually clear
# NB. this still fits the Mann Whitney U stats because that is done on rank values
ggplot(proxC.plot, aes(x = label, y = log(freq), fill = proximity)) +
  geom_col(colour = "black", position = "fill") +
  theme(axis.text.x = element_text(size = 16, angle = 45,hjust=1),
        axis.text.y = element_text(size = 16, angle = 0),
        axis.title.x=element_blank(),axis.title.y=element_text(size=18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
        scale_fill_grey(name = "proximity of chase initiation (cm)", labels = c("0-24","25-49","50-74","75-100"),start=1,end=0.25) +
        labs(y="log(proportion of encounters)")
###################

###################
# FIGURE S1, S2
save(proximity.sample.sizes,file="FigS1_SampleSize_Proximity.RData")
save(signalling.sample.sizes,file="FigS1_SampleSize_Signal.RData")
save(chase.sample.sizes,file="FigS1_SampleSize_Chase.RData")

# FIGURE 1
save(prox.plot,prox.match.plot,proxS.plot,file="Fig1_ProximityCategories.RData")

# FIGURE 2
save(sig.before5.bs.mean,sig.before.match.bs.mean,sig.after5.bs.mean,
     sig.after.match.bs.mean,file="Fig2_SignalProportion.RData")

# FIGURE 3
save(chase.match.out,chase.post.match.out,file="Fig3_ChaseDistanceMatch.RData")

# FIGURE 4
save(lm.prox.phylo,lm.sig.phylo,lm.chase.phylo,
     prox.phylo2,sig.phylo,chase.phylo2,file="Fig4_PartialCoefficients_June2022.RData")
#################




