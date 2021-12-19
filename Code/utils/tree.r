#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

### Imports ###
library(ape)
library(phytools)

### Functions ###

## Main ##
names<- c("Dicksonosteus","Lunaspis","Entelognathus","Eusthenopteron","Cowralepis",
"Eurycaraspis","Diandongpetalichthys","Meemannia","Dialipina","Wuttagoonaspis",
"Actinolepis","Romundina","Bothriolepis","Buchanosteus","Brindabellaspis",
"Moythomasia","Guiyu","Cheirolepis")

tree = ape::read.nexus(args[1])
tree1=tree$Strict
tree2=tree$Adams
ii<-sapply(names,grep,tree$tip.label)

png(args[2],
  width     = 8,
  height    = 11,
  units     = "in",
  res       = 600,
  pointsize = 9)
par(mfrow=c(1,2))
plotTree(tree1,cex=200)
add.arrow(tree1,tip=names,arrl=0.1, hedl=0.01)
plotTree(tree2,cex=200)
add.arrow(tree2,tip=names,arrl=0.1, hedl=0.01)
dev.off()