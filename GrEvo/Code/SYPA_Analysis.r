#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
### Imports ###

### Functions ###

## Main ##
GREV = read.csv("../Results/by70_sim_by4_steep_mat.csv")
SYPA = read.csv("SYPA_S3_25_Scores.csv")
library(reshape2)
GREV = melt(GREV,varnames=c('X1', 'X2'), value.name = "GREVO")
SYPA = melt(SYPA,varnames=c('X1', 'X2'), value.name = "SYPA")
MAGNA = read.csv("Dist100Mag.csv", header=TRUE)
head(MAGNA)
head(SYPA)
head(GREV)

SYGR = merge(SYPA,GREV, by = c('X', 'variable'))
SYGR = merge(SYGR,MAGNA, by = c('X', 'variable'))
colnames(SYGR)[1:2] <- c("X1", "X2")
head(SYGR)
library(tidyverse)
library(ggplot2)

SYGR <- SYGR %>% add_row(X1 = "Dicksonosteus", X2 = "Dicksonosteus", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Lunaspis", X2 = "Lunaspis", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Entelognathus", X2 = "Entelognathus", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Eusthenopteron", X2 = "Eusthenopteron", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Cowralepis", X2 = "Cowralepis", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Eurycaraspis", X2 = "Eurycaraspis", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Diandongpetalichthys", X2 = "Diandongpetalichthys", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Meemannia", X2 = "Meemannia", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Dialipina", X2 = "Dialipina", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Wuttagoonaspis", X2 = "Wuttagoonaspis", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Actinolepis", X2 = "Actinolepis", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Romundina", X2 = "Romundina", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Bothriolepis", X2 = "Bothriolepis", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Buchanosteus", X2 = "Buchanosteus", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Brindabellaspis", X2 = "Brindabellaspis", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Moythomasia", X2 = "Moythomasia", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Guiyu", X2 = "Guiyu", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) %>%
    add_row(X1 = "Cheirolepis", X2 = "Cheirolepis", GREVO = 0, SYPA = 1, MAGNA = 1, .before = 1) 
z = 1.96
SYPAA <- SYGR %>%
    group_by(GREVO) %>%
    summarise(SYPA_Mean = mean(SYPA), MAGNA_Mean = mean(MAGNA), 
    SYPA_L = mean(SYPA) - sd(SYPA)*z, 
    MAGNA_L = mean(MAGNA) - sd(MAGNA)*z, 
    SYPA_U = mean(SYPA) + sd(SYPA)*z, 
    MAGNA_U = mean(MAGNA) + sd(MAGNA)*z)
print(SYPAA)

library(scales)
#plot
col1 <- "darkgreen"
col2 <- "blue"

png("SYGGGR2.png")
plot(SYGR$GREVO, SYGR$SYPA, col = col1, pch = 20, 
xlab="GrEvo Estimated Parsimony Score", ylab="Alignment Score", ylim=c(0,1),
xlim=c(0, 18))
points(SYGR$GREVO, SYGR$MAGNA, col=col2, pch = 20)
lines(SYPAA$GREVO, SYPAA$SYPA_Mean, type = "l", col =col1)
lines(SYPAA$GREVO, SYPAA$MAGNA_Mean, type = "l", col=col2)
SYPAA = SYPAA[complete.cases(SYPAA),]
head(SYPAA)
polygon(c(SYPAA$GREVO, rev(SYPAA$GREVO)),c(SYPAA$SYPA_L, rev(SYPAA$SYPA_U)), col = alpha("green", 0.3), border=F)
#lines(SYPAA$GREVO, SYPAA$SYPA_U, col=rgb(1, 0, 0,0.5),lty=2)
#lines(SYPAA$GREVO, SYPAA$SYPA_L, col=rgb(1, 0, 0,0.5),lty=2)
polygon(c(SYPAA$GREVO, rev(SYPAA$GREVO)),c(SYPAA$MAGNA_L, rev(SYPAA$MAGNA_U)), col = alpha(col2, 0.2), border=F)
#lines(SYPAA$GREVO, SYPAA$MAGNA_U, col=rgb(0, 0, 1, 0.5),lty=2)
#lines(SYPAA$GREVO, SYPAA$MAGNA_L, col=rgb(0, 0, 1, 0.5),lty=2)

legend("topright",
c("SYPA","MAGNA++"),
fill=c("green",col2)
)
dev.off()
