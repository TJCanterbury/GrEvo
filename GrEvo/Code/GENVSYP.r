#!/usr/bin/env Rscript
data <- read.csv("GenvS3_SYPA.txt", sep=" ")
library(ggplot2)
theme_set(
  theme_bw())
colnames(data) = c("id", "variable", "value")
png("GENSYP.png")
ggplot(data, aes(variable, value))+
    geom_point(alpha=0.3)+
    stat_smooth(method="loess", se=TRUE, colour="green")+
    xlab("Number of Repeats")+
    ylab("CSS Score") +
    theme(text=element_text(size=20))

dev.off()
