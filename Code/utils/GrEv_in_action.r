#!/usr/bin/env Rscript

parsdata=read.csv("../Results/in_action.csv")
print(parsdata)
png("in_act.png",
  width     = 5,
  height    = 5,
  units     = "in",
  res       = 600,
  pointsize = 9)
plot(parsdata, type="l", xlab="Step", ylab="1 - CSS Score")
dev.off()