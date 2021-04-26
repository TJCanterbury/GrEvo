#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
file1 <- paste(args[1])
file2 <- paste(args[2])
a <- read.csv(file1, header = T)
b <- read.csv(file2, header = T)

x <- 1:nrow(a)
y <- 1:nrow(b)

data <- expand.grid(x, y, 0)
for (i in 1:nrow(data)) {
    A <- data[i, 1]
    B <- data[i, 2]
    if (a[A,4] && b[B, 4]) {
        data[i, 3] <- 10000000
    }
    else if (a[A,5] && b[B, 5]) {
        data[i, 3] <- 10000000
    }
    else if (a[A,3] && b[B, 3]) {
        data[i, 3] <- 10000
    }
    else if (a[A,2] == b[B, 2]) {
        data[i, 3] <- 10
    }
}

result_string <- paste("../Data/", args[3], ".txt", sep = "")

write.table(data, file = result_string, sep = "\t", row.names = F, col.names = F)