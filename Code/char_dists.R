#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
### Imports ###
require(tidyverse)

### Functions ###
filt_pairs <- function(data){
    require(sets)
    k <- nrow(data)
    for (i in 1:k) {
        for(j in 1:k){
            if(as.set(data[i,]) == as.set(data[j,]) && i != j){
                data[j,] <- c(0,0)
            }
        }
    }
    detach("package:sets")
    row_sub = apply(data, 1, function(row) all(row !="0" ))
    return (data[row_sub,])
}

Collect <- function(path = "../Data/"){
    # Enter data directory
    ret_dir <- getwd()
    setwd(path)

    # Gather files and file names
    temp <- list.files(pattern = "*.csv")
    myfiles = lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), read.csv)

    # list pairs
    pairs <- expand.grid(temp, temp) %>%
        filter(Var1 != Var2)
    pairs <- matrix(as.matrix(pairs), ncol = 2, dimnames = NULL)
    pairs <- filt_pairs(pairs)
    
    # return from data directory
    setwd(ret_dir)

    # return pairs
    return(list(pairs, myfiles))
}

dist_matrix <- function(a, b, names){
    x <- 1:nrow(a)
    y <- 1:nrow(b)

    data <- expand.grid(x, y, 0)
    for (i in 1:nrow(data)) {
        A <- data[i, 1]
        B <- data[i, 2]
        if (a[A,4] && b[B, 4]) {
            data[i, 3] <- 1
        }
        else if (a[A,5] && b[B, 5]) {
            data[i, 3] <- 1
        }
        else if (a[A,3] && b[B, 3]) {
            data[i, 3] <- 0.5
        }
        else if (a[A,2] == b[B, 2]) {
            data[i, 3] <- 0.25
        }
    }

    result_string <- paste("../Data/", names[1], names[2], ".txt", sep = "")
    print(result_string)
    write.table(data, file = result_string, sep = "\t", row.names = F, col.names = F)
}

Generate <- function(data){
    pairs <- data[[1]]
    l <- nrow(pairs)
    for (i in 1:l){
        names <- make.names(gsub(x = data[[1]][i,], pattern = "*.csv$", ""))
        dist_matrix(data[[2]][[names[1]]], data[[2]][[names[2]]], names)
        print(names)
    }
}

## Main ##
Data <- Collect()
Generate(Data)
