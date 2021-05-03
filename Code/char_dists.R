#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
### Imports ###
require(tidyverse)

### Functions ###
filt_pairs <- function(data){
    require(sets)
    k <- nrow(data)
    
    # Converts pairs to sets so that the order does not matter, because it doesn't
    # Then removes redundant pairs
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

    # List meaningful pairs for comparison
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
    x <- unique(a[,1])
    y <- unique(b[,1])

    data <- expand.grid(x, y, 0)
    print(data)
    print(a)
    print(b)
    for (i in 1:nrow(data)) {
        A <- match(data[i, 1], a[,1])
        B <- match(data[i, 2], b[,1])
        print(c(A,B))
        if (a[A,4] && b[B, 4]) {
            data[i, 3] <- 10000
        }
        else if (a[A,5] && b[B, 5]) {
            data[i, 3] <- 10000
        }
        else if (a[A,3] && b[B, 3]) {
            data[i, 3] <- 10
        }
        else if (a[A,2] == b[B, 2]) {
            data[i, 3] <- 1
        }
    }
    print(names)

    dist_string <- paste("../Data/", names[1], names[2], ".txt", sep = "")
    result_string <- paste("../Results/", names[1], names[2], sep = "")
    write.table(data, file = dist_string, sep = "\t", row.names = F, col.names = F)
    return(c(dist_string, result_string))
}

file_order <- function(data1, data2, names){
    #reorder pair, smallest first
    Q <- names[1]
    T <- names[2]
    dif <- nrow(data2) - nrow(data1)
    if (dif < 0){
        Q <- names[2]
        T <- names[1]
    }
    Q_str <- paste("../Data/", list.files(path = "../Data/", pattern = paste(Q, ".gw", sep = "")), sep = "")
    T_str <- paste("../Data/", list.files(path = "../Data/", pattern = paste(T, ".gw", sep = "")), sep = "")
    Q_T <- c(Q_str,T_str)
    print(Q_T)
    return(Q_T)
}

Generate <- function(data, arg){
    pairs <- data[[1]]
    l <- nrow(pairs)
    for (i in 1:l){
        names <- make.names(gsub(x = data[[1]][i,], pattern = "*.csv$", ""))
        data1 <- data[[2]][[names[1]]]
        data2 <- data[[2]][[names[2]]]
        
        dist <- dist_matrix(data1, data2, names)
        Graphs <- file_order(data1, data2, names)
        
        # Generates appropriate bash file for running MI-GRAAL
        line <- paste("./MI-GRAALRunner.py", Graphs[1], Graphs[2], dist[2], "-p", arg, "-q", dist[1], sep = " ")
        print(line)
        write(line, file="GRAAL_it.sh", append=TRUE)
    }
}

## Main ##
Data <- Collect()
Generate(Data, args[1])
 