#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
### Imports ###
require(tidyverse)
library(igraph)

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

dist_matrix <- function(a, b, names, path){
    y <- a
    x <- b
    
    # Measure difference in number of nodes
    dif = nrow(y) - nrow(x)

    #reorder so smaller graph is rows
    if (dif < 0) {
        x <- a
        y <- b
    }
    print(c(nrow(x), nrow(y)))
    # Generate empty matrix x by y
    data <- matrix(, nrow = nrow(x), ncol = nrow(y))

    # Fill matrix with similarity scores
    for (i in 1:nrow(x)) {
        for (j in 1:nrow(y)) {

            data[i, j] <- 0
            

        }
    }
    
    dist_string <- paste("../Results/", names[1], names[2], ".txt", sep = "")
    result_string <- paste("../Results/", names[1], names[2], sep = "")
    
    title = paste(nrow(x), nrow(y), sep = " ")
    #Check its existence
    if (file.exists(dist_string)) {
        #Delete file if it exists
        file.remove(dist_string)
    }
    cat(title, file=dist_string)
    cat("\n", file=dist_string, append = TRUE)
    write.table(data, file = dist_string, sep = " ", 
        row.names = F, col.names = F, append = TRUE)
    return(c(dist_string, result_string))
}

file_order <- function(data1, data2, names, path){
    #reorder pair, smallest first
    Q <- names[1]
    T <- names[2]
    
    dif <- nrow(data2) - nrow(data1)
    if (dif < 0){
        Q <- names[2]
        T <- names[1]
    }
    print(names)
    Q_str <- paste(path, Q, ".txt", sep = "")
    T_str <- paste(path, T, ".txt", sep = "")
    Q_T <- c(Q_str,T_str)
    print(Q_T)
    return(Q_T)
}

Generate <- function(data, path = "../Data/", alpha = 0.8, population = 100, generations = 100, threads = 10){
    #Check README.txt of MAGNA for details on the parameters
    
    pairs <- data[[1]]
    l <- nrow(pairs)
    shcript = "MAGNA_it.sh"
    #Check its existence
    if (file.exists(shcript)) {
        #Delete file if it exists
        file.remove(shcript)
    }
    for (i in 1:l){
        names <- make.names(gsub(x = data[[1]][i,], pattern = "*.csv$", ""))
        data1 <- data[[2]][[names[1]]]
        data2 <- data[[2]][[names[2]]]
        
        dist <- dist_matrix(data1, data2, names, path)
        Graphs <- file_order(data1, data2, names, path)
        
        # Generates appropriate bash file for running MI-GRAAL
        line <- paste("./magnapp_cli_linux64 -G", Graphs[1], "-H", 
            Graphs[2], "-o", dist[2], "-d", dist[1], "-m S3", "-p", population, "-t", threads, "-a", alpha, "-n", generations, "-f 1",  sep = " ")
        print(line)
        write(line, file="MAGNA_it.sh", append=TRUE)
    }
}
## Main ##
Data <- Collect(args[1])
Generate(Data, args[1], args[2], args[3], args[4], args[5])
