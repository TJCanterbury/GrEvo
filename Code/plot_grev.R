library(igraph)

G = read_graph("../pond/b.txt", directed=F)

Isolated = which(degree(G)==0)
G2 = delete.vertices(G, Isolated)

plot(G2)