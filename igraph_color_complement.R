igraph_color_complement<- function(S){
A<- S>0
A<- !A
G_cA<- igraph::graph_from_adjacency_matrix(A,mode="undirected", weighted=FALSE, diag=FALSE)
rm(A); rm(S);
list(membership=igraph::greedy_vertex_coloring(G_cA))
}
