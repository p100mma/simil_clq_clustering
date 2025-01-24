
library(matrixStats)
input_path= 'gene_expr_data.rds'
message("reading input")
X<- readRDS(input_path)
source('cor_matmul.R')
message("computing correlation")
S<- similarity_matrix( fastPearsonData(X) )
message("saving graph")
saveRDS(S,'BRCA_similarity_graph.rds')
message("cleaning")
rm(S)
message("done")
