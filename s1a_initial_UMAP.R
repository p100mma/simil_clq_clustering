#library(matrixStats)
library(magrittr)
library(umap)
load("trainLeukemia.RData")

train$class-> ref
X<- train[,-c(1,2, ncol(train))] %>% as.matrix()
apply(X,2,var)-> v_x
variance_cutoffs=c(100,500,1000,
		   2000,3500, 5000,
		 7500, 10000, 12500,
		 15000, 15750, 20000 ) +1   # 12

args=commandArgs(trailingOnly=TRUE)	
as.integer(args[[1]])-> n_vc
vc<- variance_cutoffs[[n_vc]]
X<- X[, v_x> sort(v_x, decreasing = TRUE)[vc] ]

set.seed(n_vc)
sample.int(23523)-> umap_seed
set.seed(umap_seed)
t0<- Sys.time()
UM<- umap(X, n_components=3, n_neighbors=45)$layout # ~ sqrt(nrow(X))
execution_time= Sys.time() - t0
attr(UM, "execution_time")<- execution_time
saveRDS(UM,sprintf("leuk_UM/n_genes=%d_UMkD.rds",n_vc))


