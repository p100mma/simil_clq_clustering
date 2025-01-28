#library(matrixStats)
library(magrittr)
library(dbscan)
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

UM<-readRDS(sprintf("leuk_UM/n_genes=%d_UMkD.rds",n_vc))

sort(kNNdist(UM, k=2*ncol(UM)-1))-> kNN_d_curve
saveRDS(kNN_d_curve, sprintf("leuk_DB_scan_kNN_curve/n_genes=%d_kNN_d_curve.rds",n_vc))

pdf(sprintf("leuk_DB_scan_kNN_curve/curve_n_genes=%d.pdf",n_vc))
plot(kNN_d_curve, main=sprintf("n_genes variant=%d_kNNcurve",n_vc))
dev.off()

#t0<- Sys.time()
#execution_time= Sys.time() - t0
