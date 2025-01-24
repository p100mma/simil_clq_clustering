X<- readRDS("gene_expr_data.rds")
library(matrixStats)
source("cor_matmul.R")
A<-abs(corfast(X))^5
A[A>1]=1; A[A<0]=0; #roundoff errors, adjust so WGCNA does not complain
library(WGCNA)
TOMdist(A)-> dTOM; rm(A)
library(fastcluster)
hclust(as.dist(dTOM), "average")-> dendro
saveRDS(dendro,"BRCA_hclust_avg_dendro.rds")
library(dynamicTreeCut)
min_cl_sizes<- c(30, 400, 860)
labels_per_ms<-list()
for (i in seq_along(min_cl_sizes))
labels_per_ms[[i]]<-cutreeDynamic(dendro=dendro, 
	      method="hybrid", distM= dTOM,
	      deepSplit=2,minClusterSize=min_cl_sizes[[i]])

library(cliqueClusteR)
labels_per_ms<-lapply(labels_per_ms, uniqueSingletonLabels)

saveRDS(labels_per_ms,"WGCNA_labels_fulldata.rds")

