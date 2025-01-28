library(FCPS)
library(dbscan)
sets2test<- c("Atom", "Chainlink", "Tetra", "Target", "TwoDiamonds")
set="TwoDiamonds"
db_eps<- c(4,0.11,0.5,0.4,0.11) # see docs/DB_scan_parameters.Rmd for derivation of eps
names(db_eps)<-sets2test
db_scores<-vector()
for (set in sets2test)
{
  db_input<- get(set)$Data
  ref<-get(set)$Cls
 knn_curve<-sort(kNNdist(db_input, k = 2*ncol(db_input) - 1))
 pdf(sprintf("%s_kNN_curve.pdf",set))
 plot(knn_curve,main=sprintf("%s kNN curve",set))
 dev.off()
 abline(h=db_eps[[set]])
 print(set)
 igraph::compare(
   dbscan(db_input, eps=db_eps[[set]], minPts = 2*ncol(db_input))$cluster,
   ref,"adjusted.rand")-> ARI
 db_scores[[length(db_scores)+1]]=ARI
}
names(db_scores)<-sets2test
saveRDS(db_scores,"DBscan_FCPS_scores.rds")
