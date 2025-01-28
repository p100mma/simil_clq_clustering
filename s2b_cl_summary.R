library(igraph)
library(magrittr)
library(cliqueClusteR)
load("trainLeukemia.RData")
train$class-> ref

strategies=c("none","complexity","mst")

stats_per_vc<- list()

for (n_vc in 1:12){
   ARI_per_strat<-list()
   time_per_strat<-list()
   thr_per_strat<-list()
   for (strat in strategies){
   ARI_per_strat[[strat]]<-list()
   time_per_strat[[strat]]<-list()
   thr_per_strat[[strat]]<-vector()
   }

   for (trial in 1:100){
    results_per_strategy<-readRDS(
	   sprintf("leuk_clusters/n_genes=%d_trial=%d_clusters.rds", n_vc,
		   trial)
	   )
	   for (strat in strategies) {
		   lapply( results_per_strategy[[strat]][1:(-1+length(results_per_strategy[[strat]]))],
			   function(x) {
			    if(any(is.na(x))) return(NA)
			    if (class(x)!="character") return( compare(ref, 
							uniqueSingletonLabels(x),
						 	"adjusted.rand")
					      		      )
			    return( compare(ref,x,"adjusted.rand") )
					}	
		           ) %>% unlist() -> ARI_per_strat[[strat]][[trial]]	
		   lapply( results_per_strategy[[strat]][1:(-1+length(results_per_strategy[[strat]]))],
			   function(x) if (any(is.na(x))) NA else attr(x, "exec_time")) %>%
		   unlist() -> time_per_strat[[strat]][[trial]]	
		   thr_per_strat[[strat]][[trial]]<- results_per_strategy[[strat]]$thr
			
			}
	   }
	for (strat in strategies) {
	     ARI_per_strat[[strat]]<-do.call(rbind, ARI_per_strat[[strat]])
	     time_per_strat[[strat]]<-do.call(rbind, time_per_strat[[strat]])
	}
	stats_per_vc[[n_vc]]= list(ARI=ARI_per_strat, time=time_per_strat, thr=thr_per_strat)	
	}

saveRDS(stats_per_vc,"leuk_clusters_results_aggregated.rds")
