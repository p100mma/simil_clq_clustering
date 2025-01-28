#library(matrixStats)
library(magrittr)
library(dbscan)
library(apcluster)
library(igraph)
library(umap)
library(cliqueClusteR)
load("trainLeukemia.RData")

train$class-> ref
X<- train[,-c(1,2, ncol(train))] %>% as.matrix()
apply(X,2,var)-> v_x
variance_cutoffs=c(100,500,1000,
		   2000,3500, 5000,
		 7500, 10000, 12500,
		 15000, 17500, 20000 ) +1   # 12

args=commandArgs(trailingOnly=TRUE)	
as.integer(args[[1]])-> n_vc

# for separate UMAP trials
set.seed(n_vc)
n_trials=500
seeds<- sample.int(3223432, n_trials)


vc<- variance_cutoffs[[n_vc]]
X<- X[, v_x> sort(v_x, decreasing = TRUE)[vc] ]

as.integer(args[[2]])-> trial_number

message( sprintf("preparing UMAP, trial = %d, n_genes= %d",trial_number,
		 n_vc))
set.seed(seeds[[trial_number]])
umap(X, n_components=3, n_neighbors=45)$layout -> UM


#run clique based for 2 variants to get unthr. similarity S
#and threshold values
message("initial clique based clusters:")

t0= Sys.time()
run_scpcss(UM, do_signif = FALSE, X_simil_fun = "euclidean", t_S = "complexity", t_S_init_steps = 30,
           t_CS_relax_method = "components", return_CS = TRUE, return_S=TRUE,
           t_CS_n_init_steps = 30,
             t_CS = "modularity")-> clq_based_comp_mod
attr(clq_based_comp_mod, "exec_time")= Sys.time()-t0


t0= Sys.time()
run_scpcss(UM, do_signif = FALSE, X_simil_fun = "euclidean", t_S = "mst", t_S_init_steps = 30,
           t_CS_relax_method = "components", return_CS = TRUE, return_S=FALSE,
           t_CS_n_init_steps = 30,
             t_CS = "modularity")-> clq_based_mst_mod
attr(clq_based_mst_mod, "exec_time")= Sys.time()-t0

t0= Sys.time()
run_scpcss(UM, do_signif = TRUE, X_simil_fun = "euclidean", t_S = "complexity", t_S_init_steps = 30,
           t_CS_relax_method = "components", return_CS = TRUE, return_S=FALSE,
           t_CS_n_init_steps = 30,
             t_CS = "modularity")-> clq_based_comp_mod_sT
attr(clq_based_comp_mod, "exec_time")= Sys.time()-t0


t0= Sys.time()
run_scpcss(UM, do_signif = TRUE, X_simil_fun = "euclidean", t_S = "mst", t_S_init_steps = 30,
           t_CS_relax_method = "components", return_CS = TRUE, return_S=FALSE,
           t_CS_n_init_steps = 30,
             t_CS = "modularity")-> clq_based_mst_mod_sT
attr(clq_based_mst_mod, "exec_time")= Sys.time()-t0

S<- clq_based_comp_mod$S 

#will contain other thresholds after first pass of tests
thresholds_S<-list(none=0, 
			      complexity=clq_based_comp_mod$t_S,
			      mst=clq_based_mst_mod$t_S	
			      )
thr_strats<- names(thresholds_S)

#initial unthresholded runs:
	#info for k-means and hierarchical clustering
n_ref<- length(unique(ref))

#DB-scan radius and k:
db_eps=0.3 #see elbow point : leuk_DB_scan_kNN_curve/joint_plot.pdf
db_n_min=2*ncol(UM)

results_per_strategy<- list()
message("working on other algorithms..")

for (strategy in thr_strats) {
	message(strategy)
	thr<- thresholds_S[[strategy]]
	this_thr_results<- rep(NA, 14) %>% as.list()
	names(this_thr_results)<- c("clq_modularity","clq_MSTDunn",
				    "clq_modularity_sT", "clq_MSTDunn_sT",
				 "infomap","label_prop","leiden",
				 "affinity_prop","hcl_single",
				 "hcl_complete","hcl_average",
				 "hcl_ward.D2","db_scan",
   				 "k_means")
	###add clq based after the loop because it cannot be used without
	###thresholding

	t0=Sys.time()
	this_thr_results[["infomap"]]=igraph::cluster_infomap(S %thr% thr %>% igr_())$membership
	attr(this_thr_results[["infomap"]],"exec_time")= Sys.time() - t0


	t0=Sys.time()
	this_thr_results[["label_prop"]]=igraph::cluster_label_prop(S %thr% thr %>% igr_())$membership
	attr(this_thr_results[["label_prop"]],"exec_time")= Sys.time() - t0


	t0=Sys.time()
	this_thr_results[["leiden"]]=igraph::cluster_leiden(S %thr% thr %>% igr_())$membership
	attr(this_thr_results[["leiden"]],"exec_time")= Sys.time() - t0
	
	
	t0=Sys.time()
	this_thr_results[["affinity_prop"]]=apcluster(S %thr% thr) %>% labels()
	attr(this_thr_results[["affinity_prop"]],"exec_time")= Sys.time() - t0
	  

	t0=Sys.time()	
  	this_thr_results[["hcl_single"]]=hclust(as.dist( S %thr% thr %>% flip_() ),method = "single") %>% cutree(k =n_ref) 
	attr(this_thr_results[["hcl_single"]],"exec_time")= Sys.time() - t0


	t0=Sys.time()	
  	this_thr_results[["hcl_complete"]]=hclust(as.dist( S %thr% thr %>% flip_() ),method = "complete") %>% cutree(k =n_ref) 
	attr(this_thr_results[["hcl_complete"]],"exec_time")= Sys.time() - t0


	t0=Sys.time()	
  	this_thr_results[["hcl_average"]]=hclust(as.dist( S %thr% thr %>% flip_() ),method = "average") %>% cutree(k =n_ref) 
	attr(this_thr_results[["hcl_average"]],"exec_time")= Sys.time() - t0


	t0=Sys.time()	
  	this_thr_results[["hcl_ward.D2"]]=hclust(as.dist( S %thr% thr %>% flip_() ),method = "ward.D2") %>% cutree(k =n_ref) 
	attr(this_thr_results[["hcl_ward.D2"]],"exec_time")= Sys.time() - t0

	if (strategy=="none") {
	#db scan and k-means only work here
	t0=Sys.time()
	this_thr_results[["db_scan"]]=dbscan(UM, eps= db_eps,
						 minPts= db_n_min)$cluster
	attr(this_thr_results[["db_scan"]],"exec_time")= Sys.time()-t0
	

	t0=Sys.time()
	this_thr_results[["k_means"]]=kmeans(UM, centers=n_ref
						)$cluster
	attr(this_thr_results[["k_means"]],"exec_time")= Sys.time()-t0


	}
	
	this_thr_results$thr = thr	
	results_per_strategy[[ length(results_per_strategy) +1 ]] = this_thr_results
	} #algo testing for END
	
names(results_per_strategy)<- thr_strats

#add clq based variants

message("adding 2 missing clique based variants")

results_per_strategy[["complexity"]][["clq_modularity"]]= clq_based_comp_mod$cluster_membership$node	
results_per_strategy[["complexity"]][["clq_modularity_sT"]]= clq_based_comp_mod_sT$cluster_membership$node	
results_per_strategy[["mst"]][["clq_modularity"]]= clq_based_mst_mod$cluster_membership$node
results_per_strategy[["mst"]][["clq_modularity_sT"]]= clq_based_mst_mod$cluster_membership$node
results_per_strategy[["mst"]][["clq_modularity_sT"]]= clq_based_mst_mod_sT$cluster_membership$node

t0=Sys.time()
results_per_strategy[["complexity"]][["clq_MSTDunn"]]= run_scpcss( UM,
do_signif = FALSE, X_simil_fun = "euclidean", t_S = "complexity", t_S_init_steps = 30,
           t_CS_relax_method = "components", return_CS = TRUE, return_S=FALSE,
           t_CS_n_init_steps = 30,
             t_CS = "mst")$cluster_membership$node
attr(results_per_strategy[["complexity"]][["clq_MSTDunn"]], "exec_time")= Sys.time()-t0

t0=Sys.time()
results_per_strategy[["complexity"]][["clq_MSTDunn_sT"]]= run_scpcss( UM,
do_signif = TRUE, X_simil_fun = "euclidean", t_S = "complexity", t_S_init_steps = 30,
           t_CS_relax_method = "components", return_CS = TRUE, return_S=FALSE,
           t_CS_n_init_steps = 30,
             t_CS = "mst")$cluster_membership$node
attr(results_per_strategy[["complexity"]][["clq_MSTDunn_sT"]], "exec_time")= Sys.time()-t0

t0=Sys.time()
results_per_strategy[["mst"]][["clq_MSTDunn"]]= run_scpcss( UM,
do_signif = FALSE, X_simil_fun = "euclidean", t_S = "mst", t_S_init_steps = 30,
           t_CS_relax_method = "components", return_CS = TRUE, return_S=FALSE,
           t_CS_n_init_steps = 30,
             t_CS = "mst")$cluster_membership$node
attr(results_per_strategy[["mst"]][["clq_MSTDunn"]], "exec_time")= Sys.time()-t0

t0=Sys.time()
results_per_strategy[["mst"]][["clq_MSTDunn_sT"]]= run_scpcss( UM,
do_signif = TRUE, X_simil_fun = "euclidean", t_S = "mst", t_S_init_steps = 30,
           t_CS_relax_method = "components", return_CS = TRUE, return_S=FALSE,
           t_CS_n_init_steps = 30,
             t_CS = "mst")$cluster_membership$node
attr(results_per_strategy[["mst"]][["clq_MSTDunn_sT"]], "exec_time")= Sys.time()-t0

attr(results_per_strategy, "seed_used")= seeds[[trial_number]]



message("saving all results")
saveRDS(results_per_strategy,
	sprintf("leuk_clusters/n_genes=%d_trial=%d_clusters.rds", n_vc,
		trial_number)
	)
message("done")
