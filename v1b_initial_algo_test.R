args=commandArgs(trailingOnly=TRUE)
as.integer(args[[1]])-> comb_number #1 to 48 over all combinations below

expand.grid(
t_S=c("complexity", "mst"),
t_CS=c("modularity", "mst"),
expansion_mode=c("basic", "average", "max", "outN", "outW", "custom" ), #custom - coloring of complement graph without caring about weights(unW variant)
t_CS_relax_method=c("components", "greedyCliqueJoin"))-> param_df


t_S=as.character(param_df$t_S[[comb_number]])
t_CS=as.character(param_df$t_CS[[comb_number]])
expansion_mode=as.character(param_df$expansion_mode[[comb_number]])
t_CS_relax_method=as.character(param_df$t_CS_relax_method[[comb_number]])

library(cliqueClusteR)
S<- readRDS("BRCA_similarity_graph.rds")
print("loaded graph")

if (expansion_mode!="custom") {
	t0= Sys.time()
	run_scpcss(S, t_S=t_S, t_CS=t_CS, expansion_mode=expansion_mode,
		   t_S_keep_init_partitions=TRUE,
		   t_CS_relax_method=t_CS_relax_method,
		   t_CS_keep_init_partitions= TRUE)-> clustering_result
	runtime= Sys.time() - t0
} else {

    source("igraph_color_complement.R")
     t0 =Sys.time()
    if (t_S!="mst"){
	    opt_res<- thr_optimizer(S=S, keep_init_partitions=TRUE, 
		     custom_partitioner=igraph_color_complement) 
	    clq_P<- opt_res$maximizer_partition
	    thr_S<- opt_res$thr
    }else{
	thr_S<- critical_mst_thr(S)	    
	clq_P<- igraph_color_complement( S %thr% thr_S)$membership
    }
    rel_res<- relax_cliques_optimally(partition=clq_P ,S=S, t_S= thr_S, relax_method= t_CS_relax_method,
				      cs_thr_objective=t_CS, keep_init_partitions=TRUE,
					do_signif=FALSE, dX.Y="hausdorff")
    clustering_result<- list( clique_membership= clq_P, t_S= thr_S, t_CS= rel_res$frac,
			      cluster_membership= rel_res$membership)
    if (t_S!="mst"){
	clustering_result$t_S_objective= opt_res$objective
	clustering_result$t_S_init_partitions= opt_res$init_partitions
	clustering_result$t_S_init_search_points= opt_res$init_search_points
	clustering_result$t_S_init_scores= opt_res$init_scores
	}

    clustering_result$t_CS_objective= rel_res$objective
    clustering_result$t_CS_init_partitions= rel_res$init_partitions
    clustering_result$t_CS_init_search_points= rel_res$init_search_points
    clustering_result$t_CS_init_scores= rel_res$init_scores
    runtime= Sys.time() - t0
}
clustering_result$runtime= runtime
print("done clustering")

result_string<- "clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
fname<- sprintf(result_string, t_S, t_CS, expansion_mode, t_CS_relax_method)
saveRDS(clustering_result, fname)
print("saved result")
