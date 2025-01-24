args=commandArgs(trailingOnly=TRUE)
as.integer(args[[1]])-> B


X<- readRDS("gene_expr_data.rds")
#X<-X[,1:900]
bs_sizes<-c(200,nrow(X))
stopifnot(B %in% 1:100)
readRDS("v_bs_200.rds")-> bs_200
readRDS("v_bs_N.rds")-> bs_N
b_IDX<- list( bs_200[[B]], bs_N[[B]])




readRDS("vcl_clq_based_variants_to_test.rds")-> param_df

library(cliqueClusteR)
library(matrixStats)
library(WGCNA)
library(dynamicTreeCut)
library(fastcluster)
source("cor_matmul.R")
source("igraph_color_complement.R")


clres_N_clq<-list()

for (i_N in seq_along(bs_sizes)){

   print("working on size:")
   print(bs_sizes[[i_N]])
   IDX_iB<- b_IDX[[i_N]]
   X_iB<- X[ IDX_iB,]

   S<- similarity_matrix( fastPearsonData(X_iB) )
   print("clique based:")
for (comb_number in 1:nrow(param_df))
	{
	t_S=as.character(param_df$t_S[[comb_number]])
	t_CS=as.character(param_df$t_CS[[comb_number]])
	expansion_mode=as.character(param_df$expansion_mode[[comb_number]])
	t_CS_relax_method=as.character(param_df$t_CS_relax_method[[comb_number]])
	print(param_df[comb_number,])
	if (expansion_mode!="custom") {
		t0= Sys.time()
		run_scpcss(S, t_S=t_S, t_CS=t_CS, expansion_mode=expansion_mode,
			   t_S_keep_init_partitions=TRUE,
			   t_CS_relax_method=t_CS_relax_method,
			   t_CS_keep_init_partitions= TRUE)-> clustering_result
		runtime= Sys.time() - t0
	} else {

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
	labs<- clustering_result$cluster_membership$node

	result_string<- "B_%d/%d_clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
	fname<- sprintf(result_string, bs_sizes[[i_N]], B, t_S, t_CS, expansion_mode, t_CS_relax_method)
	saveRDS(clustering_result, fname)
	}
   rm(S)
   print("WGCNA")
   A<-abs(corfast(X_iB))^5
   A[A>1]=1; A[A<0]=0; #roundoff errors, adjust so WGCNA does not complain
   TOMdist(A)->dTOM; rm(A)
   hclust(as.dist(dTOM), "average")-> dendro
   min_cl_sizes<- c(30, 400, 860)
   labels_per_ms<-list()
   for (i in seq_along(min_cl_sizes))
   labels_per_ms[[i]]<-cutreeDynamic(dendro=dendro, 
	      method="hybrid", distM= dTOM,
	      deepSplit=2,minClusterSize=min_cl_sizes[[i]])
   labels_per_ms<-lapply(labels_per_ms, uniqueSingletonLabels)
   rm(dTOM)
   saveRDS(labels_per_ms,sprintf("B_%d/%d_WGCNA_labels.rds",bs_sizes[[i_N]], B ))
}
      
      
      
