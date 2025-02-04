library(igraph)
library(cliqueClusteR)
N_TRIALS=100
print("loading data")
X<- readRDS("gene_expr_data.rds")
print("loaded data")
bs_sizes<-c(200,nrow(X))
readRDS("v_bs_200.rds")-> bs_200
readRDS("v_bs_N.rds")-> bs_N
print("loaded bs indexes")
readRDS("vcl_clq_based_variants_to_test.rds")-> param_df #clique alg params
print("loaded param_df") 
method_sel<- c(2,3,4,6) #those with cores
n_methods<- length(method_sel)
clq_core<- rep(NA, n_methods) %>% as.list()
result_string<- "clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
for (j  in seq_along(method_sel)){ clq_i<-method_sel[[j]]
	t_S<- as.character(param_df$t_S[[clq_i]])
	t_CS<- as.character(param_df$t_CS[[clq_i]])
	t_CS_relax_method<- as.character(param_df$t_CS_relax_method[[clq_i]])
	expansion_mode= as.character(param_df$expansion_mode[[clq_i]])
	fname<-sprintf(result_string, t_S,t_CS,expansion_mode,t_CS_relax_method)
	res_i<- readRDS(fname)
	cl_mem<- res_i$cluster_membership$node
	is_core<-attr(res_i$cluster_membership,"core")$node
	core_colors<- cl_mem
	core_colors[ !is_core ] = 0
	clq_core[[j]] = core_colors
	}
print(lapply(clq_core, length))
print("loaded_ref_clq")
saveRDS(clq_core,"v2d_2346_cores.rds")
names_methods<- paste0(method_sel,"clq")
#2x100xn_methods numbers to retrieve
ARI_matrices<- list()
for (i_N in seq_along(bs_sizes)){
   print("working on size:")
   print(bs_sizes[[i_N]])
   ARI_matrix_i_N<- matrix(nrow=N_TRIALS,
			   ncol=n_methods)
   colnames(ARI_matrix_i_N)<- names_methods
for (B in 1:N_TRIALS) {
       for (j in seq_along(method_sel)) #clq_methods
	       { comb_number<- method_sel[[j]]
	       t_S=as.character(param_df$t_S[[comb_number]])
	       t_CS=as.character(param_df$t_CS[[comb_number]])
	       expansion_mode=as.character(param_df$expansion_mode[[comb_number]])
	       t_CS_relax_method=as.character(param_df$t_CS_relax_method[[comb_number]])
	       result_string<- "B_%d/%d_clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
	       fname<- sprintf(result_string, bs_sizes[[i_N]], B, t_S, t_CS, expansion_mode, t_CS_relax_method)
	       readRDS(fname)-> res_j
		cl_mem<- res_j$cluster_membership$node
		is_core<-attr(res_j$cluster_membership,"core")$node
		core_colors<- cl_mem
		core_colors[ !is_core ] = 0
		#print(length(is_core))
		#print(any(is.na(core_colors)))
		#print(any(is.na(clq_core[[j]])))
		#print(sum(is.na(clq_core[[j]])))
		#print(length(clq_core[[j]]))
		
	       ARI_matrix_i_N[B,j]= compare(
			zeroOutSmall(uniqueSingletonLabels(core_colors),100),
			zeroOutSmall(uniqueSingletonLabels(clq_core[[j]]),100),
						      "adjusted.rand")
		}
	
		}
   ARI_matrices[[i_N]]=ARI_matrix_i_N
}
saveRDS(ARI_matrices, "v2d_aggregated_stability_cores.rds")
warnings()
