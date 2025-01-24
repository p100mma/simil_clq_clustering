library(igraph)
library(cliqueClusteR)
N_TRIALS=100
print("loading data")
X<- readRDS("gene_expr_data.rds")
print("oaded data")
bs_sizes<-c(200,nrow(X))
readRDS("v_bs_200.rds")-> bs_200
readRDS("v_bs_N.rds")-> bs_N
print("loaded bs indexes")
min_cl_sizes<- c(30, 400, 860) #WGCNA parameters
readRDS("vcl_clq_based_variants_to_test.rds")-> param_df #clique alg params
print("loaded param_df") 
n_methods<- length(min_cl_sizes) + nrow(param_df)
clq_ref<- rep(NA, length(param_df)) %>% as.list()
result_string<- "clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
for (clq_i in 1:nrow(param_df)){
	t_S<- as.character(param_df$t_S[[clq_i]])
	t_CS<- as.character(param_df$t_CS[[clq_i]])
	t_CS_relax_method<- as.character(param_df$t_CS_relax_method[[clq_i]])
	expansion_mode= as.character(param_df$expansion_mode[[clq_i]])
	fname<-sprintf(result_string, t_S,t_CS,expansion_mode,t_CS_relax_method)
	clq_ref[[clq_i]] <-uniqueSingletonLabels( readRDS(fname)$cluster_membership$node  )
	}
print("loaded_ref_clq")
ref_WG<-readRDS("WGCNA_labels_fulldata.rds")
print("load ref_WGCNA")
names_methods<- c(paste0(1:nrow(param_df),"clq"),
                  paste0(1:length(min_cl_sizes), "WGCNA"))
print(names_methods)
#2x100xn_methods numbers to retrieve
ARI_matrices<- list()
for (i_N in seq_along(bs_sizes)){

   print("working on size:")
   print(bs_sizes[[i_N]])
   ARI_matrix_i_N<- matrix(nrow=N_TRIALS,
			   ncol=n_methods)
   colnames(ARI_matrix_i_N)<- names_methods
for (B in 1:N_TRIALS) {
       for (comb_number in 1:nrow(param_df)) #clq_methods
	       {
	       t_S=as.character(param_df$t_S[[comb_number]])
	       t_CS=as.character(param_df$t_CS[[comb_number]])
	       expansion_mode=as.character(param_df$expansion_mode[[comb_number]])
	       t_CS_relax_method=as.character(param_df$t_CS_relax_method[[comb_number]])
	       result_string<- "B_%d/%d_clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
	       fname<- sprintf(result_string, bs_sizes[[i_N]], B, t_S, t_CS, expansion_mode, t_CS_relax_method)
	       readRDS(fname)-> result_cl
	       ARI_matrix_i_N[B,comb_number]= compare(
			uniqueSingletonLabels(result_cl$cluster_membership$node),
						      clq_ref[[comb_number]],
						      "adjusted.rand")
		}
	
   	wg_B<-readRDS(sprintf("B_%d/%d_WGCNA_labels.rds",bs_sizes[[i_N]], B ))
	for (wg_i in 1:length(ref_WG)){
	ARI_matrix_i_N[B, nrow(param_df) + wg_i ]= compare(ref_WG[[wg_i]], wg_B[[wg_i]], "adjusted.rand")		
		}
	}
   ARI_matrices[[i_N]]=ARI_matrix_i_N
}
saveRDS(ARI_matrices, "v2c_aggregated_stability_results.rds")
warnings()
      
      
      
