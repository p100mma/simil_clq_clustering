
X<- readRDS("gene_expr_data.rds")
set.seed(32424)
N<- nrow(X)
N_TRIALS=100
bs_list_N<- lapply(1:N_TRIALS,
		 function(i) sample(1:N, N,replace=TRUE))
bs_list_200<- lapply(1:N_TRIALS,
		 function(i) sample(1:N, 200,replace=TRUE))


saveRDS(bs_list_N, "v_bs_N.rds")
saveRDS(bs_list_200, "v_bs_200.rds")
