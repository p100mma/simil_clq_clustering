
library(igraph)

expand.grid(
   t_S=c("complexity", "mst"),
   t_CS=c("modularity", "mst"),
   expansion_mode=c("basic", "average", "max", "outN", "outW", "custom" ),
   t_CS_relax_method=c("components", "greedyCliqueJoin"))-> param_df

IDers<- vector()
results<- list()
param_df$completed<- rep(NA, nrow(param_df))
param_df$exec_time<- rep(NA, nrow(param_df))
for (i in 1:48){
t_S=as.character(param_df$t_S[[i]])
t_CS=as.character(param_df$t_CS[[i]])
expansion_mode=as.character(param_df$expansion_mode[[i]])
t_CS_relax_method=as.character(param_df$t_CS_relax_method[[i]])
IDers[[i]]<- paste("t_S.",t_S,"t_CS.",t_CS,"clq.",expansion_mode,"alg.",t_CS_relax_method)

result_string<- "clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
fname<- sprintf(result_string, t_S, t_CS, expansion_mode, t_CS_relax_method)
param_df$completed[[i]]= file.exists(fname)
if (param_df$completed[[i]]) {
res<-readRDS(fname)
results[[i]]<- res
param_df$exec_time[[i]]= res$runtime} else { results[[i]]<-NA}
}
saveRDS(param_df, "devel_method_runtime_df.rds")

ARIs<- matrix(nrow=48,ncol=48)
for (i in 1:48) for (j in 1:48) if (param_df$completed[[i]] && param_df$completed[[j]]) {
					ARIs[i,j]= compare(results[[i]]$cluster_membership$node, 
						 	   results[[j]]$cluster_membership$node, 
							   "adjusted.rand")
						}
lapply(results, function(RES) RES$cluster_membership$node) -> label_vectors
names(label_vectors)<- IDers
colnames(ARIs)<-rownames(ARIs)<- IDers
saveRDS(ARIs,"devel_method_similarity.rds")
saveRDS(label_vectors,"devel_label_vecs.rds")

