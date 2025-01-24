readRDS("v2c_aggregated_stability_results.rds")-> ARI_matrices
#readRDS("v2c_aggregated_stability_results_min300.rds")-> ARI_matrices

library(reshape2)
bs_sizes=c(200,1394)
ARI_dfs<- cbind(lapply(ARI_matrices, melt))
for (i in seq_along(ARI_dfs)) ARI_dfs[[i]]$N<- bs_sizes[[i]]
do.call(rbind,ARI_dfs)->ARI_df
colnames(ARI_df)<- c("trial","algorithm","ARI","BS_size")

readRDS("v2d_2346_cores.rds")->cores_2346
readRDS("v2d_aggregated_stability_cores_min100.rds")-> ARI_cores_matrices

ARI_dfs<- cbind(lapply(ARI_cores_matrices, melt))
for (i in seq_along(ARI_dfs)) ARI_dfs[[i]]$N<- bs_sizes[[i]]
do.call(rbind,ARI_dfs)->ARI_df_core
colnames(ARI_df_core)<- c("trial","algorithm","ARI","BS_size")
ARI_df_core$algorithm<- paste0(ARI_df_core$algorithm,"_core")
ARI_all<-rbind(ARI_df,ARI_df_core)
library(ggplot2)
pdf("v2e_stability.pdf", width=8,height=6)
ggplot(ARI_all[!(ARI_all$algorithm %in% c("2WGCNA","3WGCNA")),],
       aes(x=as.factor(BS_size), y=ARI, fill=algorithm)) +
  geom_boxplot()
dev.off()







readRDS("WGCNA_labels_fulldata.rds")-> fulldata_WGCNA
readRDS("clq_reference_modules.rds")-> fulldata_clq
ALL<-c(fulldata_clq,fulldata_WGCNA[1])
length(ALL)
names(ALL)<-c(paste0(1:6,"clq"),"WGCNA")
all2all<-matrix(nrow=length(ALL),ncol=length(ALL))
colnames(all2all)<-rownames(all2all)<-names(ALL)
library(igraph)
for (i in 1:length(ALL))
  for (j in 1:length(ALL))
    all2all[i,j]= compare(ALL[[i]],ALL[[j]],"adjusted.rand")
heatmap(all2all)->h
h$rowInd-> sim_order
melt(all2all[sim_order,sim_order])-> sim_df
colnames(sim_df)[[3]]="ARI"
pdf("v2e_pairwise_sim.pdf", height=6)
ggplot(sim_df, aes(Var1,Var2, fill=ARI)) +
  geom_tile() + scale_fill_gradient(low="blue",high="red")+
  geom_text(aes(label = round(ARI,2)),size=6, color="white")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=12),
        axis.text.y=element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11)
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))
dev.off()
  #scale_y_discrete(guide = guide_axis(n.dodge = 2))
X<-readRDS("gene_expr_data.rds")
#source("cor_matmul.R")
S<- similarity_matrix(fastPearsonData(X))
lapply(fulldata_clq, function(P)
cliqueClusteR:::`%modularity%`(S,P)) %>% unlist -> MOD_S_clq
lapply(fulldata_WGCNA, function(P)
cliqueClusteR:::`%modularity%`(S,P)) %>% unlist -> MOD_S_wgc
c(MOD_S_clq,MOD_S_wgc)-> MOD_S
MOD_S<-MOD_S[1:7]
rm(S)
#A<- corfast(X)^5
#A[A<0]=0
#lapply(fulldata_clq, function(P)
#cliqueClusteR:::`%modularity%`(A,P)) %>% unlist -> MOD_A_clq
#lapply(fulldata_WGCNA, function(P)
#  cliqueClusteR:::`%modularity%`(A,P)) %>% unlist -> MOD_A_wgc
#c(MOD_S_clq,MOD_S_wgc) -> MOD_A
#rm(A)
#round(cbind(MOD_S,MOD_A),3) %>% as.data.frame()-> mod_df
mod_df<-round(MOD_S,3) %>% as.data.frame()
cbind(algorithm=names(ALL),mod_df)->mod_df
mod_df_viz<-mod_df
mod_df_viz$x_pos=rep(1,nrow(mod_df_viz))
colnames(mod_df_viz)[[2]]="modularity"
mod_df
h_order<-c(1,4,5,7,2,6,3)
pdf("v2e_modularity.pdf",height=6,width=3.3)
ggplot(mod_df_viz, aes(x_pos,reorder(algorithm,h_order), fill=modularity)) +
  geom_tile() + scale_fill_gradient(low="lightblue",high="yellow")+
  geom_text(aes(label = modularity),size=6)+
  labs(x="modularity")+
  theme(#axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_text(size=20, face="bold"),
    axis.text=element_text(size=12),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = "left"
  )+ 
  scale_x_discrete(guide = guide_axis(n.dodge = 2),expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))
dev.off()
