
all_stuff<-readRDS("leuk_clusters_results_aggregated.rds")

library(ggplot2)

flatten_algo_summary<-function(algo_mat,thr_variant){
  
        algo<-colnames(algo_mat)
        colnames(algo_mat)<-NULL
        data.frame(algorithm= rep(algo,nrow(algo_mat)),
          score=as.vector(t(algo_mat)),
          thr=thr_variant)
}


do.call(rbind,lapply(seq_along(all_stuff), function(n_vc)
    { all_stuff[[n_vc]]-> vc_variant
        rbind(
        flatten_algo_summary(vc_variant$ARI$none, "none"),
        flatten_algo_summary(vc_variant$ARI$complexity, "complexity"),
        flatten_algo_summary(vc_variant$ARI$mst, "mst")
        )-> vc_df
        vc_df$v_cutoff<- n_vc
        vc_df
          }
  ))-> ARI_df

library(dplyr)

v_cutoff_map<- c(
  100,500,1000,
                      2000,3500, 5000,
                    7500, 10000, 12500,
                    15000, 17500, 20000
)

# replace 1:12 labels with actual gene numbers used
ARI_df$v_cutoff<- v_cutoff_map[ ARI_df$v_cutoff ]

cutoffs_left<- c(500, 1000,2000, #only those used in publication
                 5000,10000,20000)
nrow(ARI_df)

ARI_df<- ARI_df[ ARI_df$v_cutoff %in% cutoffs_left, ]
simil_based<- c("clq_modularity","clq_MSTDunn_sT",      
                "infomap","label_prop",       
                "leiden","affinity_prop",    
                "hcl_single","hcl_complete",     
                "hcl_average")


metric_based<- c("hcl_ward.D2",
                 "k_means","db_scan")

tab3<- c(metric_based, "affinity_prop")

summarise_ARI<- function(DF, subgroup, q=0.75) {
  DF[subgroup,] %>% group_by(v_cutoff,algorithm) %>%
  summarise(Q1=quantile(score,0.25),
            Qu=quantile(score,q)) %>% as.data.frame() -> df_long
  colnames(df_long)[[4]]<- paste0("Q", ifelse(q==0.75,3,q))
  reshape(df_long, direction="wide", timevar="algorithm",
          idvar="v_cutoff")
}

simil_based_pure= setdiff(simil_based, c("clq_modularity","clq_MSTDunn_sT"))
clq_based<-  c("clq_modularity","clq_MSTDunn_sT")
mask_Sall<- (ARI_df$algorithm %in% simil_based_pure)
mask_clq<- (ARI_df$algorithm %in% clq_based) & 
            (ARI_df$thr!="none") 

tab_simil_pub<- setdiff(simil_based, c("hcl_complete","hcl_average"))
pub_mask_pt1<- (ARI_df$algorithm %in% tab_simil_pub) &
                (ARI_df$thr=="complexity")
summarise_ARI(ARI_df, pub_mask_pt1) -> pubtab_s_pt1
summarise_ARI(ARI_df, (ARI_df$algorithm=="affinity_prop") &
                (ARI_df$thr=="mst") ) -> pubtab_s_pt2

colnames(pubtab_s_pt2)[2:3]= paste0(colnames(pubtab_s_pt2)[2:3],"_mst")
pubtab_s_pt2

pubtab_s<- cbind(pubtab_s_pt1, pubtab_s_pt2[,2:3])

pubtab_s

pub_mask_metric<-(ARI_df$algorithm %in% tab3 ) &
                  (ARI_df$thr=="none")

summarise_ARI(ARI_df, pub_mask_metric)-> pubtab_m

write.csv(pubtab_s,"similarity_based_methods_publication_table.csv")
write.csv(pubtab_m,"metric_based_methods_publication_table.csv")


format_tab2tex<- function(DataFrame,digits=2){
  DataFrame<-round(DataFrame,digits)
  res<- matrix(ncol=2*ncol(DataFrame),
               nrow=nrow(DataFrame))
  for (i in 1:ncol(DataFrame)){
    res[,2*(i-1)+1]= DataFrame[,i]
  }
  colnames(res)<- as.vector(rbind(colnames(DataFrame),"&"))
  
  res<-as.data.frame(res, check.names=FALSE)
  for (j in 1:ncol(res))
    if (j%%2==0)
      res[,j]="&"
  res[,ncol(res)]="\\\\"
  colnames(res)[[ncol(res)]]<-"\\\\"
  res
}

ARItab_format<- function(init_tab,q13sep=":"){
 format_tab2tex(init_tab)-> tex_format
 for (j in 1:ncol(tex_format)){
   if (j %%4==0) tex_format[,j]=q13sep
   if (j %%4==0) colnames(tex_format)[[j]]=""
   if (j %%4==3) colnames(tex_format)[[j]]= stringr::str_remove(
                                        colnames(tex_format)[[j]],
                                        "Q1.")
   if (j %%4==1) colnames(tex_format)[[j]]= stringr::str_remove(
                                        colnames(tex_format)[[j]],
                                        "Q3.+" 
                                        )
   
   
 }
  tex_format
}

ARItab_format(pubtab_s)-> tex_s_tab
ARItab_format(pubtab_m)-> tex_m_tab

name_map=list(v_cutoff="$p$",affinity_prop="\\verb|aff_prp|",
              clq_MSTDunn_sT="\\verb|clq_mdn|",
              clq_modularity="\\verb|clq_mod|",
              hcl_average="\\verb|hcl_avg|",
              hcl_single="\\verb|hcl_sng|",
              hcl_complete="\\verb|hcl_cmp|",
              label_prop="\\verb|lab_prp|",
              hcl_ward.D2="\\verb|hcl_wrd|",
              infomap="\\verb|infomap|",
              leiden="\\verb|leiden|",
              affinity_prop_mst="\\verb|aff_mst|",
              k_means="\\verb|k_means|",
              db_scan="\\verb|db_scan|")
map_tab_names<- function(tex_tab, map2use){
  cn<-colnames(tex_tab)
  new_cn<-cn
for (j in 1:ncol(tex_tab))
{
  if (cn[[j]] %in% names(map2use)) new_cn[[j]]= map2use[[cn[[j]]]] 
}
  return(new_cn)
}


colnames(tex_s_tab)<-map_tab_names(tex_s_tab,name_map)
colnames(tex_m_tab)<-map_tab_names(tex_m_tab,name_map)

write.csv(tex_s_tab,"s_compARI.csv")
#write.csv(ARItab_format(table_2),"2compARI.csv")
write.csv(tex_m_tab,"m_compARI.csv")


db_clq_df<- ARI_df[ (ARI_df$algorithm %in% c("db_scan") &
                       ARI_df$thr %in% c("none") ) |
                      ( ARI_df$algorithm %in% c("clq_modularity", 
                                                "clq_MSTDunn_sT")
                        & (ARI_df$thr=="complexity")
                      )
                    ,]

db_clq_df$v_cutoff<- db_clq_df$v_cutoff/1000

db_clq_df$algorithm[db_clq_df$algorithm=="clq_modularity"]="clq_mod(c)"
db_clq_df$algorithm[db_clq_df$algorithm=="clq_MSTDunn_sT"]="clq_mdn(c)"

pdf("ARI_clq_dbscan.pdf", width=4, height =5 )
ggplot(db_clq_df, aes(x=score, y= as.factor(v_cutoff),
                      fill=v_cutoff)) +
   geom_boxplot(colour="darkred")+
  facet_wrap(algorithm~., ncol=1) +
  labs(y="# genes (in thousands)", x="resulting ARI") +
  guides(fill="none")


dev.off()

ARI_df[mask_Sall,] %>% group_by(v_cutoff,algorithm,thr) %>%
  summarise(Q1=quantile(score,0.25),
            Q2=median(score),
            Q3=quantile(score,0.75)) %>% as.data.frame() -> iqr_df
iqr_df[,1:2] %>% unique() -> v_alg
v_alg_thr<-cbind(v_alg,best_thr=rep(NA,nrow(v_alg)),Q2=rep(NA,nrow(v_alg)),IQR=rep(NA,nrow(v_alg)))
for (i in 1:nrow(v_alg)){
  vat<- v_alg[i,]
 (iqr_df$v_cutoff==vat[[1]])&
 (iqr_df$algorithm==vat[[2]])-> sub_mask
  iqr_df[sub_mask,]->chunk
  which((chunk$Q1 > chunk$Q3[c(2,1,2)]) & (chunk$Q1 > chunk$Q3[c(3,3,1)]) )-> x
  if (length(x)) {
    (v_alg_thr$v_cutoff==chunk[x,]$v_cutoff) &
    (v_alg_thr$algorithm==chunk[x,]$algorithm) -> v_alg_mask
    stopifnot(sum(v_alg_mask)==1)
    v_alg_thr$best_thr[[which(v_alg_mask)]]= chunk[x,]$thr
    v_alg_thr$Q2[[which(v_alg_mask)]]= chunk[x,]$Q2
    v_alg_thr$IQR[[which(v_alg_mask)]]= chunk[x,]$Q3 - chunk[x,]$Q1
    
  }
}
print("unambigiously best thr per tested variant:")
v_alg_thr
ARI_df[mask_clq,]
ARI_df[mask_clq,] %>% group_by(v_cutoff,algorithm,thr) %>%
summarise(Q1=quantile(score,0.25),
          Q2=median(score),
          Q3=quantile(score,0.75)) %>% as.data.frame() -> iqr_df_2
iqr_df_2[,1:2] %>% unique() -> v_alg_2
v_alg_thr_2<-cbind(v_alg_2,best_thr=rep(NA,nrow(v_alg_2)),Q2=rep(NA,nrow(v_alg_2)),
                   IQR=rep(NA,nrow(v_alg_2)))
for (i in 1:nrow(v_alg_2)){
  vat<- v_alg_2[i,]
  (iqr_df_2$v_cutoff==vat[[1]])&
    (iqr_df_2$algorithm==vat[[2]])-> sub_mask
  iqr_df_2[sub_mask,]->chunk
  which((chunk$Q1 > chunk$Q3[c(2,1)]) )-> x
  if (length(x)) {
    (v_alg_thr_2$v_cutoff==chunk[x,]$v_cutoff) &
      (v_alg_thr_2$algorithm==chunk[x,]$algorithm) -> v_alg_mask
    stopifnot(sum(v_alg_mask)==1)
    v_alg_thr_2$best_thr[[which(v_alg_mask)]]= chunk[x,]$thr
    v_alg_thr_2$Q2[[which(v_alg_mask)]]= chunk[x,]$Q2
    v_alg_thr_2$IQR[[which(v_alg_mask)]]= chunk[x,]$Q3 - chunk[x,]$Q1
    
  }
}
print("unambigiously best thr per tested variant of simCliqNet:")
v_alg_thr_2

write.csv(v_alg_thr,"simil_based_best_thr_per_algVn_genes.csv")
write.csv(v_alg_thr_2,"clique_based_best_thr_per_algVn_genes.csv")
