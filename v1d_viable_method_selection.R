lv<-readRDS("devel_label_vecs.rds")
sm<-readRDS("devel_method_similarity.rds")
df<-readRDS("devel_method_runtime_df.rds")
library(magrittr)

sm_order<-order(colSums(sm))
smv<- sm[sm_order,sm_order] #visualization purposes
lapply(lv, function(x) max(table(x)))[sm_order] %>% 
  unlist() -> max_cl
#see v1_visual_summary.Rmd for the plots
#plot(max_cl, 
#        xlab="method (heatmap order)",
#        ylab="size of biggest cluster", 
#        main="method vs max.cl size (dataset of 8676 genes)",
#        ylim=c(0,9000))
#abline(h=8676)

colnames(smv)[(max_cl<8000) & (max_cl>2) ]
viable_methods<-cbind(max_cl_size=unname(max_cl), df[sm_order,])
viable_methods<- viable_methods[ (viable_methods$max_cl_size>2) & 
                                  (viable_methods$max_cl_size<8000),
]
viable_methods

smv[(max_cl<8000) & (max_cl>2),
    (max_cl<8000) & (max_cl>2)]-> smv_viable
short_names<-stringr::str_replace(colnames(smv_viable),"t_S.","")
short_names<-stringr::str_replace(short_names,"t_CS.","")
short_names<-stringr::str_replace(short_names,"clq.","")
short_names<-stringr::str_replace(short_names,"alg.","")
short_names<-stringr::str_replace_all(short_names," ",".")
short_names<-stringr::str_replace_all(short_names,"\\.\\.",".")
short_names<-stringr::str_replace(short_names,"\\.","")
short_names<-stringr::str_replace(short_names,"complexity","c")
short_names<-stringr::str_replace(short_names,"modularity","mo")
short_names<-stringr::str_replace(short_names,"mst","ms")
short_names<-stringr::str_replace(short_names,"greedyCliqueJoin","g")
short_names<-stringr::str_replace(short_names,"components","c")
short_names<-stringr::str_replace(short_names,"max","m")
short_names<-stringr::str_replace(short_names,"basic","b")
short_names<-stringr::str_replace(short_names,"custom","c")
short_names<-stringr::str_replace(short_names,"average","a")
short_names<-stringr::str_replace(short_names,"outW","ow")
short_names<-stringr::str_replace(short_names,"outN","on")
short_names
colnames(smv_viable)<-rownames(smv_viable)<-short_names
#see v1_visual_summary.Rmd for the plots
#library(gplots)
#heatmap.2(smv_viable, cellnote = round(smv_viable,2))

#clusters of similar variants according to heatmap
methods_by_similarity= list(
              CS_components=c(
              "c.mo.on.c","c.mo.ow.c",
              "c.mo.m.c","c.mo.a.c",
              "c.mo.b.c"),
              greedyJoin_out=c("c.mo.on.g",
                               "c.mo.ow.g"),
              greedyJoin_in=c("c.mo.b.g",
                              "c.mo.m.g",
                              "c.mo.a.g")
                      )

#find most representative method
in_deg<-lapply(methods_by_similarity, function(cl)
                colMeans(smv_viable[cl,cl]))

inv_o_deg<- lapply(methods_by_similarity, function(cl)
                colMeans( 1- smv_viable[ !(rownames(smv_viable)%in% cl),cl] ))

lapply(seq_along(inv_o_deg), function(i_cl)
          sort(in_deg[[i_cl]]+inv_o_deg[[i_cl]], decr=TRUE)
        )


#then, from greedy out based, we have chosen outN because it was more
#different to other methods from the two

#from the components cluster, we took 2 methods because 
#we have subcluster of out and in based here.
# from here we chose "basic" (best representative)
# and from out subcluster we took outW because outN was taken

#from greedy in based, we tested "max" based one

representative_variants<-c("c.mo.c.c", "c.mo.c.g",
                           "c.mo.on.g",
                           "c.mo.m.g",
                           "c.mo.ow.c",
                           "c.mo.b.g")
viable_subset_to_test_IDX<- match(representative_variants,
                                  colnames(smv_viable))

 to_test_df<- viable_methods[viable_subset_to_test_IDX,2:5]
saveRDS(to_test_df,"vcl_clq_based_variants_to_test.rds")
clq_ref<- rep(NA, length(to_test_df)) %>% as.list()
result_string<- "clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
for (clq_i in 1:nrow(to_test_df)){
        t_S<- as.character(to_test_df$t_S[[clq_i]])
        t_CS<- as.character(to_test_df$t_CS[[clq_i]])
        t_CS_relax_method<- as.character(to_test_df$t_CS_relax_method[[clq_i]])
        expansion_mode= as.character(to_test_df$expansion_mode[[clq_i]])
        fname<-sprintf(result_string, t_S,t_CS,expansion_mode,t_CS_relax_method)
        clq_ref[[clq_i]] <- readRDS(fname)$cluster_membership$node
        }
saveRDS(clq_ref,"clq_reference_modules.rds")
