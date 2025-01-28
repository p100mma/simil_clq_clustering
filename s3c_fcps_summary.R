readRDS("cliquesOnFCPS.rds")-> clq_results
clq_results<-clq_results[1:6]
names(clq_results)-> datasets
length(clq_results[[1]])
library(FCPS)
library(igraph)
library(magrittr)
library(cliqueClusteR)
data("leukemia_clusters")

expand.grid(t_S=c("complexity","mst"),
            t_CS_relax_method=c("components","greedyCliqueJoin"),
            t_CS=c("modularity","mst")) -> fcps_variants
datasets=c("leukemia","Atom","Chainlink","Tetra","Target","TwoDiamonds")
do.call(cbind,lapply(datasets, function(nm) 
        lapply(clq_results[[nm]],
               function(variant)
               {
                 ref<- if (nm!="leukemia") get(nm)$Cls else leukemia_clusters
                 compare(variant$cluster_membership$node, 
                         ref, "adjusted.rand")
               }
                 ) %>%
         unlist()
       )) %>% round(2)-> FCPS_scores
datasets[[1]]<-"Leuk_6"
colnames(FCPS_scores)<- datasets
cbind(fcps_variants,FCPS_scores)-> clq_scores
clq_num<- clq_scores[,datasets]
apply(clq_scores[,1:3],1, paste0,collapse=".")-> concat_names

readRDS("DBscan_FCPS_scores.rds")->dbscan_scores
rbind(
cbind(concat_names, clq_num),
c("db_scan",NA,dbscan_scores%>% round(2)))-> tab_publication

#tab_publication
alg_order<- c(6,5,2,1,9,3,4,8,7)
tab_publication<-tab_publication[rev(alg_order),]
new_names<- c("clq_cgMST",
              "clq_mgMST",
              "clq_mgMOD",
              "clq_cgMOD",
              "db_scan",
              "clq_ccMOD",
              "clq_mcMOD",
              "clq_ccMST",
              "clq_mcMST")
tab_publication$concat_names<- new_names

library(reshape2)
long_tab<- as.matrix(tab_publication[,-1])
rownames(long_tab)<- tab_publication$concat_names      
melt(long_tab)-> hmap_df
colnames(hmap_df)<-c("algorithm","dataset","ARI")
hmap_df$ARI<-as.numeric(hmap_df$ARI)
library(ggplot2)
pdf("fcps_results.pdf")
ggplot(hmap_df, aes(dataset,algorithm,fill=ARI))+ geom_tile() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=13),
        legend.text = element_text(size = 13.5),
        legend.title = element_text(size = 15)
        )+
  geom_text(aes(label = ARI),size=7) +
  scale_fill_gradient(low="lightblue",high="yellow") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
dev.off()

#plots of 3d point clouds: building blocks

OneDim3dRot<- function(dim.nr, angle)
{
  mat<- rbind(c(cos(angle), sin(angle)),
              c(sin(angle), cos(angle)))
  mat3D<- matrix(rep(0,3*3),ncol=3)
  diag(mat3D)<-1
  mat3D[-dim.nr,-dim.nr]<- mat
  mat3D
}
rad<- function(angle) angle*pi/180
General3dRot<- function(angles_xyz)
{
  Reduce(`%*%`,lapply(1:3, function(x) OneDim3dRot(x, angles_xyz[[x]])))
}


#project coordinates on 2D and rotate if matrix has 3 columns otherwise return matrix without changes
project2D<-function(XYZ, degs=c(10,10,0)) {
  if (ncol(XYZ)==3) XYZ %*% General3dRot(lapply(degs, rad )) %*% rbind(c(1,0),c(0,1),c(0,0)) else XYZ }

pointClouds<- datasets[-c(1:2)]
plots<-list()
do.call(rbind,lapply(pointClouds, function(dat){
 XY<- project2D(get(dat)$Data)
 ref<- get(dat)$Cls
 cbind(as.data.frame(XY),ref=ref)-> points_df
 colnames(points_df)[1:2]=c("X","Y")
 points_df$data=dat
 points_df
}
))-> long_df
pdf("FCPS_pointClouds.pdf")
ggplot(long_df, aes(X,Y,col=as.factor(ref)))+geom_point(size=2.5)+
  facet_wrap(data~.,scales = "free", )+
  guides(col="none")+ theme_classic()+
  theme( strip.text =element_text(size=14, face='bold'))
dev.off()
