library(cliqueClusteR)
library(FCPS)
library(magrittr)

sets2test<- c("leukemia","Atom", "Chainlink", "Tetra", "Target", "TwoDiamonds")

expand.grid(t_S=c("complexity","mst"),
	    t_CS_relax_method=c("components","greedyCliqueJoin"),
	    t_CS=c("modularity","mst")) -> fcps_variants
print(fcps_variants)
print(class(fcps_variants$t_S[[1]]))


results_per_set<- rep(NA, nrow(fcps_variants)) %>% as.list()
names(results_per_set)<- sets2test
for (set in sets2test) {
if (set=="leukemia"){
	data(leukemia)
	data(leukemia_clusters)
	scpcss_input<- leukemia
	ref<-leukemia_clusters} else{
	get(set)$Data -> scpcss_input	
	get(set)$Cls -> ref
	}
	method_results<- list()
	for (i in 1:nrow(fcps_variants)){
		t_S<- fcps_variants$t_S[[i]] %>% as.character()
		t_CS<- fcps_variants$t_CS[[i]] %>% as.character()
		t_CS_relax_method<- fcps_variants$t_CS_relax_method[[i]] %>% as.character()
		X_simil_fun<-NULL
		if (set!="leukemia") X_simil_fun="euclidean"
		run_scpcss(scpcss_input, do_signif = FALSE, 
			   X_simil_fun = X_simil_fun,
			   t_S = t_S, t_S_init_steps = 80,
           		    t_CS_relax_method = t_CS_relax_method, return_CS =FALSE, return_S=FALSE,
           		    t_CS_n_init_steps = 30,
             		    t_CS = t_CS)-> fcps_res_clq_i
	method_results[[i]]<-fcps_res_clq_i	
	}
	method_results-> results_per_set[[set]]
		}
saveRDS(results_per_set,"cliquesOnFCPS.rds")


