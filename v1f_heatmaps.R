args<-commandArgs(trailingOnly=TRUE)
base_nr<- as.integer(args[[1]])
hmap_size<- as.numeric(args[[2]])
library(cliqueClusteR)
clq_methods_df<-readRDS("vcl_clq_based_variants_to_test.rds")
clq_methods<- list()
clq_methods_clq<- list()
result_string<- "clusters_t_S=%s;t_CS=%s;mode=%s;join=%s.rds"
for (i in 1:nrow(clq_methods_df)) {
	readRDS(
	sprintf(result_string, 
	clq_methods_df$t_S[[i]], clq_methods_df$t_CS[[i]], clq_methods_df$expansion_mode[[i]], clq_methods_df$t_CS_relax_method[[i]])
		)-> i_file
i_file$cluster_membership$node -> clq_methods[[i]]
if (!(i %in% c(1,5)))
-attr(i_file$cluster_membership,"core")$node + 1-> clq_methods_clq[[i]]
}
readRDS("WGCNA_labels_fulldata.rds")-> WGCNA_

all_methods<- c( clq_methods, WGCNA_)
core_method<- c( clq_methods_clq, WGCNA_)


base_method<- all_methods[[base_nr]]
base_method_cores<- if (base_nr %in% c(1,5)) all_methods[[base_nr]] else core_method[[base_nr]]

set.seed(123)
#print(table(base_method))
fastTable(base_method)$count-> base_sizes
fastTable(base_method)$value-> base_lab


bs_order<- order(-base_sizes)

base_sizes<- base_sizes [ bs_order ][1:15]
base_lab  <- base_lab [ bs_order ][1:15]
print(base_sizes)
cl_distr<-  (base_sizes)^(1/3)
print(cl_distr)
cl_distr<- cl_distr/sum(cl_distr)
print(cl_distr)

potential_gene_set<- base_method[ base_method %in% base_lab ]

kept_genes_labels<- sample(base_lab, hmap_size, 
			  prob=cl_distr, replace=TRUE)

kept_gene_idx<- kept_genes_labels
kept_gene_idx[]<-NA

for (lab in unique(kept_genes_labels))
	kept_gene_idx[kept_genes_labels==lab ] = sample(
	which(base_method==lab), sum(kept_genes_labels==lab),
		replace=FALSE)

stopifnot(all(!is.na(kept_gene_idx)))








cluster_heatmap<-function(matr,labels_1,
			  row_order, main,fname,
			  col_order=row_order,
			  labels_2=labels_1) {
matr<- matr[row_order,col_order]
row_labels<- labels_1[row_order]
col_labels<- labels_2[col_order]

r_clrs<- colors()[row_labels+4]
c_clrs<- colors()[col_labels+4]

rownames(matr)<-colnames(matr)<-NULL
pdf(fname)
heatmap(matr, Colv=NA, Rowv=NA, scale="none", revC=FALSE,
	ColSideColors=c_clrs,
	RowSideColors=r_clrs,
	main=main)
dev.off()
}

X<- readRDS("gene_expr_data.rds")
X<- X[, kept_gene_idx]
library(matrixStats)
source("cor_matmul.R")
S<-corfast(X)^2

base_method<- base_method[ kept_gene_idx ]
base_method_cores<- base_method_cores[ kept_gene_idx ]



order(base_method, base_method_cores)-> r_ord
cluster_heatmap(sqrt(S),base_method, r_ord, "CliqueSimNet","4clq_S.pdf",
		labels_2=base_method_cores)
w=1
wg<- WGCNA_[[w]][ kept_gene_idx ]
order(wg)->r_ord_w
cluster_heatmap(sqrt(S),wg, r_ord_w,"WGCNA",sprintf("%d__WGCNAS.pdf",w),
		)
