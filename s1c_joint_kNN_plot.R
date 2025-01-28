
n_vc= 1:12

curves<- lapply(n_vc, function(n)

readRDS(
	sprintf("leuk_DB_scan_kNN_curve/n_genes=%d_kNN_d_curve.rds",n))
	)

range_y<- range(unlist(curves))

pdf("leuk_DB_scan_kNN_curve/joint_plot.pdf")
plot(curves[[1]], main=sprintf("kNNcurves, all variants"),
	ylim=range_y)
for (n in 2:12)
points(curves[[n]])
dev.off()
