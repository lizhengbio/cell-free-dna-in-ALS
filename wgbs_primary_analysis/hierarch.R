
o1 = read.table("merged_most_var.txt")

o2 = dist(cov(o1))
clusters <- hclust(o2)

saveRDS(clusters, "cluster.rdb")

png("clust.png")
plot(clusters)
dev.off()