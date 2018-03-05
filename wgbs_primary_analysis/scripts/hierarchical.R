# creates hierarchical clustering of ALS most variable methylation sites
# fall 2017
# author <christa.caggiano@ucsf.edu>

rm(ls=list())

# load data
require(data.table)
sites = read.table("merged_most_var.txt")

# create distance matrix of the covariance
# will probably be extremely slow
distances = dist(cov(sites))

# find clusters based
clusters <- hclust(distances)

# save clusters as R object
saveRDS(clusters, "cluster.rdb")

# output a png of the clusters
png("clust.png")
plot(clusters)
dev.off()