# Load required libraries
library(phyloseq)
library(ggClusterNet)
library(tidyverse)
library(WGCNA)
library(igraph)
library(Matrix)
library(vegan)
library(magrittr)
library(reticulate)
library(SpiecEasi)
library(ggstatsplot)
library(Hmisc)
library(dplyr)

# Read data files
otu <- read.csv("rarefied_otu_table.csv", row.names=1)
tax <- read.table("tax.txt", header=TRUE, row.names=1)
env <- read.table("env.txt", header=TRUE, row.names=1)
geo <- read.table("geo.txt", header=TRUE, row.names=1)

# Segment environmental factors into soil, climate, and plant categories
soil <- env[, 1:16]
climate <- env[, 17:20]
plant <- env[, 21]

# Read and prepare edge list for network construction
edgelist <- elist.mb
ig.spieceasi <- graph_from_edgelist(as.matrix(edgelist[,-3]), directed = FALSE)
ig.spieceasi <- set_edge_attr(ig.spieceasi, 'weight', index = E(ig.spieceasi), as.numeric(edgelist[,3]))

# Analyze ecological associations based on geographical distance
result_dl <- assigned_process(link_table_row=edgelist, OTUabd=otu, p=0.05, data=geo, cutoff=0, method="dl")

# Detect redundancy in environmental variables using hierarchical clustering
plot(varclus(as.matrix(env)))

# Evaluate ecological factors (soil) impact on network
results_soil <- lapply(names(soil), function(s) {
  assigned_process(link_table_row=edgelist, OTUabd=otu, p=0.05, data=soil[[s]], cutoff=0, method="ef")
})

# Aggregate soil links from results
soil_link <- unlist(lapply(results_soil, row.names))
unique_soil_link <- unique(soil_link)

# Detect redundancy in climate variables
plot(varclus(as.matrix(climate)))

# Evaluate ecological factors (climate) impact on network
results_climate <- lapply(names(climate), function(c) {
  assigned_process(link_table_row=edgelist, OTUabd=otu, p=0.05, data=climate[[c]], cutoff=0, method="ef")
})

# Aggregate climate links from results
climate_link <- unlist(lapply(results_climate, row.names))
unique_climate_link <- unique(climate_link)

# Evaluate ecological factors (plant) impact on network
result_ndvi <- assigned_process(link_table_row=edgelist, OTUabd=otu, p=0.05, data=plant['ndvi'], cutoff=0, method="ef")
plant_link <- row.names(result_ndvi)

# Identify and separate biotic links
total_link <- row.names(edgelist)
bi_link <- setdiff(total_link, c(unique_soil_link, unique_climate_link, plant_link))

# Save different ecological network edges to CSV files
write.csv(edgelist[unique_soil_link,], "soil_edge.csv")
write.csv(edgelist[unique_climate_link,], "climate_edge.csv")
write.csv(edgelist[plant_link,], "plant_edge.csv")
write.csv(edgelist[result_dl,], "dl_edge.csv")
write.csv(edgelist[bi_link,], "biotic_edge.csv")

# Construct complete and biotic graphs from edge lists
ig_total <- graph_from_edgelist(as.matrix(edgelist[,1:2]), directed = FALSE)
ig_total <- set_edge_attr(ig_total, 'weight', index = E(ig_total), as.numeric(edgelist[,3]))
ig_biotic <- graph_from_edgelist(as.matrix(edgelist[bi_link,1:2]), directed = FALSE)
ig_biotic <- set_edge_attr(ig_biotic, 'weight', index = E(ig_biotic), as.numeric(edgelist[bi_link,3]))

# Analyze subgraphs within the biotic network
otu_biotic <- otu[intersect(row.names(otu), vertex_attr(ig_biotic)[[1]]),]
sub_graph_ig_biotic <- lapply(names(otu_biotic), function(i) {
  sample_i <- otu_biotic[i]
  select_node <- rownames(sample_i)[which(sample_i != 0)]
  induced_subgraph(ig_biotic, select_node)
})

# Calculate network properties and save results
result_all <- sapply(sub_graph_ig_biotic, net_properties)
result_biotic <- do.call(cbind, result_all)
write.csv(result_biotic, "result_biotic.csv")

# Compute and save cohesion metrics using qcmi
cohesion_biotic <- qcmi(igraph= ig_biotic, OTU= otu_biotic, pers.cutoff=0)
re_qcmi1 <- data.frame(cohesion_biotic[4], cohesion_biotic[3])
write.csv(re_qcmi1, "re_qcmi1.csv")
